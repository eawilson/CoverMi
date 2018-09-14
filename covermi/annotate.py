from __future__ import print_function, absolute_import, division

import requests, csv, json, subprocess, tempfile, os, pdb
from collections import defaultdict

try: # python2
    from cStringIO import StringIO
except ImportError:
    try:
        from StringIO import StringIO
    except ImportError: # python3
        from io import StringIO
        basestring = str

from .include import *
from .gr import load_targets, load_principal
from .cache import keyvalcache


def ncbi_esearch(db, keyvals):
    server = "https://eutils.ncbi.nlm.nih.gov"
    ext = "/entrez/eutils/esearch.fcgi"
    terms = []
    for key, vals in keyvals:
        if isinstance(vals, basestring):
            vals = (vals,)
        terms += ["("+" OR ".join(["{}[{}]".format(val, key) for val in vals])+")"]
    params = {"db": db, "retmode": "json", "term": " AND ".join(terms)}
    r = requests.get(server+ext, params=params)
    if not r.ok: 
        r.raise_for_status()
    return r.json()["esearchresult"]["idlist"]


def ncbi_esummary(db, items):
    server = "https://eutils.ncbi.nlm.nih.gov"
    ext = "/entrez/eutils/esummary.fcgi"
    params = {"db": db, "retmode": "json", "id": ",".join(items)}
    r = requests.post(server+ext, params=params)
    if not r.ok:
        r.raise_for_status()
    if "result" not in r.json():
        pdb.set_trace()
    return list(r.json()["result"].items())


@keyvalcache
def _ncbi(refids):
    limit = 450
    for index in range(0, len(refids), limit):
        for key, snp in ncbi_esummary("snp", refids[index:index+limit]):
            if key != "uids":
                yield(key, snp)


def ncbi(annotations, **kwargs):
    refids = defaultdict(list)
    for annotation in annotations:
        for ref in annotation.get("dbsnp", ()):
            refids[ref.strip("rs")] += [annotation]
    if refids:
        for refid, snp in _ncbi(refids.keys(), **kwargs):
            if "clinical_significance" in snp:
                for annotation in refids[refid]:
                    annotation["clinical_significance"] = snp["clinical_significance"] 


@keyvalcache
def _uniprotid(genes):
    results = requests.get("http://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+[9606]%22&sort=score&columns=id,genes&format=tab")
    for row in csv.DictReader(StringIO(results.text), delimiter="\t"):
        for gene in row["Gene names"].split():
            if gene:
                yield (gene, row["Entry"])


def uniprotid(genes, **kwargs):
    for gene, uniprotid in _uniprotid(["TP53"], **kwargs):
        pass
    for keyval in _uniprotid(genes, **kwargs):
        yield keyval


@keyvalcache
def _uniprottext(uniprotids, **kwargs):
    for uniprotid in uniprotids:    
        results = requests.get("http://www.uniprot.org/uniprot/{}.txt".format(uniprotid))
        yield (uniprotid, results.text)


def uniprottext(genes, **kwargs):
    for gene, identifier in _uniprotid(genes, **kwargs):
        for identifier, text in _uniprottext([identifier], **kwargs):
            yield (gene, text)






@keyvalcache
def _ncbigene(genes):
    for offset in range(0, len(genes), 10):
        subset = set(genes[offset:offset+10])
        geneids = ncbi_esearch("Gene", [("Gene Name", subset), ("Organism", "human")])
        if geneids:
            for geneid, response in ncbi_esummary("Gene", geneids):
                if geneid != "uids":
                    if response["name"] in subset:
                        yield(response["name"], response)


def ncbigene(genes, **kwargs):
    for gene, response in _ncbigene(genes, **kwargs):
        print(list(response.items()))


IMPACTVAL = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}

IMPACT = {
"transcript_ablation": "HIGH", 	
"splice_acceptor_variant": "HIGH", 	
"splice_donor_variant": "HIGH", 	
"stop_gained": "HIGH", 	
"frameshift_variant": "HIGH", 	
"stop_lost": "HIGH", 	
"start_lost": "HIGH", 	
"transcript_amplification": "HIGH", 	
"inframe_insertion": "MODERATE", 	
"inframe_deletion": "MODERATE", 	
"missense_variant": "MODERATE", 	
"protein_altering_variant": "MODERATE", 	
"splice_region_variant": "LOW", 	
"incomplete_terminal_codon_variant": "LOW", 	
"stop_retained_variant": "LOW", 
"synonymous_variant": "LOW", 
"coding_sequence_variant": "MODIFIER", 	
"mature_miRNA_variant": "MODIFIER", 	
"5_prime_UTR_variant": "MODIFIER", 	
"3_prime_UTR_variant": "MODIFIER", 
"non_coding_transcript_exon_variant": "MODIFIER", 	
"intron_variant": "MODIFIER", 
"NMD_transcript_variant": "MODIFIER", 	
"non_coding_transcript_variant": "MODIFIER", 	
"upstream_gene_variant": "MODIFIER", 	
"downstream_gene_variant": "MODIFIER", 	
"TFBS_ablation": "MODIFIER", 	
"TFBS_amplification": "MODIFIER", 	
"TF_binding_site_variant": "MODIFIER", 	
"regulatory_region_ablation": "MODIFIER", 	
"regulatory_region_amplification": "MODIFIER", 	
"feature_elongation": "MODIFIER", 	
"regulatory_region_variant": "MODIFIER", 	
"feature_truncation": "MODIFIER", 	
"intergenic_variant": "MODIFIER", 	
}

    
#examplefilters = filters == 'PASS'
#                 vaf > 0.2
#                 maf < 0.05
#                 impact in 'MODERATE/HIGH' or (gene_symbol == 'FTL' and '5_prime_UTR_variant' in consequence_terms)


class Filter(object):
    def __repr__(self):
        return "Filter({})".format(self._filters)

    def __init__(self, filters=()):
        self._filters = filters
        self.filters = [(compile(expression, "<string>", "eval"), expression) for expression in filters]

    def annotate(self, annotations):
        num = 0
        for annotation in annotations:
            passedfilters = True
            for compiled, expression in self.filters:
                try:
                    if not eval(compiled, globals(), annotation):
                        passedfilters = False
                        break
                except NameError:
                    pass
                except Exception as e:
                    raise CoverMiException("ERROR {} while filtering with {}".format(e, expression))
            annotation["passedfilters"] = passedfilters
            if passedfilters:
                num += 1
        return num

    def filter(self, annotations):
        return _Filter(annotations, self.filters)


class _Filter(object):
    def __init__(self, annotations, filters):
        self.annotations = annotations
        self.filters = filters

    def __iter__(self):
        for annotation in self.annotations:
            passedfilters = True
            for compiled, expression in self.filters:
                try:
                    if not eval(compiled, globals(), annotation):
                        passedfilters = False
                        break
                except NameError:
                    pass
                except Exception as e:
                    raise CoverMiException("ERROR {} while filtering with {}".format(e, expression))
            if passedfilters:
                yield annotation



#class Filter(object):
#    def __init__(self, rawfilters):
#        if isinstance(rawfilters, basestring):
#            rawfilters = (rawfilters,)
#        tryblocks = []
#        for rawfilter in rawfilters:
#            newtokens = []
#            for num, val, _, _, _ in tokenize.generate_tokens(StringIO(rawfilter).readline):
#                if num == tokenize.NAME and val not in ("or", "and", "not", "in", "set"):
#                    newtokens += [(tokenize.NAME, "variant"), (tokenize.OP, ".")]
#                newtokens += [(num,  val)]
#            tryblocks += ["  try:\n    if not("+tokenize.untokenize(newtokens)+"):\n      return False\n  except AttributeError:\n    pass\n"]
#        exec("def passfilters(variant):\n"+"".join(tryblocks)+"  return True\n")
#        self.passfilters = passfilters

#    def __call__(self, variants):
#        for variant in variants:
#            if self.passfilters(variant):
#                yield variant


#@keyvalcache
#def vepweb(variants, assembly="GRCh37", species="homo_sapiens", transcripts="refseq"):
#    try:
#        url = {"GRCh37": "http://grch37.rest.ensembl.org/vep/{}/region", "GRCh38": "http://rest.ensembl.org/vep/{}/region"}[assembly]
#    except KeyError:
#        raise CoverMiException("ERROR Unknown assembly {}".format(assembly))
#    url = url.format(species)

#    data = {"canonical": True, "hgvs": True}
#    if transcripts == "refseq":
#        data["refseq"] = True

#    variants = variants.__iter__()
#    while True:
#        batch = list(zip(range(0, 300), variants))
#        if len(batch) == 0:
#            break
#        data["variants"] = ["{} {} {} {}/{} + {}\n".format(str(v.chrom)[3:], v.start, v.stop, v.ref, v.alt, i+1) for i, v in batch]
#        r = requests.post(url, headers={"Content-Type": "application/json", "Accept": "application/json"}, data=json.dumps(data))
#        if not r.ok:
#            r.raise_for_status()
#        for vep_output in r.json():
#            print(vep_output)
#            yield (batch[int(vep_output["id"])-1][1], vep_output)


def vepscript(variants, assembly=DEFAULT_ASSEMBLY, species=DEFAULT_SPECIES, transcript_source=DEFAULT_TRANSCRIPTS):
    if not hasattr(variants, "__getitem__"):
        raise TypeError("'variants' argument to vepscript must implement __getitem__")

    extra = []
    operating_system = os.name
    if operating_system == "posix":
        basedir = os.path.expanduser("~")
    elif operating_system == "nt":
        basedir = "C:\\"
        extra = ["--dir", "C:\\cache"]
    else:
        raise RuntimeError("Unable to locate vep path on {}".format(operating_system))

    with tempfile.NamedTemporaryFile(delete=False) as temp:
        for i, v in enumerate(variants):
            row = "{} {} {} {}/{} + {}\n".format(str(v.chrom)[3:], v.start, v.stop, v.ref, v.alt, i+1)
            temp.write(row)
    
    scriptpath = os.path.join(basedir, "ensembl-vep-90", "vep")
    command = ["perl", scriptpath, "--json", "--offline", "--everything", "--no_stats", "--use_given_ref", "-o", "STDOUT", "-i", temp.name, \
                                                   "--assembly", assembly, "--species", species] + extra
    if transcript_source == "refseq":
        command += ["--refseq"]
    elif transcript_source != "ensembl":
        raise CoverMiError("ERROR Unknown transcript type {}".format(transcript_source))
    vep_script = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=False)
    for vep_output in vep_script.stdout:
        vep_output = json.loads(vep_output)
        yield (variants[int(vep_output["id"])-1], vep_output)
    vep_script.wait()
    os.unlink(temp.name)


#class Reference(object):
#    #__slots__ = ("id", "db")
# 
#    def __repr__(self):
#        return "{}({}, {})".format(type(self), repr(self.db), repr(self.id))

#    def __init__(self, db, refid):
#        self.id = refid
#        self.db = db

#    @property
#    def url(self):
#        if self.db == "dbSNP":
#            return "https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs={}".format(self.id)
#        elif self.db == "PUBMED":
#            return "https://www.ncbi.nlm.nih.gov/pubmed/{}".format(self.id)
#        elif self.db == "COSMIC":
#            return "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id={}".format(self.id[4:] if self.id.startswith("COSM") else self.id)
#        return None




def score_biotype(biotype):
    if biotype == "protein_coding":
        return 2
    elif biotype == "pseudogene":
        return 0
    else:
        return 1


class Annotation(dict):
    def __init__(self, variant=None):
        self.variant = variant

    def __missing__(self, key): 
        try:
            return getattr(self.variant, key)
        except AttributeError:
            raise KeyError

    def __eq__(self, other):
        try:
            return self.variant.__eq__(other.variant) 
        except AttributeError:
            return self.variant.__eq__(variant) 

    def __hash__(self):
        return self.variant.__hash__() 


# Yield "one annotation per variant" or
#       "one annotation per gene" or
#       "one annotation per specified transcript and one annotation per variant for all others" or
#       "everything"
def vep(variants, what="routine", panel=None, targets="", principal="", assembly=DEFAULT_ASSEMBLY, species=DEFAULT_SPECIES, transcript_source=DEFAULT_TRANSCRIPTS):
    ROUTINE = 0
    ONE_PER_GENE = 1
    ONE_PER_VARIANT = 2
    EVERYTHING = 3
    what = {"routine": ROUTINE, "one_per_gene": ONE_PER_GENE, "one_per_variant": ONE_PER_VARIANT, "everything": EVERYTHING}[what]

    if panel:
        if not targets and "targets" in panel:
            targets = panel.targets
        if not principal and "principal" in panel:
            principal = panel.principal
        if not transcript_source and "transcript_source" in panel.properties:
            transcript_source = panel.properties["transcript_source"] 
        if not assembly and "assembly" in panel.properties:
            assembly = panel.properties["assembly"]
        if not species and "species" in panel.properties:
            species = panel.properties["species"]

    transcript_ids, gene_symbols, _ = load_targets(targets)
    principal = load_principal(principal)

    kwargs = {}
    if transcript_source:
        kwargs["transcript_source"] = transcript_source
    if assembly:
        kwargs["assembly"] = assembly
    if species:
        kwargs["species"] = species

    if transcript_source == "ensembl":
        transcriptsort = lambda x: (principal[x["transcript_id"]], "canonical" in x, -int(x["transcript_id"][4:]))
    elif transcript_source == "refseq":
        transcriptsort = lambda x: (x["transcript_id"].startswith("N"), principal[x["transcript_id"]], "canonical" in x, -int(x["transcript_id"][3:]))
    else:
        raise CoverMiException("ERROR Unknown transcript source {}".format(transcript_source))
        
    annotations = []
    for variant, vep_output in vepscript(variants, **kwargs):

        demographics = {}
        for colocated in vep_output.get("colocated_variants", ()):
            allele_string = colocated["allele_string"]
            identifier = colocated["id"]
            if allele_string == "COSMIC_MUTATION":
                if "cosmic" not in demographics:
                    demographics["cosmic"] = []
                demographics["cosmic"] += [identifier]

            elif allele_string == "HGMD_MUTATION":
                if "hgmd" not in demographics:
                    demographics["hgmd"] = []
                demographics["hgmd"] += [identifier]

            else:
                if identifier.startswith("rs"):
                    if "dbsnp" not in demographics:
                        demographics["dbsnp"] = []
                    demographics["dbsnp"] += [identifier]
                
                if "pubmed" in colocated:
                    demographics["pubmed"] = str(colocated["pubmed"]).split(",")

                gnomad = []
                nhlbi = []
                thousand = []
                for series, maf, allele in ((gnomad, "gnomad_afr_maf", "gnomad_afr_allele"),
                                            (gnomad, "gnomad_amr_maf", "gnomad_amr_allele"),
                                            (gnomad, "gnomad_asj_maf", "gnomad_asj_allele"),
                                            (gnomad, "gnomad_eas_maf", "gnomad_eas_allele"),
                                            (gnomad, "gnomad_fin_maf", "gnomad_fin_allele"),
                                            (gnomad, "gnomad_nfe_maf", "gnomad_nfe_allele"),
                                            (gnomad, "gnomad_oth_maf", "gnomad_oth_allele"),
                                            (gnomad, "gnomad_sas_maf", "gnomad_sas_allele"),
                                            (nhlbi, "ea_maf", "ea_allele"),
                                            (nhlbi, "aa_maf", "aa_allele"),
                                            (thousand, "afr_maf", "afr_allele"),
                                            (thousand, "amr_maf", "amr_allele"),
                                            (thousand, "asn_maf", "asn_allele"),
                                            (thousand, "eur_maf", "eur_allele"),
                                            (thousand, "eas_maf", "eas_allele"),
                                            (thousand, "sas_maf", "sas_allele")):
                    try:
                        if colocated[allele] == variant.alt:
                            series += [colocated[maf]]
                    except KeyError:
                        pass
                if gnomad:
                    demographics["gnomad"] = max(gnomad)
                if nhlbi:
                    demographics["nhlbi"] = max(nhlbi)
                if thousand:
                    demographics["thousand"] = max(thousand)
                if any((gnomad, nhlbi, thousand)):
                    demographics["maf"] = max(gnomad + nhlbi + thousand)

        if "transcript_consequences" in vep_output:
            consequence_by_gene = defaultdict(list)
            for consequence in vep_output["transcript_consequences"]:
                if "gene_symbol" not in consequence:
                    consequence["gene_symbol"] = "LOC{}".format(consequence["gene_id"])
                consequence_by_gene[consequence["gene_symbol"]] += [consequence]

            affected_genes = set(consequence_by_gene)
            selected_genes = affected_genes & gene_symbols if (what in (ROUTINE, ONE_PER_VARIANT)) else affected_genes
            if selected_genes:
                gene_match = True
            else:
                selected_genes = affected_genes
                gene_match = False
            
            selected = []
            for gene in selected_genes:
                dedup = {}
                for consequence in consequence_by_gene[gene]:
                    transcript_id, consequence["transcript_version"] = consequence["transcript_id"].split(".")
                    if transcript_id not in dedup or int(consequence["transcript_version"]) > int(dedup[transcript_id]["transcript_version"]) :
                        consequence["transcript_id"] = transcript_id
                        dedup[transcript_id] = consequence

                potentials = [consequence for consequence in dedup.values() if consequence["transcript_id"] in transcript_ids] if transcript_ids else []
                transcript_match = bool(potentials)
                if not transcript_match:
                    potentials = dedup.values()

                if len(potentials) > 1 and (what in (ONE_PER_VARIANT, ONE_PER_GENE) or (what == ROUTINE and not transcript_match)):
                    potentials = sorted(potentials, key=transcriptsort, reverse=True)[:1]

                selected += potentials

            if len(selected) > 1 and (what == ONE_PER_VARIANT or (what == ROUTINE and not gene_match)):
                selected = sorted(selected, key=lambda x: (score_biotype(x["biotype"]), -int(x["gene_id"])), reverse=True)[:1]
                selected_genes = set([selected[0]["gene_symbol"]])
                
            other_genes = list(affected_genes - selected_genes)

            for consequence in selected:
                annotation = Annotation(variant)
                annotation.update(demographics)
                annotation["other_genes"] = other_genes
                annotation["most_severe_consequence"] = vep_output["most_severe_consequence"]
                annotation["gene_symbol"] = consequence["gene_symbol"]
                annotation["gene_id"] = consequence["gene_id"]
                annotation["transcript_id"] = consequence["transcript_id"]
                annotation["impact"] = consequence["impact"]
                annotation["consequence_terms"] = consequence["consequence_terms"]
                annotation["biotype"] = consequence["biotype"]
                for key in ("hgvsc", "hgvsp"):
                    try:
                        hgvs = consequence[key]
                    except KeyError:
                        continue
                    transcript, change = hgvs.split(":")
                    annotation[key] = "{}:{}".format(transcript.split(".")[0], change)
                for key in ("sift_score", "sift_prediction", "polyphen_score", "polyphen_prediction"):
                    try:
                        annotation[key] = consequence[key]
                    except KeyError:
                        pass
                annotations += [annotation]
    return annotations

        
















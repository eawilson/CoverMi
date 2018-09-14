import requests, json, pdb, csv, os, subprocess, tempfile, copy
from itertools import izip
import tokenize, StringIO
from cache import keyvalcache


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


class FileContext(object):
    def __init__(self, fn_fp_list, mode=None):
        self.fn_fp_list = fn_fp_list
        self.mode = mode

    def __enter__(self):
        if isinstance(self.fn_fp_list, basestring):
            self.fp = open(self.fn_fp_list, self.mode)
            return self.fp
        else:
            try:
                self.fn_fp_list.seek(0)
            except AttributeError:
                pass
            return self.fn_fp_list

    def __exit__(self, type_, value, traceback):
        if hasattr(self, "fp"):
            self.fp.close()

    
examplefilters = ["filters == 'PASS'",
                  "vaf > 0.2",
                  "maf < 0.05",
                  "impact in 'MODERATE/HIGH' or (gene_symbol == 'FTL' and '5_prime_UTR_variant' in consequence_terms)",
                 ]


class Filter(object):
    def __init__(self, rawfilters):
        if isinstance(rawfilters, basestring):
            rawfilters = (rawfilters,)
        tryblocks = []
        for rawfilter in rawfilters:
            newtokens = []
            for num, val, _, _, _ in tokenize.generate_tokens(StringIO.StringIO(rawfilter).readline):
                if num == tokenize.NAME and val not in ("or", "and", "not", "in"):
                    newtokens += [(tokenize.NAME, "variant"), (tokenize.OP, ".")]
                newtokens += [(num,  val)]
            tryblocks += ["  try:\n    if not("+tokenize.untokenize(newtokens)+"):\n      return False\n  except AttributeError:\n    pass\n"]
        exec("def passfilters(variant):\n"+"".join(tryblocks)+"  return True\n")
        self.passfilters = passfilters

    def __call__(self, variants):
        for variant in variants:
            if self.passfilters(variant):
                yield variant


@keyvalcache
def vepweb(variants):
    data = {"canonical": True, "refseq": True, "hgvs": True}
    variants = variants.__iter__()
    while True:
        batch = list(izip(range(0, 300), variants))
        if len(batch) == 0:
            break
        data["hgvs_notations"] = ["{} {}".format(v.vep_input, i+1) for i, v in batch]
        r = requests.post("http://grch37.rest.ensembl.org/vep/human/hgvs", headers={"Content-Type": "application/json", "Accept": "application/json"}, data=json.dumps(data))
        if not r.ok:
            r.raise_for_status()
        for vep_output in r.json():
            yield (batch[int(vep_output["id"])-1][1], vep_output)


@keyvalcache
def vepscript(variants):
    batch = list(variants)
    with FileContext(tempfile.TemporaryFile()) as temp:
        for i, v in enumerate(batch):
            temp.write("{} {}\n".format(v.vep_input, i+1))
        temp.seek(0)
        command = ["perl", "/home/ed/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl", "--quiet", "--json", 
                   "--no_progress", "--offline", "--refseq", "--fork", "4", "--everything", "--check_existing", "--no_stats", "--check_alleles", "-o", "STDOUT"]
        vep_script = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=temp, universal_newlines=True, shell=False)
        for vep_output in vep_script.stdout:
            vep_output = json.loads(vep_output)
            yield (batch[int(vep_output["id"])-1], vep_output)
        

class ErrorCallback():
    def __init__(self):
        self.errors = {}

    def __call__(self, errtype, variant, details):
        self.errors[variant] = (errtype, details)

    def show(self):
        for error in self.errors.values():
            print error


def dummycallback(*args, **kwargs):
    pass


def vep(variants, specified_transcripts=(), error_callback=dummycallback):
    for variant, vep_output in vepscript(variants):
        references = []
        colocateds = vep_output.get("colocated_variants")
        if colocateds is not None:
            for colocated in colocateds:
                allele_string = colocated["allele_string"]
                if allele_string.endswith("_MUTATION"): #COSMIC or HGMD
                    database = allele_string.split("_")[0]
                    if database in ("COSMIC", "HGMD"):
                        references.append([database, colocated["id"]])
                    else:
                        error_callback("Unknown id in colocated_variants", variant.identifier, str(colocated))
                else:
                    value = colocated.get("id")
                    if value is not None:
                        if value.startswith("rs"):
                            references.append(["dbSNP", value])
                        elif not value.startswith("TMP_ESP_"):
                            error_callback("Unknown id in colocated_variants", variant.identifier, str(colocated))

                    value = str(colocated.get("pubmed"))
                    if value is not None:
                        references.extend([["PUBMED", subval] for subval in value.split(",")])

                    value = colocated.get("phenotype_or_disease")
                    if value is not None:
                        variant.phenotype_or_disease = bool(value)
                    
                    number_of_alleles = len(allele_string.split("/"))
                    colocated["minor_maf"] = colocated.get("minor_allele_freq")
                    
                                      # exac
                    frequency_keys = ("exac_adj", "exac_afr", "exac_amr", "exac_eas", "exac_fin", "exac", "exac_nfe", "exac_oth", "exac_sas",
                                      # 1000 genome
                                      "afr", "amr", "eas", "eur", "minor", "sas",
                                      # NHLBI-ESP
                                      "aa", "ea")
                    for key in frequency_keys:
                        allele = colocated.get(key+"_allele")
                        freq = colocated.get(key+"_maf")
                        if freq is not None:
                            if isinstance(freq, basestring):
                                if freq.count(":") == 1 and allele is None:
                                    allele, freq = freq.split(":")
                                    if allele.strip("-AGCT") != "":
                                        error_callback("Invalid allele", variant.identifier, str(colocated))
                                        freq = None
                                try:
                                    freq = float(freq)
                                except ValueError:
                                    error_callback("Non numeric freq", variant.identifier, str(colocated))
                                    freq = None

                            if allele is not None:
                                if allele == variant.alt:
                                    setattr(variant, key, freq)
                                elif number_of_alleles == 2 and allele == variant.ref and freq is not None:
                                    setattr(variant, key, 1 - freq)
                            else:
                                error_callback("No allele specified", variant.identifier, str(colocated))
                                freq = None

                    minor = getattr(variant, "minor", None)
                    exac = getattr(variant, "exac", None)
                    if minor is not None or exac is not None:
                        variant.maf = minor if minor>exac else exac
                    else:
                        frequencies = [getattr(variant, key) for key in frequency_keys if hasattr(variant, key)]
                        if len(frequencies) > 0:
                            variant.maf = sum(frequencies)/len(frequencies)

        variant.references = references
        variant.input = vep_output.get("input")

        transcript_consequences = vep_output.get("transcript_consequences")
        most_severe_consequence = vep_output.get("most_severe_consequence")
        if transcript_consequences is not None:
            thesevariants = [variant]
            matches = [c for c in transcript_consequences if c["transcript_id"].split(".")[0] in specified_transcripts]
            if len(matches) == 0:
                #print "WARNING: No canonical transcript found for {}".format(variant.identifier)#################################################################3
                matches = transcript_consequences
            if len(matches) > 1:
                variant.alternate_transcripts = [variant]
                for _ in matches[1:]:
                    thisvariant = copy.copy(variant)
                    variant.alternate_transcripts += [thisvariant]
                thesevariants = variant.alternate_transcripts

            firstpass = True
            for thisvariant, consequence in izip(thesevariants, matches):                
                variant_allele = consequence.get("variant_allele")
                if variant_allele is None or variant_allele != variant.alt:
                    error_callback("Alt allele differs", variant.identifier, str(consequence))

                for key in ("sift_prediction", "sift_score", "polyphen_prediction", "polyphen_score", "gene_symbol", "impact", "biotype", "transcript_id", "hgvsc", "hgvsp"):
                    value = consequence.get(key) 
                    if value is not None:
                        setattr(thisvariant, key, value)
                consequence_terms = consequence.get("consequence_terms", ())
                if most_severe_consequence not in consequence_terms:
                    thisvariant.most_severe_consequence = most_severe_consequence
                if consequence_terms:
                    thisvariant.consequence_terms = " ".join(consequence_terms)
                yield thisvariant

        else:
            if most_severe_consequence is not None:
                variant.consequence_terms = most_severe_consequence
                variant.impact = IMPACT[most_severe_consequence]
            yield variant


def get_from_ncbi(db, items):
    server = "https://eutils.ncbi.nlm.nih.gov"
    ext = "/entrez/eutils/esummary.fcgi"
    params = {"db": db, "retmode": "json", "id": ",".join(items)}
    r = requests.post(server+ext, params=params)
    if not r.ok:
        r.raise_for_status()
    return r.json()

cached_ncbi = {}

def Ncbi(variants):
    refids = []
    variants = list(variants)
    for variant in variants:
        for db, refid in variant.references:
            if db == "dbSNP":
                if refid not in cached_ncbi:
                    refids += [refid.strip("rs")]
    if refids:
        response = get_from_ncbi("snp", refids)
        for key, snp in response["result"].items():
            if key != "uids" and "clinical_significance" in snp:
                cached_ncbi["rs"+key] = snp["clinical_significance"]
        
    for variant in variants:
        for db, refid in variant.references:
            if db == "dbSNP" and refid in cached_ncbi:#???????????????????????????????????????????? ?not returned ids?
                variant.clinical_significance = cached_ncbi[refid]
        yield variant



from reader import Hgmd

hgmd = Hgmd("/home/ed/Desktop/qway/panels/IPD v1/IPD Panel variants (08.11.16).txt")#"/home/ed/Desktop/icbd/ICBD v1/ICBD HGMD MUTATIONS.txt")#
variants = list(hgmd.read(Hgmd.MUTATION))
ec = ErrorCallback()
variants = set(vep(variants, error_callback=ec))
#print len(variants)

y = 0
n = 0
for variant in variants:
    transcripts = [v.hgvsc.split(":")[1].replace("dup", "ins") for v in getattr(variant, "alternate_transcripts", (variant,))  if hasattr(v, "hgvsc")]
    
    if variant.name.replace("dup", "ins") in transcripts:
        y += 1        
    else:
        print transcripts, variant.row
        n += 1
print y, n









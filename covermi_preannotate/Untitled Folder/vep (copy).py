import requests, json, pdb, csv, os, subprocess, tempfile, copy
from itertools import izip
from collections import namedtuple
import tokenize, StringIO
from covermi import reader

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

class Alternate(object):

    def __init__(self, variant):
        self.variant = variant

    def __getattr__(self, name):
        return getattr(self.variant, name)


class Variant(object):

    homozygous_vaf = (.9, 1)
    heterozygous_vaf = (.35, .65)

    def __repr__(self):
        return self.identifier

    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            if attr == "qual":
                try:
                    self.qual = float(val)
                except ValueError:
                    pass
            elif attr in ("pos", "depth", "alt_depth"):
                try:
                    setattr(self, attr, int(val))
                except ValueError:
                    pass
            else:
                setattr(self, attr, intern(val) if attr in ("chrom", "ref", "alt", "filters") else val)
        self.identifier = "{}:{} {}>{}".format(self.chrom, self.pos, self.ref, self.alt)     

    @property
    def vaf(self):
        return float(self.alt_depth)/self.depth

    @property
    def start(self):
        return self.pos

#    @property
#    def stop(self): #?????
#        return self.start+1 if len(self.ref)<=len(self.alt) else self.start+len(self.alt)-1

    @property
    def ensemblref(self):
        return self.ref if len(self.ref) == len(self.alt) else self.ref[1:].ljust(1, "-")

    @property
    def ensemblalt(self):
        return self.alt if len(self.ref) == len(self.alt) else self.alt[1:].ljust(1, "-")

    @property
    def type(self):
        if len(self.ref) == len(self.alt):
            return "snp"
        return "del" if len(self.ref) > len(self.alt) else "ins"

    @property
    def zygosity(self):
        vaf = self.vaf
        if type(self).homozygous_vaf[0] <= vaf <= type(self).homozygous_vaf[1]:
            return "hemizygous" if self.chrom in ("chrX", "chrY") else "homozygous"
        elif type(self).heterozygous_vaf[0] <= vaf <= type(self).heterozygous_vaf[1]:
            return "heterozygous"
        else:
            return "unknown"
        
    def __hash__(self):
        return hash(self.identifier)

    def __eq__(self, other):
        try:
            return self.identifier == other.identifier
        except AttributeError:
            return self.identifier == other


class BaseVariant(Entry):
    __slots__ = ("pos", "ref", "alt")

    def __repr__(self):
        return "{}:{} {}/{}".format(self.chrom, self.pos, self.ref, self.alt)

    def __init__(self, chrom, pos, ref, alt, name="."):
        lenref = len(ref) if self.ref!="-" else 0
        if lenref == (len(alt) if self.alt!="-" else 0):
            start = pos
            stop = pos + lenref - 1
        else:
            start = pos - 1
            stop = pos + lenref
        super(BaseVariant, self).__init__(chrom, start, stop, name)

        self.pos = pos
        self.ref = ref
        self.alt = alt

    @property
    def identifier(self):
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def vep_input(self):
        return "{} {} {} {}/{} +".format(self.chrom[3:], self.pos, self.pos + (len(self.ref) if self.ref!="-" else 0) - 1, self.ref, self.alt)

    @property
    def vartype(self):
        if self.ref == "-":
            return "ins"
        elif self.alt == "-":
            return "del"
        elif len(self.ref) == len(self.alt):
            return "snp"
        return "del" if len(self.ref) > len(self.alt) else "ins"


class SequencedVariant(BaseVariant):
    __slots__ = ("qual", "depth", "alt_depth", "filters")
    homozygous_vaf = (.9, 1)
    heterozygous_vaf = (.35, .65)

    @property
    def vaf(self):
        return float(self.alt_depth)/self.depth

    @property
    def zygosity(self):
        vaf = self.vaf
        if type(self).homozygous_vaf[0] <= vaf <= type(self).homozygous_vaf[1]:
            return "hemizygous" if self.chrom in ("chrX", "chrY") else "homozygous"
        elif type(self).heterozygous_vaf[0] <= vaf <= type(self).heterozygous_vaf[1]:
            return "heterozygous"
        else:
            return "unknown"


class Variant(SequencedVariant):
    pass










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


def vep(variantlist):
    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/human/hgvs"
    headers={ "Content-Type": "application/json", "Accept": "application/json"}
    data = {"hgvs_notations": variantlist, "canonical": True, "refseq": True, "hgvs": True}
    r = requests.post(server+ext, headers=headers, data=json.dumps(data))
    if not r.ok:
        r.raise_for_status()
    return r.json()


def Vcf(path):
    ILLUMINA = 0
    IONTORRENT = 1
    vcfstyle = None
    with FileContext(path, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row[0].startswith("#"):
                chrom, pos, ident, ref, alts, qual, filters, info, format, sample = row[:10]
                pos = int(pos)

                if vcfstyle is None:
                    if any(key.startswith("AO=") for key in info.split(";")):
                        vcfstyle = IONTORRENT
                    elif "AD" in format.split(":"):
                        vcfstyle = ILLUMINA
                
                filters = ";".join(code for code in filters.split(";") if code != "LowVariantFreq")
                if filters == "":
                    filters = "PASS"
                if vcfstyle == IONTORRENT:
                    info = info.split(";")
                    keys = [val.split("=")[0] for val in info]
                    depth = int(info[keys.index("DP")].split("=")[1])
                    alt_depths = map(int, info[keys.index("AO")].split("=")[1].split(",")) #???????????????????????????????????????
                elif vcfstyle == ILLUMINA:
                    ad = map(int, sample.split(":")[format.split(":").index("AD")].split(","))
                    depth = max(sum(ad), 1)
                    alt_depths = ad[1:]

                if "," in alts:
                    print row

                for alt, alt_depth in zip(alts.split(","), alt_depths):
                    if alt_depth > 0 and ref != alt and alt != ".":
                        thisref = ref
                        thispos = pos
                        while thisref[:2] == alt[:2]:
                            thisref = thisref[1:]
                            alt = alt[1:]
                            thispos += 1
                        while thisref[-1] == alt[-1] and len(thisref) > 1 and len(alt) > 1:
                            thisref = thisref[:-1]
                            alt = alt[:-1]
                        yield Variant(chrom=chrom, pos=thispos, ref=thisref, alt=alt, alt_depth=alt_depth, depth=depth, filters=filters, qual=qual)


def Vcf(path, slots=False):

    def fao_fdp_depths():
        data = dict(keyval.split("=") for keyval in info.split(";"))
        yield int(data["FDP"])
        for d in map(int, data["FAO"].split(",")):
            yield d

    def ao_dp_depths():
        data = dict(keyval.split("=") for keyval in info.split(";"))
        yield int(data["DP"])
        for d in map(int, data["AO"].split(",")):
            yield d

    def ad_depths():
        ad = map(int, sample.split(":")[format.split(":").index("AD")].split(","))
        yield sum(ad)
        for d in ad[1:]:
            yield d

    with FileContext(path, "rb") as f:
        firstpass = True
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row[0].startswith("#"):

                if firstpass:
                    infokeys = [keyval.split("=")[0] for keyval in row[7].split(";")] if len(row)>=8 else ()
                    formatkeys = row[8].split(":") if len()>=10 else ()
                    if "FAO" in infokeys and "FDP" in infokeys:
                        depthsfunc = fao_fdp_depths
                    elif "AO" in infokeys and "DP" in infokeys:
                        depthsfunc = ao_dp_depths
                    elif "AD" in formatkeys:
                        depthsfunc = ad_depths
                    else:
                        depthsfunc = None

                    hasdepth = bool(depthsfunc)
                    hasfilters = bool(row[6] != ".")
                    hasqual = bool(row[6] != ".")
                    variantclass = (SequencedVariant if (hasdepth or hasfilters or hasqual) else BaseVariant) if slots else Variant
                    firstpass = False

                chrom = row[0]
                if hasfilters:
                    filters = ";".join(code for code in row[6].split(";") if code != "LowVariantFreq")
                    if filters == "":
                        filters = "PASS"
                if hasqual:
                    qual = float(row[5])
                if hasdepth:
                    depth = depthsfunc.next()
                basepos = int(row[1])
                baseref = row[3]
                alts = row[4].split(",")

                for alt in alts:
                    if alt not in (ref, "."):
                        if depthsfunc:
                            altdepth = depthsfunc.next()
                            if altdepth == 0:
                                continue

                        ref = baseref
                        pos = basepos
                        while ref[:1] == alt[:1]:
                            ref = ref[1:]
                            alt = alt[1:]
                            pos += 1
                        while ref[-1:] == alt[-1:]:
                            ref = ref[:-1]
                            alt = alt[:-1]

                        variant = variantclass(chrom, pos, ref or "-", alt or "-")#?name
                        if hasqual:
                            variant.qual = qual
                        if hasfilters:
                            variant.filters = filters
                        if hasdepth:
                            variant.depth = depth
                            variant.alt_depth = alt_depth
                        yield variant



def VepWeb(variants, **kwargs):
    variants = variants.__iter__()
    while True:
        batch = list(izip(range(0, 300), variants))
        if len(batch) == 0:
            break
        for vep_output in vep(["{} {} {} {} {} . PASS".format(v.chrom, v.pos, i+1, v.ref, v.alt) for i, v in batch]):
            for variant in annotate(batch[int(vep_output["id"])-1][1], vep_output, **kwargs):
                yield variant



def VepScript(variants, **kwargs):
    with FileContext(tempfile.TemporaryFile()) as temp:
        batch = list(variants)
        for i, v in enumerate(batch):
            temp.write("{} {} {} {} {} . PASS\n".format(v.chrom, v.pos, i+1, v.ref, v.alt))
        temp.seek(0)
        command = ["perl", "/home/ed/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl", "--quiet", "--json", 
                   "--no_progress", "--offline", "--refseq", "--fork", "4", "--everything", "--check_existing", "--no_stats", "--check_alleles", "-o", "STDOUT"]
        vep_script = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=temp, universal_newlines=True, shell=False)
        for vep_output in vep_script.stdout:
            vep_output = json.loads(vep_output)
            for variant in annotate(batch[int(vep_output["id"])-1], vep_output, **kwargs):
                yield variant

Vep = VepScript
        


    
examplefilters = ["filters == 'PASS'",
                  "vaf > 0.2",
                  "maf < 0.05",
                  "impact in 'MODERATE/HIGH' or (gene_symbol == 'FTL' and '5_prime_UTR_variant' in consequence_terms)",
                 ]
printed = set()
def dummy(*args, **kwargs):
    return

def annotate(variant, vep_output, specified_transcripts=(), error_callback=dummy):
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

        matches = [c for c in transcript_consequences if c["transcript_id"].split(".")[0] in specified_transcripts]
        if len(matches) == 0:
            print "WARNING: No canonical transcript found for {}".format(variant.identifier)
            matches = transcript_consequences
        if len(matches) > 1:
            print "WARNING: Multiple transcripts found for {}".format(variant.identifier)

        for consequence in matches:
            thisvariant = Alternate(variant)
            
#                variant_allele = consequence.get("variant_allele")
#                if variant_allele is None or variant_allele != variant["ensemblalt"]:
#                    error_callback("Alt allele differs", variant["identifier"], str(consequence))

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
#    ????????????????if not hasattr(variant, "impact"):
#    ????????????????    variant.impact = "MODIFIER"
    



#        consequences = annotation.get("transcript_consequences")
#        if consequences:
#            nmnr_consequences = [consequence for consequence in consequences if consequence["transcript_id"][:3] in ("NM_", "NR_")]
#            if not nmnr_consequences:
#                error_callback(["No NM or NR transcripts", variant["identifier"], str(consequences)])
#                continue # Need to be careful not to skip important information
#            nmnr_consequences = sorted(nmnr_consequences, key=lambda cons: {"HIGH": 3, "MODERATE": 2 , "LOW": 1, "MODIFIER": 0}[cons["impact"]])

#            specified_index = 0
#            canonical_index = 0
#            for index, consequence in enumerate(nmnr_consequences, start=1):
#                transcript_id = consequence.get("transcript_id")
#                if transcript_id is not None and transcript_id.split(".")[0] in specified_transcripts:
#                    specified_index = index
#                if consequence.get("canonical"):
#                    canonical_index = index

#            index = specified_index or canonical_index or index
#            consequence = consequences[index-1]
#            if not(specified_index or canonical_index) and len(nmnr_consequences) > 1:
#                error_callback(["Multiple transcripts", variant["identifier"], str(consequences)])
#                variant["multiple_transcripts"] = True

#            variant_allele = consequence.get("variant_allele")
#            if variant_allele is None or variant_allele != variant["ensemblalt"]:
#                error_callback(["Alt allele differs", variant["identifier"], str(consequence)])

#            for key in ("sift_prediction", "sift_score", "polyphen_prediction", "polyphen_score", "gene_symbol", "impact", "biotype", "transcript_id", "hgvsc", "hgvsp"):
#                value = consequence.get(key) 
#                if value is not None:
#                    variant[key] = value
#            consequence_terms = consequence.get("consequence_terms")
#            if consequence_terms is not None:
#                variant["consequence_terms"] = " ".join([term.replace(" variant", "").replace(" ", "_") for term in consequence_terms])

#            variant["references"] = references
#            if pass_filters(variant, userfilters):
#                yield (variant, unexpected)
#            elif unexpected:
#                yield (None, unexpected)



def VepTest(variants):
    with FileContext(tempfile.TemporaryFile()) as temp:
        batch = list(variants)
        for i, v in enumerate(batch):
            temp.write("{} {} {} {} {} . PASS\n".format(v.chrom, v.pos, i+1, v.ref, v.alt))
        temp.seek(0)
        command = ["perl", "/home/ed/ensembl-tools-release-87/scripts/variant_effect_predictor/variant_effect_predictor.pl", "--quiet", "--json", 
                   "--no_progress", "--offline", "--refseq", "--fork", "4", "--everything", "--check_existing", "--no_stats", "--check_alleles", "-o", "STDOUT"]
        vep_script = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=temp, universal_newlines=True, shell=False)
        for vep_output in vep_script.stdout:
            vep_output = json.loads(vep_output)
            yield (batch[int(vep_output["id"])-1], vep_output)


#from covermi import Panel
#if __name__ == "__main__":
#    panel = Panel("/home/ed/Desktop/icbd/ICBD v1")
#    specified_transcripts = set(transcript.transcript.split()[1] for transcript in panel.transcripts)
#    for vep in Vep(Vcf("/home/ed/Desktop/icbd/test.vcf"), specified_transcripts=specified_transcripts):
#        print vep


   #annotate("/home/ed/Desktop/vcf/GOSH14DW1015PTBMMARCH.vcf")
    
#    with open("/media/ed/Transcend/Ion Torrent/G184155V_VCF_v1/TSVC_variants.vcf", "wb") as f:
        
#        for variant in VepScript("/home/ed/Desktop/qway/miseq/TSVC_variants.vcf"):
#            print variant
#            n += 1
#        print n
#            json.dump(variant, f)
#            f.write("\n\n")


#    scriptf = open("scriptvep.txt", "wt")
#    webf = open("webvep.txt", "wt")
#    w = {}
#    s = {}
#    for web, script in izip(VepWeb("/home/ed/Desktop/qway/miseq/G167124B.vcf"), VepScript("/home/ed/Desktop/qway/miseq/G167124B.vcf")):
#        w[web[0]["identifier"]] = web
#        s[script[0]["identifier"]] = script
#        scriptf.write(str(web)+"\n")
#        webf.write(str(script)+"\n")

#    scriptf.close()
#    webf.close()

#    for key in s.keys():
#        if  len(set(w[key][1].keys()) ^ set(s[key][1].keys())) > 1:
#            print w[key]
#            print s[key]











#    rundir = "/home/ed/Desktop/qway/cache/16MISEQ025"
#    for fn in os.listdir(rundir):
#        if fn.endswith(".vcf"):
#            for variant, unexpected in annotate(os.path.join(rundir, fn), userfilters=[], specified_transcripts=[]):#"NM_001203248"])
#                print fn        


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



import sys, os, tkFileDialog, Tkinter, csv, pdb

rootwindow = Tkinter.Tk()
rootwindow.withdraw()

print("Please select a folder of vcfs")
path = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a folder of vcfs')
if not bool(path):
    sys.exit()
path = os.path.abspath(path)
print("{0} folder selected".format(os.path.basename(path)))

with open(os.path.join(path, "canonical.txt")) as f:
    canonical = set(row[1] for row in csv.reader(f, delimiter="\t"))      

for vcffile in os.listdir(path):
    if vcffile.endswith(".vcf"):
        print vcffile
        vcfpath = os.path.join(path, vcffile)
        variants = list( VepScript(Filter(["qual > 180", "filters == 'PASS'", "vaf >= 0.04"])(Vcf(vcfpath)), specified_transcripts=canonical) )
        
        with open(os.path.splitext(vcfpath)[0]+"_annotated.tsv", "wb") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Variant", "Depth", "Alt Depth", "VAF", "Qual", "Gene", "HGVSc", "HGVSp", "Impact", "Consequences", "Most Severe Consequence", "Input"])
            for v in variants:
                writer.writerow([getattr(v, name, "") for name in ("identifier", "depth", "alt_depth", "vaf", "qual", "gene_symbol", "hgvsc", "hgvsp", "impact", "consequence_terms", "most_severe_consequence", "input")])








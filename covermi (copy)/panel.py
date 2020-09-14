# panel is a dict containing all of the data structures that define a panel
# 	Amplicons:		genomic range
# 	Exons:			genomic range
# 	Transcripts:		genomic range
# 	Depth:
#	Variants_Disease:	genomic range
#	Variants_Gene:		genomic range
#	Variants_Mutation:	genomic range
#	Filenames:		dict of all files in the panel directory with filetype as the key
#       Options:		dict of all options, including depth

# the following are only included for a design panel
#	 AllTranscripts: 	genomic range
#	 AllExons:		genomic range
#	 Excluded:		list of excluded amplicons
from __future__ import print_function, absolute_import, division

import os, re, pdb
from functools import wraps

from .gr import Gr, bed, vcf, illuminamanifest, variants, reference, load_targets, load_principal, HETERO_LL, HETERO_UL, HOMO_LL
from .include import *

                                                   #Eleven columns - refflat
REGEXPS = (("reference",     "refseq",  re.compile(GENE_SYMBOL+"\t"+REFSEQ_TRANSCRIPT+"\tchr.+?\t[+-]\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9,]+\t[0-9,]+$")), 
           ("reference",     "ensembl", re.compile("#!genome-build")), #Ensembl gtf  
           ("targets",       "",        re.compile(GENE_SYMBOL+"$")), #Single column - targets
           ("targets",       "refseq",  re.compile(GENE_SYMBOL+" +"+REFSEQ_TRANSCRIPT+"$")), #Single column - targets
           ("targets",       "ensembl", re.compile(GENE_SYMBOL+" +"+ENSEMBL_TRANSCRIPT+"$")), #Single column - targets
#           ("canonical",     "refseq",  re.compile(GENE_SYMBOL+"\t"+REFSEQ_TRANSCRIPT+"$")), #Two column - canonical
#           ("canonical",     "ensembl", re.compile(GENE_SYMBOL+"\t"+ENSEMBL_TRANSCRIPT+"$")), #Two column - canonical
           ("manifest",      "",        re.compile("\\[Header\\]$")), #Manifest
           ("variants",      "",        re.compile("HGMD ID.+Disease")), #HGMD
           ("variants",      "",        re.compile("Gene name.+Primary site.+Primary histology")), #Cosmic
           ("designstudio",  "",        re.compile("chr[0-9XYM]+\t[0-9]+\t[0-9]+")), #Six+ column Bedfile - design
           ("diseases",      "",        re.compile("#diseases")), #Disease_Names
           ("vcf",           "",        re.compile("##fileformat=VCFv[0-9\.]+$")), #VCF
           #("properties",    "",        re.compile("[^#].*=.+")), # CoverMi panel options
           ("properties",    "",        re.compile("#covermi")), # CoverMi panel options
           ("principal",     "refseq",  re.compile(GENE_SYMBOL+"\t[0-9]+\t"+REFSEQ_TRANSCRIPT+".+\t(PRINCIPAL|ALTERNATIVE):")),
           ("principal",     "ensembl", re.compile(GENE_SYMBOL+"\tENSG[0-9]+\t"+ENSEMBL_TRANSCRIPT+".+\t(PRINCIPAL|ALTERNATIVE):")),
          )

ASSEMBLIES = (("hg19",   "GRCh37"),
              ("grch37", "GRCh37"),
              ("hg38",   "GRCh38"),
              ("grch38", "GRCh38"),
             )


class StrObj(str):
    pass


def Panels(path, **kwargs):
    panels = {}
    for name in os.listdir(path):
        panel_path = os.path.join(path, name)
        if os.path.isdir(panel_path):
            try:
                panels[name] = Panel(panel_path, **kwargs)
            except CoverMiException:
                pass
    return panels


def cached(func):
    @wraps(func)
    def wrapped(self, *args, **kwargs):
        cachedname = "_"+func.__name__
        if not hasattr(self, cachedname):
            setattr(self, cachedname, func(self, *args, **kwargs))
        return getattr(self, cachedname)
    return wrapped


class Panel(object):

    def __init__(self, path, required=(), verbose=True, splice_site_buffer=0):
        self.path = os.path.abspath(path)
        self._eprint = eprint if verbose else lambda *args: None
        self.splice_site_buffer = splice_site_buffer
        self.files = {}
        self.name = os.path.basename(path)
        for fn in os.listdir(path):
            full_path = os.path.join(path, fn)
            if os.path.isfile(full_path):
                filetype = None
                loops = 0
                with open(full_path, "rt") as f:
                    for testrow in f.read(1000).split("\n"): # Don't get screwed by really big binary files
                        testrow = testrow.strip()
                        loops += 1
                        for thisfiletype, thissource, regexp in REGEXPS:
                            match = regexp.match(testrow)
                            if match:
                                if match and filetype not in (thisfiletype, None):
                                    raise CoverMiException("ERROR. File {} matches both {} and {} format".format(os.path.basename(full_path), thisfiletype, filetype))
                                filetype = StrObj(thisfiletype)
                                filetype.transcript_source = thissource
                        if loops == 2 or match:
                            break
                if filetype is not None:
                    if filetype in self.files:
                        raise CoverMiException("ERROR. in {} panel - {} and {} are both of {} type".format(self.name, os.path.basename(self.files[filetype]), fn, filetype))
                    self.files[filetype] = os.path.abspath(full_path)

        if not self.files:
            raise CoverMiException("ERROR. {} panel is empty".format(self.name))
        for filetype in required:
            if filetype not in self.files:
                raise CoverMiException("ERROR. {} panel does not contain a {} file".format(self.name, filetype))


    @property
    @cached
    def targets(self):
        if "targets" in self.files:
            self._eprint("Loading targets file: {}".format(os.path.basename(self.files["targets"])))
            return load_targets(self.files["targets"])
        else:
            raise CoverMiException("ERROR. {} panel has no targets file".format(self.name))


    @property
    @cached
    def principal(self):
        if "principal" in self.files:
            self._eprint("Loading principal transcripts list: {}".format(os.path.basename(self.files["principal"])))
            return load_principal(self.files["principal"])
        else:
            raise CoverMiException("ERROR. {} panel has no principal transcripts file".format(self.name))


    @property
    @cached
    def amplicons(self):
        if "manifest" in self.files and "designstudio" in self.files:
            raise CoverMiException("ERROR. {} panel contains both manifest and design studio bedfile".format(self.name))
        elif "manifest" in self.files:
            self._eprint("Loading manifest: {}".format(os.path.basename(self.files["manifest"])))
            return Gr(illuminamanifest(self.files["manifest"]))
        elif "designstudio" in self.files:
            self._eprint("Loading design studio bedfile: {}".format(os.path.basename(self.files["designstudio"])))
            return Gr(bed(self.files["designstudio"]))
        else:
            raise CoverMiException("ERROR. {} panel does not contain a manifest or design studio bedfile".format(self.name))


    @property
    @cached
    def offtarget(self):
        if "manifest" in self.files:
            self._eprint("Loading manifest: {}".format(os.path.basename(self.files["manifest"])))
            return Gr(illuminamanifest(self.files["manifest"], ontarget=False, offtarget=True))
        else:
            raise CoverMiException("ERROR. {} panel does not contain a manifest".format(self.name))
        

    @property
    @cached
    def expected(self):
        return self.amplicons.combined_with(self.offtarget)


    @property
    @cached
    def unexpected(self):
        return self.expected.inverted


    def _reference(self, what, everything=False):
        if "reference" in self.files:
            self._eprint("Loading reference file: {}".format(os.path.basename(self.files["reference"])))
            kwargs = {"splice_site_buffer": self.splice_site_buffer}
            if "principal" in self:
                kwargs["principal"] = self.principal

            if everything:
                gr = Gr(reference(self.files["reference"], what, **kwargs))
            elif "targets" in self.files:
                gr = Gr(reference(self.files["reference"], what, targets=self.files["targets"], **kwargs))
            elif "amplicons" in self:
                self._eprint("WARNING. {} panel contains no targets file. Loading all genes touching an amplicon".format(self.name))
                if what == "transcripts":
                    gr = self.alltranscripts.touched_by(self.amplicons)
                elif what == "codingregions":
                    gr = self.allcodingregions.touched_by(self.amplicons)
                elif what == "exons":
                    gr = self.allexons.touched_by(self.alltranscripts.touched_by(self.amplicons))
                elif what == "codingexons":
                    gr = self.allcodingexons.touched_by(self.allcodingregions.touched_by(self.amplicons))
                else:
                    raise ValueError("Bad argumnent {} to _reference in Panel".format(waht))
            else:
                self._eprint("WARNING. {} panel contains no targets or amplicons files. Loading all genes".format(self.name))
                gr = getattr(self, "all{}".format(what))
            return gr
        else:
            raise CoverMiException("ERROR. {} panel does not contain a reference file".format(self.name))


    @property
    @cached
    def transcripts(self):
        return self._reference("transcripts")


    @property
    @cached
    def codingregions(self):
        return self._reference("codingregions")


    @property
    @cached
    def exons(self):
        return self._reference("exons") 


    @property
    @cached
    def codingexons(self):
        return self._reference("codingexons") 


    @property
    @cached
    def targeted_transcripts(self):
        try:
            return self.transcripts.overlapped_by(self.amplicons)
        except CoverMiException:
            return self.transcripts


    @property
    @cached
    def targeted_codingregions(self):
        try:
            return self.codingregions.overlapped_by(self.amplicons)
        except CoverMiException:
            return self.codingregions


    @property
    @cached
    def targeted_exons(self):
        try:
            return self.exons.overlapped_by(self.amplicons)
        except CoverMiException:
            return self.exons


    @property
    @cached
    def targeted_codingexons(self):
        try:
            return self.codingexons.overlapped_by(self.amplicons)
        except CoverMiException:
            return self.codingexons


    @property
    @cached
    def alltranscripts(self):
        return self._reference("transcripts", everything=True) 
 

    @property
    @cached
    def allcodingregions(self):
        return self._reference("codingregions", everything=True)


    @property
    @cached
    def allexons(self):
        return self._reference("exons", everything=True) 


    @property
    @cached
    def allcodingexons(self):
        return self._reference("codingexons", everything=True) 


    def _variants(self, what):
        if "variants" in self.files:
            self._eprint("Loading variants file: {}".format(os.path.basename(self.files["variants"])))
            kwargs = {}
            if "diseases" in self.files:
                self._eprint("Loading diseases file: {}".format(os.path.basename(self.files["diseases"])))
                kwargs["diseases"] = self.files["diseases"]
            if what != "disease":
                kwargs["genes"] = [entry.gene for entry in self.transcripts]

            gr = Gr(variants(self.files["variants"], what, **kwargs))

            if what != "disease" and gr and "transcripts" in self:
                variants_in_correct_location = gr.touched_by(self.transcripts).components * 100 // gr.components
                if variants_in_correct_location < 98:
                    self._eprint("WARNING. Only {}% of variants are within targeted genes. ?Correct reference genome".format(variants_in_correct_location))
            return gr
        else:
            raise CoverMiException("ERROR. {} panel does not contain a variants file".format(self.name))

    @property
    @cached
    def variants_disease(self):
        return self._variants("disease")


    @property
    @cached
    def variants_gene(self):
        return self._variants("gene")


    @property
    @cached
    def variants_mutation(self):
        return self._variants("mutation")


    @property
    @cached
    def targeted_variants_disease(self):
        try:
            return self.variants_disease.subranges_covered_by(self.amplicons)
        except CoverMiException:
            return self.variants_disease


    @property
    @cached
    def targeted_variants_gene(self):
        try:
            return self.variants_gene.subranges_covered_by(self.amplicons)
        except CoverMiException:
            return self.variants_gene


    @property
    @cached
    def targeted_variants_mutation(self):
        try:
            return self.variants_mutation.subranges_covered_by(self.amplicons)
        except CoverMiException:
            return self.variants_mutation


    @property
    def depth(self):
        try:
            return int(self.properties["depth"])
        except KeyError:
            raise CoverMiException("ERROR {} panel has no depth".format(self.name))
        except ValueError:
            raise CoverMiException("ERROR {} panel has a non-numeric depth".format(self.name))


    @property
    @cached
    def properties(self):
        properties = {"name": self.name}

        def update(key, val): # raise error if different values for same property
            if val:
                if key in properties and properties[key] != val:
                    raise CoverMiException("ERROR. {} panel has values of {} and {} for {}".format(properties.name, properties[key], val, key))
                properties[key] = val

        for filetype, filepath in self.files.items():
            filename = os.path.basename(filepath).lower()
            for text, assembly in ASSEMBLIES:
                if text in filename:
                    update("assembly", assembly)
            update("transcript_source", filetype.transcript_source)

        if "properties" in self.files:
            self._eprint("Loading properties file: {}".format(os.path.basename(self.files["properties"])))
            with open(self.files["properties"], "rU") as f:
                key = None
                for row in f:
                    row = row.rstrip()
                    lsrow = row.lstrip()
                    if row and lsrow[0] != "#":
                        if row == lsrow:
                            pos = lsrow.find("=")
                            if pos < 1 or pos == len(lsrow) - 1:
                                self._eprint("Malformed properties file: {}".format(row))
                                continue
                            else:
                                key = lsrow[:pos].lower()
                                update(key, lsrow[pos+1:])
                        elif key:
                            if isinstance(properties[key], list):
                                properties[key] += [lsrow]
                            else:
                                properties[key] = [properties[key], lsrow]
                        else:
                            self._eprint("Malformed properties file: {}".format(row))
                            continue

        for key, val in (("hetero_ll", HETERO_LL), ("hetero_ul", HETERO_UL), ("homo_ll", HOMO_LL)):
            if key not in properties:
                update(key, val)

        return properties


    def __contains__(self, attr):
        try:
            test = getattr(self, attr)
        except CoverMiException:
            return False
        return True


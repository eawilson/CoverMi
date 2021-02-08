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
import os
import re
import pdb
import gzip
from itertools import chain

from .gr import Gr, bed, gff3



def gzopen(fn, *args, **kwargs):
    return (gzip.open if fn.endswith(".gz") else open)(fn, *args, **kwargs)



REFSEQ_TRANSCRIPT = r"[NX][MR]_[0-9]+(\.[0-9])?"
ENSEMBL_TRANSCRIPT = r"ENST[0-9]{11}"
GENE_SYMBOL = r"[A-Z][A-Zorf0-9-\.]+"

#Eleven columns - refflat
REGEXPS = (("reference", re.compile("##gff-version 3")),
           ("reference", re.compile("chr[0-9a-zA-Z]+\t.+\t.+\t[0-9]+\t[0-9]+\t.+\t[+-\\.]+\t.+\t[^\t]+$")), # gff file
           ("targets",   re.compile("chr[0-9a-zA-Z]+\t[0-9]+\t[0-9]+")), # Bedfile
           ("principal", re.compile(f"{GENE_SYMBOL}\t[0-9]+\t{REFSEQ_TRANSCRIPT}.+\t(PRINCIPAL|ALTERNATIVE):")),
           ("principal", re.compile(f"{GENE_SYMBOL}\tENSG[0-9]+\t{ENSEMBL_TRANSCRIPT}.+\t(PRINCIPAL|ALTERNATIVE):")),
           ("names",     re.compile(f"{GENE_SYMBOL}$")), #Single column
           ("names",     re.compile(f"{GENE_SYMBOL}( {REFSEQ_TRANSCRIPT})+$")), #Single column
           ("names",     re.compile(f"{GENE_SYMBOL}( {ENSEMBL_TRANSCRIPT})+$")), #Single column
          )



class Panel(object):
    def __init__(self, panel_path):
        self._data = {}
        self.path = panel_path
        self.paths = {}
        for fn in os.listdir(panel_path):
            path = os.path.join(panel_path, fn)
            if os.path.isfile(path):
                matchedfiletype = ""
                with gzopen(path, "rt") as f:
                    try:
                        for testrow in f.read(1000).splitlines()[:2]: # Don't get screwed by really big binary files
                            testrow = testrow.strip()
                            for filetype, regexp in REGEXPS:
                                if regexp.match(testrow):
                                    if not matchedfiletype:
                                        try:
                                            self.paths[filetype].append(os.path.abspath(path))
                                        except KeyError:
                                            self.paths[filetype] = [os.path.abspath(path)]
                                    elif filetype != matchedfiletype:
                                        raise RuntimeError(f"file {fn} matches both {filetype} and {matchedfiletype} formats")
                                    matchedfiletype = filetype                                    
                    except UnicodeDecodeError:
                        pass
                    
        if not self.paths:
            raise RuntimeError("panel {} is empty".format(os.path.basename(panel_path)))
    
    
    def __repr__(self):
        return "{}({})".format(type(self).__name__, repr(self.path))

    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, key))
    
    
    def __getitem__(self, key):
        try:
            return self._data[key]
        except KeyError as e:
            pass
        
        if key == "names":
            val = set()
            for path in self.paths.get("names", ()):
                with open(path) as f:
                    for row in f:
                        name = row.strip()
                        if name:
                            val.add(name)
            if not val:
                sep = re.compile("[,;]")
                for target in self.targets:
                    for name in sep.split(target.name):
                        val.add(name.strip())
        
        elif key == "targets":
            if "targets" in self.paths:
                val = Gr(bed(self.paths["targets"]))
            else:
                val = self.exons
        
        elif key in ("transcripts", "codingregions", "exons", "codingexons"):
            val = Gr(gff3(self.paths["reference"], key, names=self["names"], principal=self.paths.get("principal", "")))
        
        else:
            raise KeyError(key)

        self._data[key] = val
        return val
    
    
    def __contains__(self, key):
        if key == "names":
            return "names" in self.paths or "targets" in self.paths
        
        elif key == "targets":
            return "targets" in self.paths or ("names" in self.paths and "reference" in self.paths)
        
        elif key in ("transcripts", "codingregions", "exons", "codingexons"):
            return "reference" in self.paths and "names" in self


    def get(self, key, *args):
        num_args = len(args) + 1
        if num_args > 2:
            raise TypeError(f"get expected at most 2 arguments, got {num_args}")
        elif num_args == 2:
            default = args[0]
        else:
            default = none
        
        try:
            return self[key]
        except KeyError:
            return default




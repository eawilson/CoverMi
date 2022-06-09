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
from collections import UserDict
from functools import partial

from .gr import Gr, bed, gff3, cosmic, gzopen

__all__ = ("Panel",)



REFSEQ_TRANSCRIPT = r"[NX][MR]_[0-9]+(\.[0-9]+)?"
ENSEMBL_TRANSCRIPT = r"ENST[0-9]{11}(\.[0-9]+)?"
GENE_SYMBOL = r"[A-Z][A-Z0-9orf_#@-]+"

GFF = "##gff-version 3"
BED = "chr[0-9a-zA-Z]+\t[0-9]+\t[0-9]+"
APPRIS_REFSEQ = f"{GENE_SYMBOL}\t[0-9]+\t{REFSEQ_TRANSCRIPT}.+\t(PRINCIPAL|ALTERNATIVE):"
APPRIS_ENSEMBL = f"{GENE_SYMBOL}\tENSG[0-9]+\t{ENSEMBL_TRANSCRIPT}.+\t(PRINCIPAL|ALTERNATIVE):"
COSMIC = "Gene name\tAccession Number\t"
GENE = f"{GENE_SYMBOL} *$"
GENE_REFSEQ_TRANSCRIPT = f"{GENE_SYMBOL} +{REFSEQ_TRANSCRIPT}[^\t]*$"
GENE_ENSEMBL_TRANSCRIPT = f"{GENE_SYMBOL} +{ENSEMBL_TRANSCRIPT}[^\t]*$"



REGEXPS = (("reference", re.compile(GFF)),
           ("amplicons",   re.compile(BED)),
           ("principal", re.compile(f"{APPRIS_REFSEQ}|{APPRIS_ENSEMBL}")),
           ("names",     re.compile(f"{GENE}|{GENE_REFSEQ_TRANSCRIPT}|{GENE_ENSEMBL_TRANSCRIPT}")),
           ("variants",  re.compile(COSMIC)),
          )



def identify(path):
    with gzopen(path, "rt") as f_in:
        try:
                # Don't get screwed by really big binary files
            contents = f_in.read(1000).splitlines()[:2]
        except UnicodeDecodeError:
            return []
    return [filetype for filetype, regexp in REGEXPS if any(regexp.match(row) for row in contents)]



class Panel(UserDict):    
    def __init__(self, *paths):
        super().__init__()
        
        self.paths = {}
        
        for path in paths:
            if os.path.isfile(path):
                self.add(path)
            elif os.path.isdir(path):
                for fn in os.listdir(path):
                    file_path = os.path.join(path, fn)
                    if os.path.isfile(file_path):
                        self.add(file_path)


    def add(self, path):
        filetypes = identify(path)
        
        if len(filetypes) == 1:
            filetype = filetypes.pop()
            self.paths[filetype] = os.path.abspath(path)
            self.clear()
            return filetype
        
        elif len(filetypes) > 1:
            raise RuntimeError(f"panel file {path} matches multiple file types")


    def __repr__(self):
        return "{}({})".format(type(self).__name__, ", ".join(sorted(self.paths.values())))


    def __missing__(self, key):
        val = None
        
        if key == "names":
            if "names" in self.paths:
                with open(self.paths["names"]) as f_in:
                    val = set(row.strip() for row in f_in if row.strip())
            
            elif "amplicons" in self:
                val = set(entry.name for entry in self["amplicons"])
        
        elif key == "amplicons":
            if "amplicons" in self.paths:
                val = Gr(bed(self.paths["amplicons"]))
        
        elif key == "targets":
            if "amplicons" in self:
                val = self["amplicons"]
            elif "exons" in self:
                val = self["exons"]
        
        elif key == "variants":
            if "variants" in self.paths:
                val = Gr(cosmic(self.paths["variants"]))
        
        else:
            all_genes = False
            if key.startswith("all"):
                all_genes = True
                key = key[3:]
            if key in ("transcripts", "codingregions", "exons", "codingexons"):
                if "reference" in self.paths and "names" in self:
                    val = Gr(gff3(self.paths["reference"], key, names=self["names"] if not all_genes else None, principal=self.paths.get("principal")))
        
        if val is None:
            raise KeyError(key)
        
        self[key] = val
        return val
    
    
    def __contains__(self, key):
        if key == "names":
            return "names" in self.paths or "amplicons" in self
        
        elif key == "amplicons":
            return "amplicons" in self.paths
        
        elif key == "targets":
            return "amplicons" in self or "exons" in self
        
        elif key == "variants":
            return "variants" in self.paths
        
        elif key in ("transcripts", "codingregions", "exons", "codingexons", "alltranscripts", "allcodingregions", "allexons", "allcodingexons"):
            return "reference" in self.paths and "names" in self
        
        return False

















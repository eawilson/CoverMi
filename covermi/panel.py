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
import os, re, pdb

from .gr import Gr, bed, reference

REFSEQ_TRANSCRIPT = r"[NX][MR]_[0-9]+(\.[0-9])?"
ENSEMBL_TRANSCRIPT = r"ENST[0-9]{11}"
GENE_SYMBOL = r"[A-Z][A-Zorf0-9-\.]+"

#Eleven columns - refflat
REGEXPS = (("reference", re.compile(GENE_SYMBOL+"\t"+REFSEQ_TRANSCRIPT+"\tchr.+?\t[+-]\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9]+\t[0-9,]+\t[0-9,]+$")), 
           ("reference", re.compile("#!genome-build")), #Ensembl gtf  
           ("targets",   re.compile("chr[0-9a-zA-Z]+\t[0-9]+\t[0-9]+")), # Bedfile
           ("principal", re.compile(GENE_SYMBOL+"\t[0-9]+\t"+REFSEQ_TRANSCRIPT+".+\t(PRINCIPAL|ALTERNATIVE):")),
           ("principal", re.compile(GENE_SYMBOL+"\tENSG[0-9]+\t"+ENSEMBL_TRANSCRIPT+".+\t(PRINCIPAL|ALTERNATIVE):")),
           ("names",     re.compile(GENE_SYMBOL+"$")), #Single column
           ("names",     re.compile(GENE_SYMBOL+" +"+REFSEQ_TRANSCRIPT+"$")), #Single column
           ("names",     re.compile(GENE_SYMBOL+" +"+ENSEMBL_TRANSCRIPT+"$")), #Single column
          )



class Panel(object):
    def __init__(self, panel_path):
        self.path = panel_path
        self.paths = {}
        self._data = {}
        for fn in os.listdir(panel_path):
            path = os.path.join(panel_path, fn)
            if os.path.isfile(path):
                filetype = None
                loops = 0
                with open(path, "rt") as f:
                    try:
                        for testrow in f.read(1000).split("\n"): # Don't get screwed by really big binary files
                            testrow = testrow.strip()
                            loops += 1
                            for thisfiletype, regexp in REGEXPS:
                                match = regexp.match(testrow)
                                if match:
                                    if match and filetype not in (thisfiletype, None):
                                        sys.exit("ERROR. File {} matches both {} and {} format".format(fn, thisfiletype, filetype))
                                    filetype = thisfiletype
                            if loops == 2 or match:
                                break
                    except UnicodeDecodeError:
                        pass
                if filetype is not None:
                    if filetype in self.paths:
                        sys.exit("ERROR. in {} panel - {} and {} are both of {} type".format(self.name, os.path.basename(self.paths[filetype]), fn, filetype))
                    self.paths[filetype] = os.path.abspath(path)

        if not self.paths:
            sys.exit("ERROR. {} panel is empty".format(self.name))
    
    
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
        
        try:
            if key == "names":
                try:
                    with open(self.paths["names"], "rt") as f:
                        val = set(f.read().splitlines())
                except KeyError:
                    val = self["targets"].names
            
            elif key == "targets":
                try:
                    val = Gr(bed(self.paths["targets"]))
                except KeyError as e:
                    val = self.exons
            
            elif key in ("transcripts", "codingregions", "exons", "codingexons"):
                val = Gr(reference(self.paths["reference"], key, names=self["names"], principal=self.paths.get("principal", "")))

            self._data[key] = val
            return val
        except Exception as e:
            raise e from None
    
    
    def __contains__(self, key):
        if key == "names":
            return "names" in self.paths or "targets" in self.paths
        
        elif key == "targets":
            return "targets" in self.paths or ("names" in self.paths and "reference" in self.paths)
        
        elif key in ("transcripts", "codingregions", "exons", "codingexons"):
            return "reference" in self.paths and "names" in self
        
        return False


    def get(self, key):
        return self[key]




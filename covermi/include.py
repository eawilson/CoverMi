from __future__ import print_function, absolute_import, division
import sys

class WarnOnce(object):
    def __init__(self):
        self.dedup = set()

    def __call__(self, text):
        if text not in self.dedup:
            self.dedup.add(text)
            print(text, file=sys.stderr)

class CoverMiException(Exception):
    pass

REFSEQ_TRANSCRIPT = r"[NX][MR]_[0-9]+(\.[0-9])?"
ENSEMBL_TRANSCRIPT = r"ENST[0-9]{11}"
GENE_SYMBOL = r"[A-Z][A-Zorf0-9-\.]+"

DEFAULT_ASSEMBLY = "GRCh37"
DEFAULT_TRANSCRIPTS = "refseq"
DEFAULT_SPECIES = "homo_sapiens"

__version__ = "3.1.5"

eprint = WarnOnce()

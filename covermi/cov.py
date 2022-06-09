import pdb
import struct
import sys
import gzip
from io import BufferedReader
from collections import Counter, defaultdict
import collections.abc

from .gr import Gr, Entry



MAX_LEN = 2**29-1 # As defined by max size supported by bai indexes.



def bisect_left(a, x):
    lo = 0
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid][1] < x: lo = mid+1
        else: hi = mid
    return lo



def zero_coverage():
    return [[0, MAX_LEN+1, 0]]



def multiply(gr, n):
    for entry in gr:
        for i in range(n):
            yield entry



def perfect_coverage(gr, n=1000):
    return Cov(reads=multiply(gr, n))



class Cov(collections.abc.Mapping):
    def __init__(self, path=None, reads=None):
        chrom_depths = defaultdict(Counter)
        
        if path:
            with Bam(path) as bam:
                for row in bam.coverage(): # Is zero/one based counting correct here
                    depths = chrom_depths[row[0]]
                    depths[row[1]] += 1
                    depths[row[2]+1] -= 1
        
        for entry in (reads or ()): # Is zero/one based counting correct here
            depths = chrom_depths[entry.chrom]
            depths[entry.start] += 1
            depths[entry.stop+1] -= 1
            
        self._data = defaultdict(zero_coverage)
        for chrom, depths in chrom_depths.items():
            start = 0
            depth = 0
            cov = []
            for pos, change in sorted(depths.items()):
                cov.append((start, pos-1, depth))
                start = pos
                depth += change
            if depth != 0: # Should never ever happen under any circumstance.
                raise RuntimeError("Fatal error in initial coverage calculation algorithm")
            cov.append((start, MAX_LEN+1, depth))
            self._data[chrom] = cov
    
    
    def __getitem__(self, key):
        return self._data[key]
    
    
    def __iter__(self):
        return iter(self._data)
    
    def __len__(self):
        return len(self._data)


    def calculate(self, region, depth, name=""):
        total = bool(name)
        results = defaultdict(CoverageInfo)
        for entry in region:
            if entry.stop >= entry.start:
                start = entry.start
                stop = entry.stop
            else: # if an insertion then calculate from base before to base afterwards
                start = entry.stop
                stop = entry.start
            if not total:
                name = entry.name
            info = results[name]
            info.name = name

            allcovered = True
            cchrom = self._data[entry.chrom]
            #              leftmost element where coverage.stop >= start
            for i in range(bisect_left(cchrom, start), len(cchrom)):
                cstart, cstop, cdepth = cchrom[i]
                if cstart > stop:
                    break
                elif cstop >= start:
                    bases = min(stop, cstop) - max(start, cstart) + 1
                    if cdepth >= depth:
                        info.bases_covered += bases
                        info.depth_covered += bases * cdepth
                        info.range_covered.add(Entry(entry.chrom, max(start, cstart), min(stop, cstop), entry.name, entry.strand))
                    else:
                        info.bases_uncovered += bases
                        info.depth_uncovered += bases * cdepth
                        info.range_uncovered.add(Entry(entry.chrom, max(start, cstart), min(stop, cstop), entry.name, entry.strand))
                        allcovered = False
            if allcovered:
                info.components_covered += 1
            else:
                info.components_uncovered += 1

        results = sorted(results.values())
        for info in results:
            info.depth_covered = info.depth_covered // max(info.bases_covered, 1)
            info.depth_uncovered = info.depth_uncovered // max(info.bases_uncovered, 1)
            info.range_covered = info.range_covered.merged
            info.range_uncovered = info.range_uncovered.merged
        return results[0] if total else results



class CoverageInfo(object):
    def __repr__(self):
        return "{} {}%".format(self.name, self.percent_covered)

    def __init__(self):
        self.name = ""
        self.depth_covered = 0
        self.depth_uncovered = 0
        self.bases_covered = 0
        self.bases_uncovered = 0
        self.range_covered = Gr()
        self.range_uncovered = Gr()
        self.components_covered = 0
        self.components_uncovered = 0

    def __lt__(self, other):
        return self.name < other.name
    
    @property
    def depth(self):
        return ((self.depth_covered * self.bases_covered) + (self.depth_uncovered * self.bases_uncovered)) // (self.bases or 1)

    @property
    def percent_covered(self):
        return float(self.bases_covered*100) / (self.bases or 1)

    @property
    def percent_uncovered(self):
        return 100 - self.percent_covered

    @property
    def range(self):
        return self.range_covered.combined_with(self.range_uncovered).merged

    @property
    def bases(self):
        return self.bases_covered + self.bases_uncovered

    @property
    def components(self):
        return self.components_covered + self.components_uncovered

    @property
    def percent_components_covered(self):
        return float(self.components_covered*100) / (self.components or 1)

    @property
    def percent_components_uncovered(self):
        return 100 - self.percent_components_covered

    @property
    def completely_covered(self):
        return not(self.incompletely_covered)

    @property
    def incompletely_covered(self):
        return bool(self.bases_uncovered)



class Bam(object):
    def __init__(self, path):
        """ May raise IOError or RuntimeError
        """
        self.path = path
        self.bam = None
        self.bai = None

        try:
            self.bam = BufferedReader(gzip.open(path, "rb"))
            if self.bam.read(4) != b"BAM\1":
                self.bam.close()
                raise RuntimeError(f"{path} is not a BAM file!")
            
            len_header_text = struct.unpack("<i", self.bam.read(4))[0]
            header_text = self.bam.read(len_header_text)
            num_ref_seq = struct.unpack("<i", self.bam.read(4))[0]
            chr_2_ref = {}
            self.ref_2_chr = [None] * num_ref_seq
            for x in range(0, num_ref_seq):
                len_ref_name = struct.unpack("<i", self.bam.read(4))[0]
                ref_name = self.bam.read(len_ref_name - 1)
                chrom = ref_name.decode("utf-8")
                self.ref_2_chr[x] = chrom
                chr_2_ref[chrom] = x
                self.bam.read(5)
            
        except struct.error:
            raise RuntimeError(f"{path} has a truncated header")
    
    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        self.close()


    def close(self):
        self.bam.close()


    def coverage(self):
        duplicate = 0x400
        secondary = 0x100
        unmapped = 0x4
        bad = duplicate | secondary | unmapped
        try:
            while True:
                read = self.bam.read(36)
                if len(read) == 0:
                    break
                
                block_size, ref_id, pos, bin_mq_nl, flag_nc, len_seq, next_ref_id, next_pos, len_template = struct.unpack("<iiiIIiiii", read)
                flag = flag_nc >> 16#// 0x10000
                if (ref_id == -1) or (flag & bad):
                    self.bam.read(block_size-32)
                else:
                    len_read_name = bin_mq_nl & 0xFF
                    n_cigar_op = flag_nc & 0xFFFF
                    direction = "-" if flag & 0x10 else "+"
                    start = pos + 1

                    read_name = self.bam.read(len_read_name - 1)
                    self.bam.read(1)

                    cigar_bytes = n_cigar_op * 4
                    length = 0
                    for cigar in struct.unpack("<" + "I" * n_cigar_op, self.bam.read(cigar_bytes)):
                        cigar_op = cigar & 0xF
                        if cigar_op in (0, 2, 7, 8):
                            length += cigar // 0x10
                        elif cigar_op == 3: # skip an intron
                            if length:
                                yield (self.ref_2_chr[ref_id], start, start + length - 1, direction)
                            start += length + (cigar//0x10)
                            length = 0
                    if length:
                        yield (self.ref_2_chr[ref_id], start, start + length - 1, direction)

                    self.bam.read(block_size - 32 - len_read_name - cigar_bytes)

        except struct.error:
            raise RuntimeError("{} is truncated".format(self.path))

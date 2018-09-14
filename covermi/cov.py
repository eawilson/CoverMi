# { chr_number : [ [start, stop, depth], [start, stop, depth], ... ] }
from __future__ import print_function, absolute_import, division

import csv, struct, os, sys, pdb
from collections import defaultdict, namedtuple
from io import BufferedReader
from itertools import islice, chain
from collections import Counter

from .gr import Gr, Chrom, Entry, MAX_LEN, PLUS, MINUS
from .include import CoverMiException, eprint

try: # python2
    from string import maketrans
    range = xrange
    PY3 = False
except ImportError: # python3
    maketrans = str.maketrans
    basestring = str
    PY3 = True

try:
    from Bio import bgzf
    BGZF = True
except ImportError:
    import gzip        
    BGZF = False


def bisect_left(a, x):
    lo = 0
    hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid][1] < x: lo = mid+1
        else: hi = mid
    return lo


def fake_paired_end_reads(amplicons, depth):
    for entry in amplicons.all_entries:
        third = (entry.stop - entry.start + 1) // 3
        for n in range(0, depth):
            yield (entry.chrom, entry.start, entry.stop-third, MINUS if entry.strand==MINUS else PLUS)
            yield (entry.chrom, entry.start+third, entry.stop, PLUS if entry.strand==MINUS else MINUS)


def extend(data, item):
    if len(data) and item[2] == data[-1][2]:
        data[-1][1] = item[1]
    else:
        data += [item]


def normalised(data, depth):
    newdata = defaultdict(zero_coverage)
    for key, values in data.items():
        newchrom = []
        for entry in values:
            extend(newchrom, [entry[0], entry[1], entry[2]>=depth])
        newdata[key] = newchrom
    return newdata


def combined(mydata, otherdata):
    newdata = defaultdict(zero_coverage)
    for key in set(mydata.keys()+otherdata.keys()):
         newchrom = []
         otherchrom = otherdata[key].__iter__()
         otherentry = [None, -1, None]
         for myentry in mydata[key]:
             while myentry[0] > otherentry[1]:
                 otherentry = next(otherchrom)
             while otherentry[1] < myentry[1]:
                 extend(newchrom, [max(myentry[0], otherentry[0]), otherentry[1], myentry[2]+otherentry[2]])
                 otherentry = next(otherchrom)
             extend(newchrom,[max(myentry[0], otherentry[0]), myentry[1], myentry[2]+otherentry[2]])
         newdata[key] = newchrom
    return newdata


class CumCov(object):

    def __init__(self, depth):
        self._data = defaultdict(zero_coverage)
        self.number = 0
        self.depth = depth


    def add(self, other):
        self.number += 1
        self._data = combined(self._data, normalised(other.data, self.depth))
        return self


    @property
    def data(self):
        newdata = defaultdict(zero_coverage)
        for key, chrom in self._data.items():
            newdata[key] = [(entry[0], entry[1], entry[2]*100//self.number) for entry in chrom]
        return newdata


def zero_coverage(): return [[0, MAX_LEN+1, 0]]


class AmpliconInfoList(list):
    def failure_rate(self, depth):
        return sum(i.failed(depth) for i in self)*100.0/len(self)


class Cov(object):

    def __init__(self, data=None, amplicons=(), locations=(), print_progress=False, q30=False, gc=False):
        self.data = defaultdict(zero_coverage)

        if data is None:
            return

        if isinstance(data, basestring):
            bampath = data
            bam = Bam(bampath, index=bool(locations))
            data = bam.coverage(locations=locations, print_progress=print_progress)
        else:
            bam = None

        fr_depth = {}
        self.amplicon_info = AmpliconInfoList()
        self.ontarget = 0
        self.offtarget = 0
        for amplicon in amplicons:
            amplicon_info = AmpliconInfo(amplicon)
            self.amplicon_info += [amplicon_info]
            fr_depth[(amplicon.start, PLUS)] = amplicon_info
            fr_depth[(amplicon.stop, MINUS)] = amplicon_info

        for chrom, start, stop, strand in data:
            depth = 1
            cov = self.data[chrom]
            upperx = -1
            lowerx = -1
            for x in range(len(cov) - 1, -1, -1):
                if stop < cov[x][0]:
                    continue
                if stop >= cov[x][1] and start <= cov[x][0]:
                    cov[x][2] += depth
                    if start == cov[x][0]:
                        break
                else:
                    if stop < cov[x][1]:
                        upperx = x
                    if start > cov[x][0]:
                        lowerx = x
                        break
            if upperx != -1:
                cov.insert(upperx, [cov[upperx][0], stop, cov[upperx][2] + depth])
                cov[upperx+1][0] = stop + 1
            if lowerx != -1:
                cov.insert(lowerx+1, [start, cov[lowerx][1], cov[lowerx][2] + depth * (upperx!=lowerx)])
                cov[lowerx][1] = start - 1
                if upperx == lowerx:
                    cov[lowerx][2] -= depth

            if amplicons:
                try:
                    fr_depth[(start, PLUS)].f_depth += 1
                    self.ontarget += 1
                except KeyError:
                    try:
                        fr_depth[(stop, MINUS)].r_depth += 1
                        self.ontarget += 1
                    except KeyError:
                        self.offtarget += 1                        

        self.amplicon_info.sort()
        if bam:
            self.unmapped = bam.unmapped
            self.mapped = bam.mapped
            self.reads = self.mapped + self.unmapped
            self.percent_unmapped = self.unmapped *100.0 / (self.reads or 1)
        if amplicons:
            self.percent_offtarget = self.offtarget * 100.0 / (self.ontarget + self.offtarget or 1)

        if q30 or gc:
            if bam:
                bam.close()
                bam = Bam(bampath, index=bool(locations))
                q30_tot = 0
                total = 0
                gc_tot = 0
                at_tot = 0

                for read in bam.read():
                    if q30:
                        quality = read.unclipped_quality
                        q30_tot += len([qual for qual in quality if qual >= 30])
                        total += len(quality)
                    if gc:
                        count = Counter(read.unclipped_sequence)
                        gc_tot += count["G"] + count["C"]
                        at_tot += count["A"] + count["T"]
                if q30:
                    self.q30 = q30_tot * 100.0 / total if total>0 else 0
                if gc:
                    self.gc = float(gc_tot) / ((at_tot+gc_tot) or 1)
        if bam: bam.close()


#    def load(self, path):
#        self.__init__(data=())
#        with open(path, "rb") as f:
#            for line in csv.reader(f, delimiter="\t"):
#                self.data[line[0]].append([int(value) for value in line[1:4]])


#    def save(self, path):
#        with open(path, "wb") as f:
#            writer = csv.writer(f, delimiter="\t")
#            for chrom in self.data:
#                for entry in self.data[chrom]:
#                    writer.writerow([chrom]+entry)


#    def save_range(self, path, gr): # Used to dump a standard region to disk for diagnostic purposes
#        with open(path, "wb") as f:
#            writer = csv.writer(f, delimiter="\t")
#            for entry in gr:
#                for start, stop, depth in self.data[entry.chrom]:
#                    if stop >= entry.start:
#                        if start > entry.stop:
#                            break
#                        writer.writerow([entry.chrom, max(start, entry.start), min(stop, entry.stop), depth])


    def as_range(self, min_depth=1):
        gr = Gr()
        for chrom, data in self.data.items():
            start = None
            for cstart, cstop, depth in data:
                if depth < min_depth and start is not None:
                    gr.add(Entry(chrom, start, stop))
                    start = None
                elif depth >= min_depth:
                    if start is None:
                        start = cstart
                    stop = cstop
        gr.sort()
        return gr


    def calculate(self, gr1, min_depth, exons=False, total=False):
        if isinstance(gr1, Entry): gr1 = Gr(gr1)
        results = defaultdict(CoverageInfo)
        for entry in gr1.all_entries:
            if entry.stop >= entry.start: # if an insertion then calculate from base before to base afterwards
                start = entry.start
                stop = entry.stop
            else:
                start = entry.stop
                stop = entry.start
            name = entry.name if (not total) else "Total"
            name = name if (not exons) else "{0} e{1}".format(name, entry.exon)
            info = results[name]
            info.name = name
            info.diseases = getattr(entry, "diseases", "")

            allcovered = True
            cchrom = self.data[entry.chrom]
            bisect = bisect_left(cchrom, start) # leftmost coverage.stop >= entry.start
            for cstart, cstop, cdepth in cchrom[bisect:]:
                if cstart > stop:
                    break
                elif cstop >= start:
                    bases = min(stop, cstop) - max(start, cstart) + 1
                    if cdepth>=min_depth:
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
                info.weighted_components_covered += getattr(entry, "weight", 1)
            else:
                info.components_uncovered += 1
                info.weighted_components_uncovered += getattr(entry, "weight", 1)

        results = sorted(results.values())
        for info in results:
            info.depth_covered = info.depth_covered // max(info.bases_covered, 1)
            info.depth_uncovered = info.depth_covered // max(info.bases_uncovered, 1)
            info.range_covered = info.range_covered.merged
            info.range_uncovered = info.range_uncovered.merged
        return results if (not total) else results[0]


class AmpliconInfo(object):

    def __repr__(self):
        return "{} forward={} reverse={}".format(self.entry, self.f_depth, self.r_depth)

    def __init__(self, entry):
        self.entry = entry
        self.f_depth = 0
        self.r_depth = 0

    def __lt__(self, other):
        return self.entry < other.entry

    @property
    def name(self):
        return self.entry.name

    @property
    def chrom(self):
        return self.entry.chrom

    @property
    def min_depth(self):
        return min(self.f_depth, self.r_depth)

    @property
    def max_depth(self):
        return max(self.f_depth, self.r_depth)

    @property
    def mean_depth(self):
        return (self.f_depth + self.r_depth) // 2

    @property
    def ratio(self):
        return float(self.min_depth) // max(self.max_depth, 1)

    @property
    def gr(self):
        return Gr(self.entry)

    def failed(self, depth):
        return (self.f_depth < depth) or (self.r_depth < depth)


class CoverageInfo(object):

    def __repr__(self):
        return "{} {}%".format(self.name, self.percent_covered)

    def __init__(self):
        self.name = ""
        self.diseases = ""
        self.depth_covered = 0
        self.depth_uncovered = 0
        self.bases_covered = 0
        self.bases_uncovered = 0
        self.range_covered = Gr()
        self.range_uncovered = Gr()
        self.components_covered = 0
        self.components_uncovered = 0
        self.weighted_components_covered = 0
        self.weighted_components_uncovered = 0

    def __lt__(self, other):
        return self.name < other.name
    
    @property
    def depth(self):
        return ((self.depth_covered * self.bases_covered) + (self.depth_uncovered + self.bases_uncovered)) // (self.bases or 1)

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
    def weighted_components(self):
        return self.weighted_components_covered + self.weighted_components_uncovered

    @property
    def percent_weighted_components_covered(self):
        return float(self.weighted_components_covered*100) / (self.weighted_components or 1)

    @property
    def percent_weighted_components_uncovered(self):
        return 100 - self.percent_weighted_components_covered

    @property
    def completely_covered(self):
        return not(self.incompletely_covered)

    @property
    def incompletely_covered(self):
        return bool(self.bases_uncovered)





def reg2bins(beg, end): # accepts 1 based start and stop
    beg -= 1
    end -= 2
    bins = [0]

    for k in range(1 + (beg>>26),  2 + (end>>26)): bins += [k]
    for k in range(9 + (beg>>23), 10 + (end>>23)): bins += [k]
    for k in range(73 + (beg>>20), 74 + (end>>20)): bins += [k]
    for k in range(585 + (beg>>17), 586 + (end>>17)): bins += [k]
    for k in range(4681 + (beg>>14), 4682 + (end>>14)): bins += [k]
    return bins


def reg2bin(beg, end): # accepts 1 based start and stop
    beg -= 1
    end -= 1
    if ((beg>>14) == (end>>14)):
        return ((1<<15)-1)//7 + (beg>>14)
    if ((beg>>17) == (end>>17)):
        return ((1<<12)-1)//7 + (beg>>17)
    if ((beg>>20) == (end>>20)):
        return ((1<<9 )-1)//7 + (beg>>20)
    if ((beg>>23) == (end>>23)):
        return ((1<<6 )-1)//7 + (beg>>23)
    if ((beg>>26) == (end>>26)):
        return ((1<<3 )-1)//7 + (beg>>26)
    return 0


class Bai(object):

    def __init__(self, fn, chr_2_ref):
        try:
            self.chr_2_ref = chr_2_ref
            with open(fn, "rb") as bai:
                if bai.read(4) != b"BAI\1":
                    raise CoverMiException("{} is not a BAI file!".format(os.path.basename(fn)))
                
                n_ref = struct.unpack("<i", bai.read(4))[0] # number of reference sequences
                self.bins = [defaultdict(tuple) for x in range(0, n_ref)]
                self.intervals = [() for x in range(0, n_ref)]
                for ref_id in range(0, n_ref):

                    n_bin = struct.unpack("<i", bai.read(4))[0] # number of bins for that reference sequence
                    for distinct_bin in range(0, n_bin):
                    
                        bin, n_chunk = struct.unpack("<Ii", bai.read(8)) # bin_number, number of chunks within that bin
                        if n_chunk:
                            if bin == 37450: # pseudo bin for unmapped reads
                                bai.read(32)
                            else:
                                chunks = [None] * n_chunk
                                for chunk in range(0, n_chunk):
                                    chunks[chunk] = struct.unpack("<QQ", bai.read(16))
                                self.bins[ref_id][bin] = chunks

                    n_intv = struct.unpack("<i", bai.read(4))[0] # length of linear index
                    if n_intv:
                        intervals = [None] * n_intv
                        for intv in range(0, n_intv):
                            intervals[intv] = struct.unpack("<Q", bai.read(8))[0]
                        self.intervals[ref_id] = intervals

        except IOError:
            raise CoverMiException("IOError reading {}".format(os.path.basename(fn)))

        except struct.error:
            raise CoverMiException("{} is truncated".format(os.path.basename(fn)))
        


    def chunks(self, locations):
        chunks = set()
        for loc in locations:
            ref_id = self.chr_2_ref[loc.chrom]
            intervals = self.intervals[ref_id]
            interval = (loc.start-1)//16384
            #lower_bound = intervals[interval] if interval < len(intervals) else intervals[-1]
            for bin in reg2bins(loc.start, loc.stop):
                for vstartstop in self.bins[ref_id][bin]:
                    #if vstartstop[1] >= lower_bound:
                        chunks.add(vstartstop)
        return sorted(chunks)


BamStats = namedtuple("BamStats", "mapq unmapped duplicate")
byte2base = "=ACMGRSVTWYHKDBN"

SWAPSTRAND = maketrans("ATGC", "TACG")
def invert(nucleotides):
    return nucleotides[::-1].translate(SWAPSTRAND)


Cigar = namedtuple("Cigar", "code length")


class BamFlags(object):
    __slots__ = ("flags",)

    def __repr__(self):
        return "Bamflags(flags={}, unmapped={}, reverse={})".format(bin(self.flags), self.unmapped, self.reverse_complement)

    def __init__(self, flags):
        self.flags = flags

    @property
    def multiple_segments(self):
        return bool(self.flags & 0x1)

    @property
    def all_segments_aligned(self):
        return bool(self.flags & 0x2)

    @property
    def unmapped(self):
        return bool(self.flags & 0x4)

    @property
    def next_segment_unmapped(self):
        return bool(self.flags & 0x8)

    @property
    def reverse_complement(self):
        return bool(self.flags & 0x10)

    @property
    def next_segment_reverse_complement(self):
        return bool(self.flags & 0x20)

    @property
    def first_segment(self):
        return bool(self.flags & 0x40)

    @property
    def last_segment(self):
        return bool(self.flags & 0x80)

    @property
    def secondary_alignment(self):
        return bool(self.flags & 0x100)

    @property
    def not_passing_filters(self):
        return bool(self.flags & 0x200)

    @property
    def not_passing_filters(self):
        return bool(self.flags & 0x200)

    @property
    def duplicate(self):
        return bool(self.flags & 0x400)

    @property
    def supplementary_alignment(self):
        return bool(self.flags & 0x800)


class BamRead(object):
    __slots__ = ("chrom", "pos", "flags", "len_read_name", "n_cigar_op", "len_seq", "data", "_cigar")

    def __repr__(self):
        return "BamRead(chrom={}, pos='{}', flags={}, cigar={}, sequence='{}', quality='{}')".format(self.chrom, self.pos, repr(self.flags), self.cigar,
                                                                                                           self.unclipped_sequence, self.unclipped_quality)

    def __init__(self, chrom, pos, flags, len_read_name, n_cigar_op, len_seq, data):
        self.chrom = chrom
        self.pos = pos
        self.flags = BamFlags(flags)
        self.len_read_name = len_read_name
        self.n_cigar_op = n_cigar_op
        self.len_seq = len_seq
        self.data = data

    @property
    def name(self):
        return self.data[:self.len_read_name-1]

    @property
    def strand(self):
        return "-" if self.flags.reverse_complement else "+"

    @property
    def cigar(self):
        try: return self._cigar
        except AttributeError: pass
        start = self.len_read_name
        stop = start + (self.n_cigar_op * 4)
        self._cigar = [Cigar(cigar & 0xF, cigar // 0x10) for cigar in struct.unpack("<" + "I" * self.n_cigar_op, self.data[start:stop])]
        return self._cigar

    @property
    def sequence(self):
        cigar = self.cigar
        if not cigar:
            start = stop = 0
        else:
            if cigar[0].code == 4: start = cigar[0].length
            elif cigar[0].code == 5 and cigar[1].code == 4: start = cigar[1].length
            else: start = 0
            if cigar[-1].code == 4: stop = cigar[-1].length
            elif cigar[-1].code == 5 and cigar[-2].code == 4: stop = cigar[-2].length
            else: stop = 0
        return self._sequence(startclip=start, stopclip=stop)

    @property
    def quality(self):
        cigar = self.cigar
        if not cigar:
            start = stop = 0
        else:
            if cigar[0].code == 4: start = cigar[0].length
            elif cigar[0].code == 5 and cigar[1].code == 4: start = cigar[1].length
            else: start = 0
            if cigar[-1].code == 4: stop = cigar[-1].length
            elif cigar[-1].code == 5 and cigar[-2].code == 4: stop = cigar[-2].length
            else: stop = 0
        return self._quality(startclip=start, stopclip=stop)

    @property
    def raw_sequence(self):
        return self._sequence() if self.strand == "+" else invert(self._sequence())

    @property
    def raw_quality(self):
        return self._quality() if self.strand == "+" else self._quality[-1::-1]

    @property
    def unclipped_sequence(self):
        return self._sequence()

    @property
    def unclipped_quality(self):
        return self._quality()

    def _sequence(self, startclip=0, stopclip=0):
        len_seq = self.len_seq - stopclip
        sequence_bytes = ((len_seq + 1) // 2)
        start = self.len_read_name + (self.n_cigar_op * 4)
        stop = start + sequence_bytes
        sequence_list = [(byte2base[byte // 0x10], byte2base[byte & 0xF]) for byte in struct.unpack("<" + "B" * sequence_bytes, self.data[start:stop])]
        return "".join(islice(chain(*sequence_list), startclip, len_seq))

    def _quality(self, startclip=0, stopclip=0):
        start = self.len_read_name + (self.n_cigar_op * 4) + ((self.len_seq + 1) // 2) + startclip
        stop = start + self.len_seq - stopclip
        return struct.unpack("<" + "B" * (stop - start), self.data[start:stop])

    @property
    def start(self):
        return self.pos

    @property
    def stop(self):
        length = 0
        for op in self.cigar:
            if op.code in (0, 2, 7, 8):
                length += op.length
        return self.start + length - 1


class Bam(object):

    COVERAGE = 1
    STATS = 2

    Q30 = 1
    READ_LENGTH = 2

    def __init__(self, fn, index=True):
        self.bam = None
        self.bai = None

        try:
            self.bam = bgzf.open(fn, "rb") if BGZF else BufferedReader(gzip.open(fn, "rb"))
            if self.bam.read(4) != b"BAM\1":
                self.bam.close()
                raise CoverMiException("{} is not a BAM file!".format(os.path.basename(fn)))
            
            len_header_text = struct.unpack("<i", self.bam.read(4))[0]
            header_text = self.bam.read(len_header_text)
            num_ref_seq = struct.unpack("<i", self.bam.read(4))[0]
            chr_2_ref = {}
            self.ref_2_chr = [None] * num_ref_seq
            for x in range(0, num_ref_seq):
                len_ref_name = struct.unpack("<i", self.bam.read(4))[0]
                ref_name = self.bam.read(len_ref_name - 1)
                if PY3:
                    ref_name = ref_name.decode("utf-8")
                try:
                    chrom = Chrom(ref_name)
                except KeyError:
                    chrom = ref_name
                self.ref_2_chr[x] = chrom
                chr_2_ref[chrom] = x
                self.bam.read(5)

        except IOError:
            if self.bam: self.bam.close()
            raise CoverMiException("IOError reading {}".format(os.path.basename(fn)))

        except (IOError, struct.error):
            if self.bam: self.bam.close()
            raise CoverMiException("{} has a truncated header".format(os.path.basename(fn)))
             
        if index:
            if BGZF:
                try:
                    self.bai = Bai(fn+".bai", chr_2_ref)
                    eprint("Using index file {}".format(os.path.basename(fn+".bai")))
                except CoverMiException as e:
                    eprint(str(e))
            else:
                eprint("Biopython not installed therefore cannot use bam indexes")

    def __enter__(self):
        return self


    def __exit__(self, type, value, traceback):
        self.close()


    def close(self):
        self.bam.close()


    def coverage(self, locations=(), print_progress=False):
        self.unmapped = 0
        self.mapped = 0
        if locations and self.bai:
            chunks = self.bai.chunks(locations)
        else:
            chunks = ((0, 0),)
            print_progress = False
        try:
            for chunk_no, chunk in enumerate(chunks):
                vstart, vstop = chunk
                if vstart:
                    self.bam.seek(vstart)
                if print_progress:
                    sys.stdout.write("Reading chunk {} of {}\r".format(chunk_no+1, len(chunks)))
                    sys.stdout.flush()
                while True:
                    if vstop and self.bam.tell() > vstop:
                        break
                    read = self.bam.read(36)
                    if len(read) == 0:
                        break
                    
                    block_size, ref_id, pos, bin_mq_nl, flag_nc, len_seq, next_ref_id, next_pos, len_template = struct.unpack("<iiiIIiiii", read)
                    flag = flag_nc >> 16#// 0x10000
                    unmapped = (ref_id == -1) or (flag & 0x4)
#                    duplicate = flag & 0x400
#                    secondary = flag & 0x100
#                    supplementary = flag & 0x800

                    if unmapped: # unmapped read
                        self.bam.read(block_size-32)
                        self.unmapped += 1
                    else:
                        self.mapped += 1
                        len_read_name = bin_mq_nl & 0xFF
                        n_cigar_op = flag_nc & 0xFFFF
                        direction = MINUS if flag & 0x10 else PLUS
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
            if print_progress:
                print()

        except (IOError, struct.error):
            self.bam.close()
            raise CoverMiException("{} is truncated".format(os.path.basename(fn)))


    def read(self, locations=()):
        try:
            for vstart, vstop in self.bai.chunks(locations) if locations and self.bai else ((0, 0),):
                if vstart:
                    self.bam.seek(vstart)
                while True:
                    if vstop and bam.tell() > vstop:
                        break
                    read = self.bam.read(36)
                    if len(read) == 0:
                        break
                    
                    block_size, ref_id, pos, bin_mq_nl, flag_nc, len_seq, next_ref_id, next_pos, len_template = struct.unpack("<iiiIIiiii", read)

                    yield BamRead(chrom=self.ref_2_chr[ref_id], pos=pos + 1, flags=flag_nc // 0x10000, len_read_name=bin_mq_nl & 0xFF, 
                                  n_cigar_op=flag_nc & 0xFFFF, len_seq=len_seq, data=self.bam.read(block_size - 32))

        except (IOError, struct.error):
            self.bam.close()
            raise CoverMiException("{} is truncated".format(fn))


    def info(self, what, locations=()):
        if what not in (Bam.Q30, Bam.READ_LENGTH):
            raise CoverMiException("Bam: Unknown option: {}".format(what))
        passing_q30 = 0
        total = 0
        try:
            for vstart, vstop in self.bai.chunks(locations) if locations and self.bai else ((0, 0),):
                if vstart:
                    self.bam.seek(vstart)
                while True:
                    if vstop and bam.tell() > vstop:
                        break
                    read = self.bam.read(36)
                    if len(read) == 0:
                        break
                    
                    block_size, ref_id, pos, bin_mq_nl, flag_nc, len_seq, next_ref_id, next_pos, len_template = struct.unpack("<iiiIIiiii", read)

                    flag = flag_nc // 0x10000
                    unmapped = (ref_id == -1) or flag & 0x4
                    duplicate = flag & 0x400
                    secondary = flag & 0x100
                    supplementary = flag & 0x800
                    n_cigar_op = flag_nc & 0xFFFF

                    if unmapped or duplicate or secondary or supplementary or n_cigar_op==0:
                        self.bam.read(block_size-32)
                        continue

                    len_read_name = bin_mq_nl & 0xFF
                    self.bam.read(len_read_name)

                    cigar_bytes = n_cigar_op * 4
                    cigar_string = [(cigar & 0xF, cigar // 0x10) for cigar in struct.unpack("<" + "I" * n_cigar_op, self.bam.read(cigar_bytes))]

                    start = 1 if cigar_string[0][0] == 5 else 0
                    if cigar_string[start][0] == 4:
                        soft_clipped = cigar_string[start][1]
                        start += 1
                    else:
                        soft_clipped = 0
                    length = sum([cigar[1] for cigar in cigar_string[start:] if cigar[0] in (0, 2, 7, 8)])
                    if what == Bam.READ_LENGTH: return length

                    self.bam.read((len_seq+1)//2)
                    quality = struct.unpack("<" + "B" * len_seq, self.bam.read(len_seq))[soft_clipped:length+soft_clipped]
                    for qual in quality:
                        if qual >= 30: passing_q30 += 1
                    total += len(quality)

                    self.bam.read(block_size - 32 - len_read_name - cigar_bytes - ((len_seq+1)/2) - len_seq)

        except (IOError, struct.error):
            self.bam.close()
            raise CoverMiException("{} is truncated".format(fn))

        return passing_q30 * 100.0 / total if total>0 else 0

    def sequence(self, locations=(), options=()):
        try:
            for vstart, vstop in self.bai.chunks(locations) if locations and self.bai else ((0, 0),):
                if vstart:
                    self.bam.seek(vstart)
                while True:
                    if vstop and bam.tell() > vstop:
                        break
                    read = self.bam.read(36)
                    if len(read) == 0:
                        break
                    
                    block_size, ref_id, pos, bin_mq_nl, flag_nc, len_seq, next_ref_id, next_pos, len_template = struct.unpack("<iiiIIiiii", read)

                    flag = flag_nc // 0x10000
                    unmapped = (ref_id == -1) or flag & 0x4
                    reverse = flag_nc & 0x10
                    duplicate = flag & 0x400
                    secondary = flag & 0x100
                    supplementary = flag & 0x800
                    n_cigar_op = flag_nc & 0xFFFF

#                    if unmapped or duplicate or sinvertecondary or supplementary or n_cigar_op==0:
#                        self.bam.read(block_size-32)invert
#                        continue

                    len_read_name = bin_mq_nl & 0xFF
                    self.bam.read(len_read_name)
                    cigar_bytes = n_cigar_op * 4
                    cigar_string = [(cigar & 0xF, cigar // 0x10) for cigar in struct.unpack("<" + "I" * n_cigar_op, self.bam.read(cigar_bytes))]

#                    start = 1 if cigar_string[0][0] == 5 else 0
#                    if cigar_string[start][0] == 4:
#                        soft_clipped = cigar_string[start][1]
#                        start += 1
#                    else:
#                        soft_clipped = 0
#                    length = sum([cigar[1] for cigar in cigar_string[start:] if cigar[0] in (0, 2, 7, 8)])
#                    if what == Bam.READ_LENGTH: return length
                    sequence_bytes = (len_seq+1)//2
                    sequence = struct.unpack("<" + "B" * sequence_bytes, self.bam.read(sequence_bytes))
                    sequence_list = []
                    for byte in sequence:
                        sequence_list += [byte2base[byte // 0x10], byte2base[byte & 0xF]]

                    sequence = "".join(sequence_list[:len_seq+1])
                    yield sequence if not reverse else invert(sequence)
                    quality = struct.unpack("<" + "B" * len_seq, self.bam.read(len_seq))#[soft_clipped:length+soft_clipped]
#                    for qual in quality:
#                        if qual >= 30: passing_q30 += 1
#                    total += len(quality)

                    self.bam.read(block_size - 32 - len_read_name - cigar_bytes - sequence_bytes - len_seq)

        except (IOError, struct.error):
            self.bam.close()
            raise CoverMiException("{} is truncated".format(fn))


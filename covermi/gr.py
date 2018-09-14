from __future__ import print_function, absolute_import, division

import csv, sys, re, os, pdb
from itertools import islice, repeat
from collections import defaultdict, Counter, namedtuple
from operator import attrgetter
from copy import copy
from gzip import GzipFile

try: # python2
    from itertools import imap as map
    from itertools import izip as zip
    from string import maketrans
    range = xrange
    CSV_READ_MODE = "rb"
    PLUS = "+"
    MINUS = "-"
    DOT = "."
except ImportError: # python3
    maketrans = str.maketrans
    basestring = str
    CSV_READ_MODE = "rt"
    PLUS = sys.intern("+")
    MINUS = sys.intern("-")
    DOT = sys.intern(".")
    
from .include import *



SPLICE_SITE_BUFFER = 5
MAX_LEN = 2**29-1 # As defined by max size supported by bai indexes.
HETERO_LL = .35
HETERO_UL = .65
HOMO_LL = .9

SWAPSTRAND = maketrans("ATGC", "TACG")
def invert(nucleotides):
    return nucleotides[::-1].translate(SWAPSTRAND)


class FileContext(object):
    def __init__(self, filename_or_iterable, mode="rt"):
        self.filename_or_iterable = filename_or_iterable
        self.mode = mode

    def __enter__(self):
        if isinstance(self.filename_or_iterable, basestring):
            self.fp = open(self.filename_or_iterable, self.mode)
            return self.fp
        return self.filename_or_iterable

    def __exit__(self, type_, value, traceback):
        if hasattr(self, "fp"):
            self.fp.close()


class Translate(dict):
    def __missing__(self, key):
        return key
    
    
class Chrom(int):
    _instances = {}

    def __new__(cls, val):
        val = STR_2_CHR[val]
        try:
            return cls._instances[val]
        except KeyError:
            cls._instances[val] = super(Chrom, cls).__new__(cls, val)
            return cls._instances[val]

    def __repr__(self):
        return "Chrom({})".format(self)

    def __str__(self):
        return CHR_2_STR[self] 

    def __add__(self, other):
        return CHR_2_STR[self]+other

    def __radd__(self, other):
        return other+CHR_2_STR[self]


STR_2_CHR =      dict(( ("chr1", 1), ("1", 1), (1, 1),
                        ("chr2", 2), ("2", 2), (2, 2),
                        ("chr3", 3), ("3", 3), (3, 3),
                        ("chr4", 4), ("4", 4), (4, 4),
                        ("chr5", 5), ("5", 5), (5, 5),
                        ("chr6", 6), ("6", 6), (6, 6),
                        ("chr7", 7), ("7", 7), (7, 7),
                        ("chr8", 8), ("8", 8), (8, 8),
                        ("chr9", 9), ("9", 9), (9, 9),
                        ("chr10", 10), ("10", 10), (10, 10),
                        ("chr11", 11), ("11", 11), (11, 11),
                        ("chr12", 12), ("12", 12), (12, 12),
                        ("chr13", 13), ("13", 13), (13, 13),
                        ("chr14", 14), ("14", 14), (14, 14),
                        ("chr15", 15), ("15", 15), (15, 15),
                        ("chr16", 16), ("16", 16), (16, 16),
                        ("chr17", 17), ("17", 17), (17, 17),
                        ("chr18", 18), ("18", 18), (18, 18),
                        ("chr19", 19), ("19", 19), (19, 19),
                        ("chr20", 20), ("20", 20), (20, 20),
                        ("chr21", 21), ("21", 21), (21, 21),
                        ("chr22", 22), ("22", 22), (22, 22),
                        ("chrX", 23), ("chr23", 23), ("23", 23), ("X", 23), (23, 23),
                        ("chrY", 24), ("chr24", 24), ("24", 24), ("Y", 24), (24, 24),
                        ("chrM", 25), ("chr25", 25), ("25", 25), ("M", 25), ("MT", 25), (25, 25),
                        ))

CHR_2_STR =      dict(( (1, "chr1"), (2, "chr2"), (3, "chr3"), (4, "chr4"), (5, "chr5"), (6, "chr6"), (7, "chr7"), (8, "chr8"), (9, "chr9"), (10, "chr10"), 
                        (11, "chr11"), (12, "chr12"), (13, "chr13"), (14, "chr14"), (15, "chr15"), (16, "chr16"), (17, "chr17"), (18, "chr18"), (19, "chr19"), (20, "chr20"), 
                        (21, "chr21"), (22, "chr22"), (23, "chrX"), (24, "chrY"), (25, "chrM"), 
                        ))

nucleotide = {"A": "R", "G": "R", "C": "Y", "T": "Y"}# puRine: A, G,  pYrimadine: T, C

def _properties(cls, self):
        properties = []
        for slot in cls.__slots___:
            try:
                properties += ["{}={}".format(slot, repr(getattr(self, slot)))]
            except AttributeError:
                pass
        return properties


class PickleSafe(object):

    def __getstate__(self):
        data = {slot: getattr(self, slot) for slot in dir(self) if not (slot.startswith("__") or hasattr(type(self), slot))}
        return data

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)


class Entry(object):
    __slots___ = ("chrom", "start", "stop", "name", "strand")

    def __repr__(self):
        return "{}({})".format(type(self).__name__, ", ".join(self._properties()))

    def _properties(self):
        return list(map(lambda attr:repr(getattr(self, attr)), Entry.__slots___))

    def __str__(self):
        return self.location

    def __init__(self, chrom, start, stop, name=DOT, strand=DOT):
        self.chrom = Chrom(chrom)
        self.start = start
        self.stop = stop
        self.name = DOT if name==DOT else name
        if strand == PLUS:
            self.strand = PLUS
        elif strand == MINUS:
            self.strand = MINUS
        else:
            self.strand = DOT

    @property
    def identifier(self):
        return (self.chrom, self.start, self.stop)

    @property
    def location(self):
        return "{}:{}-{}".format(self.chrom, self.start, self.stop)

    def __eq__(self, other):
        try:
            return self.identifier == other.identifier
        except AttributeError:
            return self.identifier == other

    def __hash__(self):
        return hash(self.identifier)

    def __lt__(self, other):
        return self.identifier < other.identifier


class Transcript(Entry):
    __slots__ = ("transcript", "gene")

    def __init__(self, chrom, start, stop, name, strand, gene, transcript):
        super(Transcript, self).__init__(chrom, start, stop, name, strand)
        self.gene = gene
        self.transcript = transcript

    def _properties(self):
        return super(Transcript, self)._properties() + _properties(Transcript, self)


class Exon(Transcript):
    __slots___ = ("exon",)

    def __init__(self, chrom, start, stop, name, strand, gene, transcript, exon):
        super(Exon, self).__init__(chrom, start, stop, name, strand, gene, transcript)
        self.exon = exon

    def _properties(self):
        return super(Exon, self)._properties() + _properties(Exon, self)


class Variant(Entry):
    __slots___ = ("ref", "alt", "prefix")

    def _properties(self):
        return super(Variant, self)._properties() + _properties(Variant, self)

    def __str__(self):
        return "{}:{} {}/{}".format(self.chrom, self.start, self.ref, self.alt)

    def __init__(self, chrom, pos, ref, alt, name=DOT, strand=DOT):
        if ref == alt:
            raise RuntimeError("Alt allele cannot be equal to ref allele")
        pos = int(pos)
        prefix = ""
        while ref[:1] == alt[:1]:
            prefix += ref[:1]
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
        while ref[-1:] == alt[-1:]:
            ref = ref[:-1]
            alt = alt[:-1]

        self.ref = ref or MINUS
        self.alt = alt or MINUS
        self.prefix = prefix
        super(Variant, self).__init__(chrom, pos, pos + (0 if self.ref==MINUS else len(ref)) - 1, name, strand)

    @property
    def pos(self):
        return self.start

    @property
    def identifier(self):
        return (self.chrom, self.pos, self.ref, self.alt)


    @property
    def vartype(self):
        if self.ref == MINUS:
            return "ins"
        elif self.alt == MINUS:
            return "del"
        elif len(self.ref) == len(self.alt) == 1:
            return "snp"
        else:
            return "delins"

    @property
    def substitution(self):
        try:
            return "transition" if nucleotide[self.ref]==nucleotide[self.alt] else "transversion"
        except KeyError:
            return None

    @property
    def as_vcf(self):
        return (self.chrom, self.pos - len(self.prefix), DOT, self.prefix + self.ref, self.prefix + self.alt, DOT, DOT, DOT, DOT)



class SequencedVariant(Variant):
    __slots___ = ("qual", "depth", "alt_depth", "filters", "zygosity")

    def _properties(self):
        return super(SequencedVariant, self)._properties() + _properties(SequencedVariant, self)

    @property
    def vaf(self):
        return float(self.alt_depth)/self.depth

    @property
    def as_vcf(self):
        try:
            depth = "{},{}".format(self.depth, self.alt_depth)
        except AttributeError:
            depth = DOT
        return (self.chrom, self.pos - len(self.prefix), DOT, self.prefix + self.ref, self.prefix + self.alt, getattr(self, "qual", DOT), getattr(self, "filters", DOT), DOT,  depth)


class HgmdVariant(Variant):
    __slots___ = ("diseases", "weight", "hgvsc", "transcript", "gene", "disease")

    def _properties(self):
        return super(HgmdVariant, self)._properties() + _properties(HgmdVariant, self)


####################################################################################################
# Gr                                                                                               #
####################################################################################################

def bisect_left(a, start, lo, hi):
    """Return the index where to insert item x in list a, assuming a is sorted.
    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    x = start - a.maxlength + 1
    while lo < hi:
        mid = (lo+hi)//2
        if a[mid].start < x: lo = mid+1
        else: hi = mid
    return lo


class Iterate(object):
    def __init__(self, a):
        self.a = a
        self.lena = len(a)
        self.lo = 0

    def yield_overlapping(self, start, stop):
        #find_ge(a, x, lo, hi) 'Find leftmost item greater than or equal to x'
        self.lo = bisect_left(self.a, start, self.lo, self.lena)
        for index in range(self.lo, self.lena):
            entry = self.a[index]
            if entry.start > stop:
                break
            if entry.stop >= start:
                yield entry


class GrList(list, PickleSafe):
    __slots___ = ("maxlength", "maxstop")

    def __init__(self):
        super(GrList, self).__init__()
        self.maxlength = 0
        self.maxstop = 0
        

class Gr(object):

    def __repr__(self):
        return "Gr([{}])".format(", ".join(entry.__repr__() for entry in self))


    def __init__(self, data=None):
        self.data = defaultdict(GrList)
        self._tuple = None
        self._hash = None
        if data is not None:
            if hasattr(data, "chrom"):
                data = (data,)
            for entry in data:
                self.add(entry)
            self.sort()


    def add(self, entry):
        grlist = self.data[entry.chrom]
        grlist.append(entry)
        length = entry.stop - entry.start + 1
        if length > grlist.maxlength:
            grlist.maxlength = length
        if entry.stop > grlist.maxstop:
            grlist.maxstop = entry.stop
        return self


    def sort(self):
        for chrom in self.data.values():
            chrom.sort(key=attrgetter("start", "stop", "name"))


    @property
    def sorted_by_name(self):
        for entry in sorted([entry for entry in self], key=attrgetter("name", "start", "stop")):
            yield entry


    @property
    def all_entries(self):
        return self.__iter__()


    def __iter__(self):
        for key, chrom in sorted(self.data.items()):
            for entry in chrom:
                yield entry


    @property
    def as_tuple(self):
        if not self._tuple: self._tuple = tuple(entry.identifier for entry in self)
        return self._tuple


    def __hash__(self):
        if not self._hash: self._hash = hash(self.as_tuple)
        return self._hash


    def __eq__(self, other):
        if isinstance(other, Entry): other = Gr(other)
        return self.as_tuple == other.as_tuple


    @property
    def inverted(self):
        new = type(self)()
        for key in sorted(CHR_2_STR.keys()):
            last_stop = 0
            if key in self.data:
                for entry in self.data[key]:
                    if last_stop+1 < entry.start-1:
                        new.add(Entry(key, last_stop+1, entry.start-1))
                    if entry.stop > last_stop:
                        last_stop = entry.stop
            if last_stop+1 < MAX_LEN:
                new.add(Entry(key, last_stop+1, MAX_LEN))
        return new


    @property
    def merged(self): 
        new = type(self)()
        for key, chrom in self.data.items():
            new.add(chrom[0])
            nchrom = new.data[key]
            nchrom.maxstop = chrom.maxstop
            for entry in islice(chrom, 1, len(chrom)):
                if entry.start-1 <= nchrom[-1].stop:
                    if entry.stop > nchrom[-1].stop:
                        nchrom[-1] = Entry(key, nchrom[-1].start, entry.stop)
                        length = entry.stop - nchrom[-1].start + 1
                        if length > nchrom.maxlength:
                            nchrom.maxlength = length
                else:
                    nchrom.append(entry)
        return new


    def overlapped_by(self, other):
        if isinstance(other, Entry): other = Gr(other)

        def a_overlapped_by_b(a, b):
            if b.start <= a.start and b.stop >= a.stop:
                entry = a
            else:
                entry = copy(a)
                if b.start > entry.start:
                    entry.start = b.start
                if b.stop < entry.stop:
                    entry.stop = b.stop
            return entry
    
        new = type(self)()
        for key, chrom in self.data.items():
            if key in other.data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other.data[key])
                for a in iterateself.yield_overlapping(other.data[key][0].start, other.data[key].maxstop):
                    entry = None
                    for b in iterateother.yield_overlapping(a.start, a.stop):
                        if entry is None:
                            entry = a_overlapped_by_b(a, b)
                        else:
                            if b.start <= entry.stop + 1:
                                if b.stop > entry.stop:
                                    if b.stop < a.stop:
                                        entry.stop = b.stop
                                        continue
                                    else:
                                        entry.stop = a.stop
                                        break
                            else:
                                new.add(entry)
                                entry = a_overlapped_by_b(a, b)
                        if entry.stop == a.stop:
                            break
                    if entry is not None:
                        new.add(entry)
        new.sort()
        return new


    def not_touched_by(self, other):
        if isinstance(other, Entry): other = Gr(other)
        new = type(self)()
        for key, chrom in self.data.items():
            if key not in other.data:
                for a in chrom:
                    new.add(a)
            else:
                iterateother = Iterate(other.data[key])
                for a in chrom:
                    nottouching = True
                    for b in iterateother.yield_overlapping(a.start, a.stop):
                            nottouching = False
                            break
                    if nottouching:
                        new.add(a)
        return new


    def touched_by(self, other):
        if isinstance(other, Entry): other = Gr(other)
        new = type(self)()
        for key, chrom in self.data.items():
            if key in other.data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other.data[key])
                for a in iterateself.yield_overlapping(other.data[key][0].start, other.data[key].maxstop):
                    for b in iterateother.yield_overlapping(a.start, a.stop):
                        new.add(a)
                        break
        return new

    
    def subranges_covered_by(self, other):
        if isinstance(other, Entry): other = Gr(other)
        new = type(self)()
        for key, chrom in self.data.items():
            if key in other.data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other.data[key])
                for a in iterateself.yield_overlapping(other.data[key][0].start, other.data[key].maxstop):
                    laststop = a.start - 1
                    for b in iterateother.yield_overlapping(a.start, a.stop):
                        if b.start > laststop + 1:
                            break
                        if b.stop > laststop:
                            laststop = b.stop
                    if laststop >= a.stop:
                        new.add(a)
        return new


    def combined_with(self, other):
        if isinstance(other, Entry): other = Gr(other)
        new = type(self)()
        for gr in (self, other):
            for entry in gr:
                new.add(entry)
        new.sort()
        return new


    def subset(self, function):
        return type(self)(entry for entry in self if function(entry))


    @property
    def names(self):
        return set([entry.name for entry in self])


    def __len__(self):
        return self.components


    @property
    def components(self):
        return sum([len(chrom) for chrom in self.data.values()])


    @property
    def weighted_components(self):
        return sum([getattr(entry, "weight", 1) for entry in self])


    @property
    def base_count(self):
        return sum([entry.stop - entry.start + 1 for entry in self])


    @property
    def bases(self):
        return self.base_count


    @property
    def locations_as_string(self):
        return ", ".join([entry.location for entry in self])


    @property
    def names_as_string(self):   ###########################???
        namedict = defaultdict(list)
        for entry in self:
            try:
                namedict[entry.name].append(entry.exon)
            except AttributeError:
                pass
        namelist = []
        for name, numbers in sorted(namedict.items()):
            numbers.sort()
            exons = []
            index = 0
            while index<=len(numbers)-1:
                start = numbers[index]
                stop = start
                for index2 in range(index+1, len(numbers)+1):
                    if index2 == len(numbers) or numbers[index2] != stop+1:
                        break
                    stop += 1
                exons.append("e{0}{1}".format(start, "" if (start==stop) else "-{0}".format(stop)))
                index = index2
            namelist.append("{0} {1}".format(name, ",".join(exons)).strip())
        return ", ".join(namelist)


    @property
    def is_empty(self):
        return self.components == 0


    def save(self, f, filetype="bed"): #Save type(self)() object in bedfile format, START POSITION IS ZERO BASED
        try:
            writer = csv.writer(f, delimiter="\t")
        except TypeError:
            with open(f, "wb") as realfile:
                self.save(realfile, filetype=filetype)
            return
        
        if filetype == "bed":
            for entry in self:
                writer.writerow([entry.chrom, entry.start-1, entry.stop, entry.name, entry.strand])

        elif filetype == "vcf":
            writer.writerow(["##fileformat=VCFv4.2"])
            writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
            try:
                for entry in self:
                    writer.writerow(entry.as_vcf)
            except AttributeError:
                raise CoverMiException("Cannot write a non-variant to a vcf file")
        else:
            raise CoverMiException("Unknown filetype: {}".format(filetype))

            
            
    @property
    def first(self):
        for entry in self:
            return entry


####################################################################################################
# Readers                                                                                          #
####################################################################################################


# Bedfile format 
#    0.  Chromosome
#    1.  Start (zero based)
#    2.  End   (one based)
#    3.  Name
#    4.  Score
#    5.  Strand



def bed(path):
    CHROM = 0
    START = 1
    STOP = 2
    NAME = 3
    STRAND = 5
    with FileContext(path, CSV_READ_MODE) as f:
        for row in csv.reader(f, delimiter="\t"):
            yield Entry(row[CHROM], int(row[START])+1, int(row[STOP]), row[NAME] if len(row)>NAME else DOT, row[STRAND] if len(row)>STRAND else DOT)



# Ilumina Manifest file format:
#    [Header]
#
#    [Probes]
#    0.  Name-ID
#    1.  Region-ID
#    2.  Target-ID 
#    3.  Species
#    4.  Build-ID
#    5.  Chromosome
#    6.  Start position (including primer, zero based)
#    7.  End position (including primer, zero based)
#    8.  Strand
#    9.  Upstream primer sequence (len() will give primer length)
#    10  Upstream hits
#    11. Downstream  primer sequence (len() will give primer length)
#
#    [Targets]
#    0.  Target-ID
#    1.  Target-ID
#    2.  Target number (1 = on target, >1 = off target)
#    3.  Chromosome
#    4.  Start position (including primer, zero based)
#    5.  End position (including primer, zero based)
#    6.  Probe strand
#    7.  Sequence (in direction of probe strand)





# offtarget probes may be wrong
def illuminamanifest(path, ontarget=True, offtarget=False, targets=True, probes=False):
    HEADER = 1
    PROBES = 2
    TARGETS = 3
    section = HEADER
    heading_row = False
    probelen = {}
    rename_offtarget = {}
    with FileContext(path, CSV_READ_MODE) as f:
        for row in csv.reader(f, delimiter="\t"):

            if heading_row:
                col = {heading: index for index, heading in enumerate(row)}
                if section == PROBES:
                    targetid = col["Target ID"]
                    ulso = col["ULSO Sequence"]
                    dlso = col["DLSO Sequence"]
                elif section == TARGETS:
                    chromosome = col["Chromosome"]
                    targeta = col["TargetA"]
                    startpos = col["Start Position"]
                    endpos = col["End Position"]
                    probestrand = col["Probe Strand"]
                    targetnumber = col["Target Number"]
                heading_row = False

            elif row[0].startswith("["):
                section = {"Header": HEADER, "Probes": PROBES, "Targets": TARGETS}[row[0][1:row[0].index("]")]]
                heading_row = section!=HEADER

            elif section == PROBES:
                probelen[row[targetid]] = (len(row[ulso]), len(row[dlso]))

            elif section == TARGETS:
                chrom = row[chromosome]
                name = row[targeta]
                start = int(row[startpos])
                stop = int(row[endpos])
                strand = row[probestrand]
                target_number = row[targetnumber]

                if target_number == "1": 
                    standardised_name = "{0}:{1}-{2}".format(Chrom(chrom), start, stop)
                    rename_offtarget[name] = standardised_name
                else: # Off-target
                    standardised_name = rename_offtarget[name]

                if (target_number == "1" and ontarget) or (target_number != "1" and offtarget):
                    if targets and not probes:
                        yield Entry(chrom, start+probelen[name][strand==PLUS], stop-probelen[name][strand==MINUS], standardised_name, strand)
                    elif targets and probes:
                        yield Entry(chrom, start, stop, standardised_name, strand)
                    elif probes and not targets:
                        yield Entry(chrom, start, start+probelen[name][strand==PLUS]-1, standardised_name, strand)
                        yield Entry(chrom, stop-probelen[name][strand==MINUS]+1, stop, standardised_name, strand)



def variants(path, what, genes=(), diseases=()):
    MUTATION = 0
    GENE = 1
    DISEASE = 2
    what = {"mutation": MUTATION, "gene": GENE, "disease": DISEASE}[what]

    if genes:
        _, genes, _ = load_targets(genes)

    translate = Translate()
    with FileContext(diseases, CSV_READ_MODE) as f:
        for row in csv.reader(f, delimiter="="):
            if len(row) == 2:
                val = row[1].strip()
                if val:
                    translate[row[0].strip()] = val

    def cosmicparser():
        if "Histology subtype" in reader.fieldnames:
            histsubheading = "Histology subtype"
        elif "Histology subtype 1" in reader.fieldnames:
            histsubheading = "Histology subtype 1"
        else:
            raise CoverMiException("Unable to identify HISTOLOGY SUBTYPE heading in COSMIC")

        for row in reader:
            gene = row["Gene name"]
            if gene in genes or not genes:
                primarysite = row["Primary site"]
                histology = row[histsubheading]
                if histology == "NS":
                    histology = row["Primary histology"]
                disease = "{} {}".format(primarysite, histology)
                if row["Mutation genome position"]:
                    chrom, startstop = row["Mutation genome position"].split(":")
                    start, stop = startstop.split(MINUS)
                    try:
                             # chrom  pos         stop       hgvsc                 strand                  disease  gene  hgvsp               transcript
                        yield (chrom, int(start), int(stop), row["Mutation CDS"],  row["Mutation strand"], disease, gene, row["Mutation AA"], row["Accession Number"])
                    except ValueError:
                        pass


    def hgmdparser():
        for row in reader:
            gene = row["Gene Symbol"]
            if gene in genes or not genes:
                hgvsc = row["hgvs"]
                if hgvsc != "null":
                    try:
                             # chrom              pos                           stop                        hgvsc  strand         disease         gene  hgvsp transcript
                        yield (row["chromosome"], int(row["coordinate start"]), int(row["coordinate end"]), hgvsc, row["strand"], row["Disease"], gene, "",   None)
                    except ValueError:
                        pass


    dedup = {}
    with FileContext(path, CSV_READ_MODE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "HGMD ID" in reader.fieldnames:
            parser = hgmdparser
            weight = False
        else:
            parser = cosmicparser
            weight = True

        
        for chrom, pos, stop, hgvsc, strand, disease, gene, hgvsp, transcript_id in parser():
            disease = translate[disease]
            if disease and hgvsc and "?" not in hgvsc:
                refalt = hgvsc.split(DOT)[-1].lstrip("0123456789_-+*()")
                if ">" in refalt:
                    ref, alt = refalt.split(">")
                    if alt.isdigit() or ref == alt: continue
                else:
                    action = refalt[:3]
                    if action == "del":
                        refalt = refalt[3:].split("ins")
                        ref = refalt[0]
                        alt = refalt[1] if len(refalt)==2 else MINUS
                    elif action in ("ins", "dup"):
                        ref = MINUS
                        alt = refalt[3:]
                        pos += 1
                    elif action == "inv":
                        ref = refalt[3:]
                        alt = invert(ref)
                    else:
                        eprint("Unknown HGVS notation: {}".format(hgvsc))
                        continue

                if ref == "" or ref.isdigit(): ref = "N" * (stop - pos + 1)
                if alt.isdigit(): alt = "N"*int(alt)
                ref = ref.strip()
                alt = alt.strip()

                if strand == MINUS:
                    ref, alt = map(invert, (ref, alt))

                if (ref+alt).strip("AGCTURYNWSMKBHDV-") != "":
                    eprint("Unknown HGVS notation: {}".format(hgvsc))
                    continue

                name = gene if what==GENE else (disease if what==DISEASE else "{} {}".format(gene, hgvsp if hgvsp not in ("", "p.?") else hgvsc))
                variant = HgmdVariant(chrom, pos, ref, alt, name, strand)
                variant.gene = gene
                variant.disease = disease
                variant.hgvsc = hgvsc
                if transcript_id:
                    variant.transcript = transcript_id
                identifier = variant.identifier if what!=DISEASE else (variant.identifier, variant.name) # possible to have identical variants with different diseases
                try:
                    variant = dedup[identifier]
                    if what == MUTATION:
                        variant.diseases = "; ".join(set(variant.diseases.split("; ") + [disease]))
                    if weight:
                        variant.weight += 1
                    continue
                except KeyError:
                    pass
                dedup[identifier] = variant
                if what == MUTATION:
                    variant.diseases = disease
                if weight:
                    variant.weight = 1
                yield variant



def firstwhere(f, lower, upper, func):
    #print ("lower={} upper={}".format(lower, upper))

    f.seek(lower)
    row = f.readline()
    if func(row):
        return lower
    firstrowends = f.tell()
    lastrowstarts = ceiling = upper
    while firstrowends != lastrowstarts:
        mid = (lower+upper) // 2
        f.seek(mid)
        f.readline()
        testrowstarts = f.tell()
        row = f.readline()
        #print (row)
        testrowends = f.tell()
        if func(row):
            #print ("DOWN")
            upper = mid
            lastrowstarts = testrowstarts
        elif testrowends >= ceiling:
            return ceiling
        else:
            #print ("UP")
            lower = mid + 1
            firstrowends = testrowends
        #print ("lower={}, upper={}, mid={}, firstrowends={}, lastrowstarts={}".format(lower, upper, mid, firstrowends, lastrowstarts))
    #print ("break")
    f.seek(lastrowstarts)
    #print( f.readline().strip() )
    return lastrowstarts



def vcf(path, locations=(None,), zygosity_data={}):
    hetero_ll = zygosity_data.get("hetero_ll", HETERO_LL)
    hetero_ul = zygosity_data.get("hetero_ul", HETERO_UL)
    homo_ll = zygosity_data.get("homo_ll", HOMO_LL)

    def fao_fdp_depths():
        data = dict(keyval.split("=") for keyval in row[7].split(";"))
        yield int(data["FDP"])
        for d in map(int, data["FAO"].split(",")):
            yield d

    def ao_dp_depths():
        data = dict(keyval.split("=") for keyval in row[7].split(";"))
        yield int(data["DP"])
        for d in map(int, data["AO"].split(",")):
            yield d

    def ad_depths():
        ad = list(map(int, row[9].split(":")[row[8].split(":").index("AD")].split(",")))
        yield sum(ad)
        for d in ad[1:]:
            yield d

    def dp4_dp_depths():
        dp4 = list(map(int, row[7].split("DP4=")[1].split(";")[0].split(",")))
        yield sum(dp4)
        yield sum(dp4[2:4])
        
    def strelka_depths():
        lookup = {key: val for key, val in zip(row[8].split(":"), row[9].split(":"))}
        try:
            yield sum(map(int, (chain(lookup[key].split(",")) for key in ("CU", "GU", "TU", "AU"))))
            key = row[4]+"U"
        except KeyError:
            yield int(lookup["DP"]) + int(lookup["DP2"]) 
            key = "TIR"
        yield sum(map(int, lookup[key].split(",")))

    with FileContext(path, "rt") as f:
        row = "#"
        while row.startswith("#"):
            lower = f.tell()
            row = f.readline()
        if row == "":
            return
        offset = lower

        row = row.strip().split("\t")
        infokeys = [keyval.split("=")[0] for keyval in row[7].split(";")] if len(row)>=8 else ()
        formatkeys = row[8].split(":") if len(row)>=10 else ()
        if "FAO" in infokeys and "FDP" in infokeys:
            depthsfunc = fao_fdp_depths
        elif "AO" in infokeys and "DP" in infokeys:
            depthsfunc = ao_dp_depths
        elif "AD" in formatkeys:
            depthsfunc = ad_depths
        elif "DP4" in infokeys and "DP" in infokeys:
            depthsfunc = dp4_dp_depths
        elif "TU" in formatkeys or "TIR" in formatkeys:
            depthsfunc = strelka_depths
        else:
            depthsfunc = None
        
        hasdepth = bool(depthsfunc)
        hasfilters = bool(row[6] != DOT)
        hasqual = bool(row[5] != DOT)
        variantclass = SequencedVariant if (hasdepth or hasfilters or hasqual) else Variant
        
        f.seek(0, 2)
        vcfend = f.tell()

        index = {}
        while lower < vcfend:
            f.seek(lower)
            row = f.readline()
            chrom = row[0:row.find("\t")]
            upper = firstwhere(f, lower, vcfend, lambda row: row[0:row.find("\t")]!=chrom)
            index[Chrom(chrom)] = (lower, upper)
            lower = upper
        
        for location in locations:
            if location :
                if location.chrom not in index:
                    continue
                offset = firstwhere(f, index[location.chrom][0], index[location.chrom][1], lambda row: int(row.split("\t")[1]) >= location.start)
            f.seek(offset)

            while True:
                row1 = f.readline()
                    
                row = row1.strip().split("\t")
                if row == [""]:
                    break

                chrom = Chrom(row[0])
                pos = int(row[1])
                if location and (chrom != location.chrom or pos > location.stop):
                    break

                if hasfilters:
                    filters = ";".join(code for code in row[6].split(";") if code != "LowVariantFreq")
                    if filters == "":
                        filters = "PASS"
                if hasqual:
                    qual = float(row[5])
                if hasdepth:
                    getdepth = depthsfunc()
                    depth = next(getdepth)
                ref = row[3]
                alts = row[4].split(",")
                name = row[2]
                names = repeat(name) if name == DOT else name.split(",")

                for name, alt in zip(names, alts):
                    if alt not in (ref, DOT):
                        if hasdepth:
                            altdepth = next(getdepth)
                            if altdepth == 0:
                                continue

                        variant = variantclass(chrom, pos, ref, alt, name=row[2])
                        if hasqual:        
                            variant.qual = qual
                        if hasfilters:
                            variant.filters = filters
                        if hasdepth:
                            variant.depth = depth
                            variant.alt_depth = altdepth
                            if homo_ll <= variant.vaf:
                                variant.zygosity = "hemizygous" if chrom in (Chrom(23), Chrom(24), Chrom(25)) else "homozygous"
                            elif hetero_ll <= variant.vaf <= hetero_ul:
                                variant.zygosity = "heterozygous"
                            else:
                                variant.zygosity = "unknown"

                        yield variant



#def samsung(path, slots=False):
#    variantclass = Sequencedvariant if slots else Variant
#    with open(path, "rb") as f:
#        reader = csv.DictReader(f, delimiter="\t")
#        for row in reader:
#            chrom = row["Chr"]
#            pos = int(row["Pos"])
#            ref = row["Ref"]
#            alt = row["Alt"]
#            variant = variantclass(chrom, pos, ref, alt, name="{}:{} {}>{}".format(chrom, pos, ref, alt))
#            variant.alt_depth = int(row["TotAltCnt"])
#            variant.depth = int(row["TotReadCnt"])
#            yield variant


def load_principal(path): # default dict, returns 2 if principal transcript, 1 if alternative and 0 otherwise

    if isinstance(path, defaultdict):
        return path

    principal = defaultdict(int)
    if path:
        with FileContext(path, "rt") as f:
            for row in f:
                row = row.split()
                principal[row[2].split(DOT)[0]] = row[4].startswith("PRINCIPAL") + 1
    return principal


Targets = namedtuple("Targets", "transcript_ids gene_symbols transcript_source")
def load_targets(path):

    if isinstance(path, Targets):
        return path

    gene_re = re.compile(GENE_SYMBOL+"$")
    refseq_re = re.compile(REFSEQ_TRANSCRIPT+"$")
    ensembl_re = re.compile(ENSEMBL_TRANSCRIPT+"$")

    source = None
    needed_transcript_ids = {}
    needed_gene_symbols = set()
    if path:
        with FileContext(path, "rt") as f:
            for row in f:
                cols = row.split()
                numcols = len(cols)
                if numcols:
                    if not gene_re.match(cols[0]):
                        eprint("WARNING Invalid gene symbol: {}".format(cols[0]))
                        continue
                    if numcols == 1:
                        needed_gene_symbols.add(cols[0])
                    elif numcols == 2:
                        if refseq_re.match(cols[1]):
                            thissource = "refseq"
                        elif ensembl_re.match(cols[1]):
                            thissource = "ensembl"
                        else:
                            eprint("WARNING Invalid transcript id: {}".format(col[0]))
                            continue
                        if not source:
                            source = thissource
                        elif thissource != source:
                            eprint("WARNING Mixed ensembl and refseq transcripts, excluding {}".format(col[0]))
                            continue
                        needed_transcript_ids[cols[1].split(DOT)[0]] = cols[0]
                        needed_gene_symbols.add(cols[0])
                    else:
                        eprint("WARNING Malformed targets file: {}".format(row))
    return Targets(needed_transcript_ids, needed_gene_symbols, source)


class Everything(set):
    def __contains__(self, key):
        return True

    def __nonzero__(self): # python2
        return True
    
    def __bool__(self): # python3
        return True
    

class TranscriptList(list):
    def add(self, item):
        self += [item]


class GeneDict(dict):
    def add(self, item):
        gene = item[0]
        try:
            self[gene] += [item]
        except KeyError:
            self[gene] = [item]


# fix splice site buffer
def reference(path, what, principal=None, targets=None, splice_site_buffer=0):
    TRANSCRIPTS = 0
    EXONS = 1
    CODINGREGIONS = 2
    CODINGEXONS = 3
    ENSEMBL = 0
    REFSEQ = 1
    what = {"transcripts": TRANSCRIPTS, "exons": EXONS, "codingregions": CODINGREGIONS, "codingexons": CODINGEXONS}[what]

    source = None
    found_transcript_ids = set()
    found_gene_symbols = set()
    transcript_locations = Counter()
    transcripts = TranscriptList()
    genes = GeneDict()
    sortorder = {}
    entries = []
    transcript = None

    needed_transcript_ids, needed_gene_symbols, targetsource = load_targets(targets)
    if not needed_transcript_ids and not needed_gene_symbols:
        needed_gene_symbols = Everything()

    principal = load_principal(principal if needed_gene_symbols else None)

    with FileContext(path, "rt") as f:
        for row in f:
            row = row.strip("\n ;").split("\t")

            if source is None:
                source = ENSEMBL if row[0].startswith("#!genome-build") else REFSEQ
                if targetsource and source != {"refseq": REFSEQ, "ensembl": ENSEMBL}[targetsource]:
                    raise CoverMiException("ERROR Mixed ensembl and refseq transcripts")

            if source == REFSEQ:
                try:
                    chrom = Chrom(row[2])
                except KeyError:
                    continue

                gene = row[0]
                transcript = row[1]
                if transcript in needed_transcript_ids:
                    if needed_transcript_ids[transcript] not in (gene, None):
                        eprint("WARNING Excluding {} as it pairs with {} not {}".format(transcript, gene, needed_transcript_ids[transcript]))
                        del(needed_transcript_ids[transcript])
                        continue
                    if gene in needed_gene_symbols:
                        needed_gene_symbols.remove(gene)
                        if gene in genes:
                            del(genes[gene])
                    found_transcript_ids.add(transcript)
                    found = transcripts
                elif gene in needed_gene_symbols:
                    found_gene_symbols.add(gene)
                    found = genes
                else:
                    continue

                transcript_locations[transcript] += 1
                strand = row[3]
                transcript_len = 0
                coding_len = 0
                entries = []

                if what == TRANSCRIPTS:
                    entries += [Transcript(chrom, int(row[4])+1-splice_site_buffer, int(row[5])+splice_site_buffer, gene, strand, gene, transcript)]

                elif what == CODINGREGIONS:
                    start = int(row[6])
                    stop = int(row[7])
                    if stop > start: # Otherwise no coding region (start equal to stop due to bed 0, 1 based format)
                        entries += [Transcript(chrom, start+1-splice_site_buffer, stop+splice_site_buffer, gene, strand, gene, transcript)]

                if what == EXONS or found is genes:
                    exon_numbers = range(1,int(row[8])+1) if (strand==PLUS) else range(int(row[8]),0,-1)            
                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
                        start = int(start)+1
                        stop = int(stop)
                        transcript_len += stop - start + 1
                        if what == EXONS:
                            entries += [Exon(chrom, start-splice_site_buffer, stop+splice_site_buffer, gene, strand, gene, transcript, exon)]
                
                if what == CODINGEXONS or found is genes:
                    exon_numbers = range(1,int(row[8])+1) if (strand==PLUS) else range(int(row[8]),0,-1)
                    codingstart = int(row[6])+1
                    codingstop = int(row[7])
                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
                        start = int(start)+1
                        stop = int(stop)
                        if stop >= codingstart and start <= codingstop:
                            coding_len += stop - start + 1
                            if what == CODINGEXONS:
                                entries += [Exon(chrom, 
                                                 (start if start>codingstart else codingstart)-splice_site_buffer,
                                                 (stop if stop<codingstop else codingstop)+splice_site_buffer, 
                                                 gene, strand, gene, transcript, exon)]

                found.add((gene, transcript, entries))
                if found is genes:
                    sortorder[transcript] = (transcript.startswith("N"), principal[transcript], coding_len, transcript_len, -int(transcript[3:]))

            else:
                if row[0][0] != "#":
                    chrom, source, feature, start, stop, score, strand, phase, keyvals = row
                    keyvals = {key: val for key, val in [keyval.split() for keyval in keyvals.split("; ")]}
                    start = int(start)
                    stop = int(stop)

                    if feature == "transcript":
                        if entries:
                            found.add((gene, transcript, entries))
                            entries = []
                            if found is genes:
                                sortorder[transcript] = (source == "ensembl_havana", principal[transcript], coding_len, transcript_len, -int(transcript[4:]))

                        try:
                            chrom = Chrom(chrom)
                        except KeyError:
                            transcript = None
                            continue
                        gene = keyvals["gene_name"].strip("\"")
                        transcript = keyvals["transcript_id"].strip("\"")
                        transcript_len = 0
                        coding_len = 0
                        cds_startstop = None
                        if transcript in needed_transcript_ids:
                            if needed_transcript_ids[transcript] not in (gene, None):
                                eprint("WARNING Excluding {} as it pairs with {} not {}".format(transcript, gene, needed_transcript_ids[transcript]))
                                del(needed_transcript_ids[transcript])
                                transcript = None
                                continue
                            if gene in needed_gene_symbols:
                                del(needed_gene_symbols[gene])
                                if gene in genes:
                                    del(genes[gene])
                            found_transcript_ids.add(transcript)
                            found = transcripts
                        elif gene in needed_gene_symbols:
                            found_gene_symbols.add(gene)
                            found = genes
                        else:
                            transcript = None
                            continue

                        if what == TRANSCRIPTS:
                            entries += [Transcript(chrom, start-splice_site_buffer, stop+splice_site_buffer, gene, strand, gene, transcript)]

                    elif transcript:
                        if feature == "exon":
                            if what == EXONS:
                                entries += [Exon(chrom, start-splice_site_buffer, stop+splice_site_buffer, gene, 
                                                 strand, gene, transcript, int(keyvals["exon_number"].strip("\"")))]
                            transcript_len += stop - start + 1

                        elif feature == "CDS":
                            if what == CODINGEXONS:
                                entries += [Exon(chrom, start-splice_site_buffer, stop+splice_site_buffer, gene, 
                                                 strand, gene, transcript, int(keyvals["exon_number"].strip("\"")))]
                            elif what == CODINGREGIONS and cds_startstop is None:
                                cds_startstop = start if strand == PLUS else stop
                            coding_len += stop - start + 1

                        elif feature == "stop_codon":
                            if what == CODINGREGIONS:
                                if strand == PLUS:
                                    cds_start = cds_startstop
                                    cds_stop = stop
                                else:
                                    cds_start = start
                                    cds_stop = cds_startstop
                                entries += [Transcript(chrom, int(cds_start)-splice_site_buffer, int(cds_stop)+splice_site_buffer, gene, strand, gene, transcript)]
                            elif what == CODINGEXONS:
                                if strand == PLUS and entries[-1].stop+1 == start:
                                    entries[-1].stop = stop
                                elif strand == MINUS and entries[-1].start-1 == stop:
                                    entries[-1].start = start
                                else:
                                    entries += [Exom(chrom, start-splice_site_buffer, stop+splice_site_buffer, gene, strand, gene, transcript, entries[-1].exon+1)]

    # Warning duplicate code to ensure cached entries are processed
    if source == ENSEMBL and entries:
        found.add((gene, transcript, entries))
        entries = []
        if found is genes:
            sortorder[transcript] = (source == "ensembl_havana", principal[transcript], coding_len, transcript_len, -int(transcript[4:]))

    notfound = list(needed_gene_symbols - set(found_gene_symbols)) + list(set(needed_transcript_ids) - set(found_transcript_ids))
    if notfound:
        eprint("WARNING {} not found in reference file".format(", ".join(notfound)))

    gene_transcripts = {}
    multi_transcript_genes = set()
    for gene, transcript, entries in transcripts:
        try:
            if gene_transcripts[gene] != transcript:
                multi_transcript_genes.add(gene)
        except KeyError:
            gene_transcripts[gene] = transcript

    for gene, transcript, entries in transcripts:
        for entry in entries:
            if gene in multi_transcript_genes:
                entry.name = "{} {}".format(entry.name, transcript)
            if transcript_locations[transcript] > 1:
                entry.name = "{} {}".format(entry.name, entry.location)
            yield entry

    for candidates in genes.values():
        candidates = sorted(candidates, key=lambda x: sortorder[x[1]], reverse=True)
        winning_transcript = candidates[0][1]
        for gene, transcript, entries in candidates:
            if transcript != winning_transcript:
                break
            for entry in entries:
                if transcript_locations[transcript] > 1:
                    entry.name = "{} {}".format(entry.name, entry.location)
                yield entry













#def reference(path, what, canonical=(), genes=(), splice_site_buffer=SPLICE_SITE_BUFFER):
#    TRANSCRIPTS = 0
#    EXONS = 1
#    CODINGREGIONS = 2
#    CODINGEXONS = 3
#    ENSEMBL = 0
#    REFSEQ = 1
#    what = {"transcripts": TRANSCRIPTS, "exons": EXONS, "codingregions": CODINGREGIONS, "codingexons": CODINGEXONS}[what]

#    with open(path, "rt") as f:
#        testline = f.read(1000).split("\n")[0]
#    source = ENSEMBL if testline.startswith("#!genome-build") else REFSEQ

#    neededtranscripts = {}
#    neededgenes = set()
#    excludedgenes = set()
#    foundgenes = {} # store first transcript found
#    foundtranscripts = Counter() # count number of copies of each transcript found
#    multitranscriptgenes = set()
#    geteverything = False
#    entries = []
#    
#    transcript_re = re.compile((REFSEQ_TRANSCRIPT if source==REFSEQ else ENSEMBL_TRANSCRIPT) + "$")
#    if genes:
#        with FileContext(genes, "rt") as f:
#            for row in f:
#                cols = row.split()
#                numcols = len(cols)
#                if numcols > 0:
#                    firstcol = cols[0]
#                    if numcols == 1:
#                        if transcript_re.match(firstcol):
#                            neededtranscripts[firstcol.split(DOT)[0]] = None
#                        else:
#                            neededgenes.add(firstcol)
#                    elif numcols == 2:
#                        neededtranscripts[cols[1].split(DOT)[0]] = firstcol
#                    else:
#                        eprint("WARNING Malformed targets file: {}".format(row))

#        if neededgenes:
#            with FileContext(canonical, "rt") as f:
#                for row in f:
#                    cols = row.split()
#                    if len(cols) == 2:
#                        gene, transcript = col
#                        if gene in neededgenes:
#                            neededgenes.remove(gene)
#                            neededtranscripts[transcript] = gene
#                    else:
#                        eprint("WARNING Malformed canonical file: {}".format(row))

#    else:
#        with FileContext(canonical, "rt") as f:
#            for row in f:
#                cols = row.split()
#                if len(cols) == 2:
#                    gene, transcript = col
#                    neededtranscripts[transcript] = gene
#                    excludedgenes.add(gene)
#                else:
#                    eprint("WARNING Malformed canonical file: {}".format(row))

#        if not neededtranscripts:
#             geteverything = True

#    with FileContext(path, "rb") as f:
#        if source == REFSEQ:
#            for row in csv.reader(f, delimiter="\t"):
#                try:
#                    chrom = Chrom(row[2])
#                except KeyError:
#                    continue

#                transcript = row[1]
#                gene = row[0]
#                if transcript in neededtranscripts:
#                    if neededtranscripts[transcript] not in (gene, None):
#                        eprint("WARNING Excluding {} as it pairs with {} not {}".format(transcript, gene, neededtranscripts[transcript]))
#                        del(neededtranscripts[transcript])
#                        continue
#                elif gene in neededgenes or (excludedgenes and gene not in excludedgenes) or geteverything:
#                    pass
#                else:
#                    continue

#                try:
#                    if foundgenes[gene] != transcript:
#                        multitranscriptgenes.add(gene)
#                except KeyError:
#                    foundgenes[gene] = transcript
#                foundtranscripts[transcript] += 1

#                strand = row[3]
#                identifier = (gene, transcript)

#                if what == TRANSCRIPTS:
#                    entries += [Transcript(chrom, int(row[4])+1-splice_site_buffer, int(row[5])+splice_site_buffer, identifier, strand, gene, transcript)]

#                elif what == CODINGREGIONS:
#                    start = int(row[6])
#                    stop = int(row[7])
#                    if stop > start: # Otherwise no coding region (start equal to stop due to bed 0, 1 based format)
#                        entries += [Transcript(chrom, start+1-splice_site_buffer, stop+splice_site_buffer, identifier, strand, gene, transcript)]

#                elif what == EXONS:
#                    exon_numbers = range(1,int(row[8])+1) if (strand==PLUS) else range(int(row[8]),0,-1)            
#                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
#                        entries += [Exon(chrom, int(start)+1-splice_site_buffer, int(stop)+splice_site_buffer, identifier, strand, gene, transcript, exon)]

#                elif what == CODINGEXONS:
#                    exon_numbers = range(1,int(row[8])+1) if (strand==PLUS) else range(int(row[8]),0,-1)
#                    codingstart = int(row[6])+1
#                    codingstop = int(row[7])
#                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
#                        start = int(start)+1
#                        stop = int(stop)
#                        if stop >= codingstart and start <= codingstop:
#                            entries += [Exon(chrom, 
#                                             (start if start>codingstart else codingstart)-splice_site_buffer,
#                                             (stop if stop<codingstop else codingstop)+splice_site_buffer, 
#                                             identifier, strand, gene, transcript, exon)]

#        else:
#            for row in f:
#                if row[0] != "#":
#                    chrom, source, feature, start, stop, score, strand, phase, keyvals = row.strip("\n ;").split("\t")
#                    if feature != "gene":
#                        keyvals = {key: val for key, val in [keyval.split() for keyval in keyvals.split("; ")]}

#                        if feature == "transcript":
#                            try:
#                                chrom = Chrom(chrom)
#                            except KeyError:
#                                transcript = None
#                                continue
#                            gene = keyvals["gene_name"].strip("\"")
#                            transcript = keyvals["transcript_id"].strip("\"")
#                            identifier = (gene, transcript)
#                            cds_startstop = None

#                            if transcript in neededtranscripts:
#                                if neededtranscripts[transcript] not in (gene, None):
#                                    eprint("WARNING Excluding {} as it pairs with {} not {}".format(transcript, gene, neededtranscripts[transcript]))
#                                    del(neededtranscripts[transcript])
#                                    transcript = None
#                            elif not(gene in neededgenes or (excludedgenes and gene not in excludedgenes) or geteverything ):
#                                transcript = None

#                            if transcript:
#                                try:
#                                    if foundgenes[gene] != transcript:
#                                        multitranscriptgenes.add(gene)
#                                except KeyError:
#                                    foundgenes[gene] = transcript
#                                foundtranscripts[transcript] += 1

#                                if what == TRANSCRIPTS:
#                                    entries += [Transcript(chrom, int(start)-splice_site_buffer, int(stop)+splice_site_buffer, identifier, strand, gene, transcript)]                        

#                        elif transcript:
#                            if what == EXONS and feature == "exon":
#                                entries += [Exon(chrom, int(start)-splice_site_buffer, int(stop)+splice_site_buffer, identifier, strand, gene, transcript, int(keyvals["exon_number"].strip("\"")))]

#                            if what == CODINGEXONS and feature == "CDS":
#                                entries += [Exon(chrom, int(start)-splice_site_buffer, int(stop)+splice_site_buffer, identifier, strand, gene, transcript, int(keyvals["exon_number"].strip("\"")))]

#                            elif what == CODINGREGIONS and feature == "CDS" and cds_startstop is None:
#                                cds_startstop = start if strand == PLUS else stop

#                            elif what == CODINGREGIONS and feature == "stop_codon":
#                                if strand == PLUS:
#                                    cds_start = cds_startstop
#                                    cds_stop = stop
#                                else:
#                                    cds_start = start
#                                    cds_stop = cds_startstop
#                                entries += [Transcript(chrom, int(cds_start)-splice_site_buffer, int(cds_stop)+splice_site_buffer, identifier, strand, gene, transcript)]



#    notfound = list(neededgenes - set(foundgenes)) + list(set(neededtranscripts) - set(foundtranscripts))
#    if notfound:
#        eprint("WARNING {} not found in refflat file".format(", ".join(notfound)))

#    for entry in entries:
#        gene, transcript = entry.name
#        entry.name = gene if gene not in multitranscriptgenes else "{} {}".format(gene, transcript)
#        if foundtranscripts[transcript] > 1:
#            entry.name = "{} {}".format(entry.name, entry.location)
#        yield entry


#    


class Fasta(object):
    def __init__(self, path):
        self.fastas = {}
        regexp = re.compile(r"\.chromosome\.([0-9XYMT]+)\.")
        for fn in os.listdir(path):
            if os.path.isfile(os.path.join(path, fn)):
                match = regexp.search(fn)
                if match:
                    try:
                        self.fastas[Chrom(match.group(1))] = os.path.join(path, fn)
                    except KeyError:
                        pass
        if not self.fastas:
            raise CoverMiException("No fasta file found within {}".format(path))

    def __call__(self, entry):
        fasta = self.fastas[entry.chrom]
        with (GzipFile if fasta.endswith(".gz") else open)(fasta, "rb") as f:
            offset = 0
            for line in f:
                if line.startswith(">"):
                    offset += len(line)
                else:
                    linebytes = len(line)
                    linebases = len(line.strip())
                    break
            startbyte = offset + (((entry.start-1) // linebases) * linebytes) + ((entry.start-1) % linebases)
            stopbyte = offset + (((entry.stop-1) // linebases) * linebytes) + ((entry.stop-1) % linebases)
            f.seek(startbyte)
            sequence = f.read(stopbyte - startbyte + 1)
        return "".join(sequence.split()).upper()



def gel_csv(path):
    igv_regexp = re.compile(">(.+):(.+);(.+)>(.+)</a>")
    with FileContext(path, CSV_READ_MODE) as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            try:
                igv_coords = row["GRCh38 coordinates;ref/alt allele"]
                match = igv_regexp.search(igv_coords)
                chrom, pos, ref, alt = match.groups()
                alt_depth, depth = [int(reads) for reads in row["Alt allele/total read depth"].split("/")]
                name = row["cDNA and protein change"].split(";")[-1]
            except (ValueError, KeyError, AttributeError):
                eprint("Bad imput in GEL CSV file: {}".format(os.path.basenane(path)))
                continue
            variant = SequencedVariant(chrom, pos, ref, alt, name, DOT)
            variant.depth = depth
            variant.alt_depth = alt_depth
            yield variant






def yield_vcf(path, locations=(None,)):

    with FileContext(path, "rt") as f:
        row = "#"
        while row.startswith("#"):
            lower = f.tell()
            row = f.readline()
        if row == "":
            return

        offset = lower
        f.seek(0, 2)
        vcfend = f.tell()

        index = {}
        while lower < vcfend:
            f.seek(lower)
            row = f.readline()
            chrom = row[0:row.find("\t")]
            upper = firstwhere(f, lower, vcfend, lambda row: row[0:row.find("\t")]!=chrom)
            index[Chrom(chrom)] = (lower, upper)
            lower = upper
        
        for location in locations:
            if location :
                if location.chrom not in index:
                    continue
                offset = firstwhere(f, index[location.chrom][0], index[location.chrom][1], lambda row: int(row.split("\t")[1]) >= location.start)
            f.seek(offset)

            while True:
                row1 = f.readline()
                    
                row = row1.strip().split("\t")
                if row == [""]:
                    break

                chrom = Chrom(row[0])
                pos = int(row[1])
                if location and (chrom != location.chrom or pos > location.stop):
                    break

                filters = ";".join(code for code in row[6].split(";") if code != "LowVariantFreq")
                if filters == "":
                    filters = "PASS"
                        
                if filters == "PASS":
                    yield row1






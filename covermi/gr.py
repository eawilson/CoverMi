import csv, sys, pdb
from itertools import islice, repeat, chain
from collections import defaultdict, namedtuple
import collections.abc
from copy import copy



__all__ = ["Entry", "Gr", "bed", "appris", "reference"]


class Entry(object):
    __slots___ = ("chrom", "start", "stop", "name", "strand")

    def __init__(self, chrom, start, stop, name=".", strand="."):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.name = name
        self.strand = strand

    def __repr__(self):
        return "{}({}, {}, {}, {}, {})".format(type(self).__name__, repr(self.chrom), repr(self.start), repr(self.stop), repr(self.name), repr(self.strand))

    def __str__(self):
        return "{}:{}-{}".format(self.chrom, self.start, self.stop)
    
    def __eq__(self, other):
        return self._tuple == other._tuple
    
    def __lt__(self, other):
        return self._tuple < other._tuple
    
    def __hash__(self):
        return hash(self._tuple)
    
    def __iter__(self):
        yield self
    
    @property
    def _tuple(self):
        return (self.chrom, self.start, self.stop, self.name, self.strand)



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



class GrList(list):
    def __init__(self):
        super().__init__()
        self.maxlength = 0
        self.maxstop = 0

        

class Gr(collections.abc.Mapping):
    def __init__(self, iterable=()):
        self._data = defaultdict(GrList)
        with self as gr:
            for entry in iterable:
                gr.add(entry)


    def __repr__(self):
        return "{}({})".format(type(self).__name__, ", ".join(repr(entry) for entry in self))


    def __len__(self):
        return sum([len(chrom) for chrom in self._data.values()])


    def __iter__(self):
        for key, chrom in sorted(self._data.items()):
            yield from chrom


    def __getitem__(self, key):
        return self._data[key]


    def __enter__(self):
        return self


    def __exit__(self, type_, value, traceback):
        self.sort()
    
    
    def keys(self):
        return self._data.keys()
        

    def add(self, entry):
        grlist = self._data[entry.chrom]
        grlist.append(entry)
        length = entry.stop - entry.start + 1
        if length > grlist.maxlength:
            grlist.maxlength = length
        if entry.stop > grlist.maxstop:
            grlist.maxstop = entry.stop
        return self


    def sort(self):
        for chrom in self._data.values():
            chrom.sort()


    @property
    def merged(self): 
        new = type(self)()
        for key, chrom in self._data.items():
            new.add(chrom[0])
            nchrom = new._data[key]
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
        for key, chrom in self._data.items():
            if key in other._data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other._data[key])
                for a in iterateself.yield_overlapping(other._data[key][0].start, other._data[key].maxstop):
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
        for key, chrom in self._data.items():
            if key not in other._data:
                for a in chrom:
                    new.add(a)
            else:
                iterateother = Iterate(other._data[key])
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
        for key, chrom in self._data.items():
            if key in other._data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other._data[key])
                for a in iterateself.yield_overlapping(other._data[key][0].start, other._data[key].maxstop):
                    for b in iterateother.yield_overlapping(a.start, a.stop):
                        new.add(a)
                        break
        return new

    
    def covered_by(self, other):
        if isinstance(other, Entry): other = Gr(other)
        new = type(self)()
        for key, chrom in self._data.items():
            if key in other._data:
                iterateself = Iterate(chrom)
                iterateother = Iterate(other._data[key])
                for a in iterateself.yield_overlapping(other._data[key][0].start, other._data[key].maxstop):
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
        new = type(self)()
        for entry in chain(self, other):
            new.add(entry)
        new.sort()
        return new


    @property
    def names(self):
        return set([entry.name for entry in self])


    @property
    def bases(self):
        bases = 0
        for entry in self:
            bases += entry.stop - entry.start + 1
        return bases



def bed(path):
    with open(path, "rt") as f:
        for row in csv.reader(f, delimiter="\t"):
            row[1] = int(row[1]) + 1
            row[2] = int(row[2])
            row[3] = " ".join(row[3].split())
            yield Entry(*row)



class MissingDict(dict):
    def __init__(self, func):
        self._func = func
    
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, func.__name__, super().__repr__())
    
    def __missing__(self, key):
        return self._func()



def appris(path):
    # returns 2 if principal transcript, 1 if alternative and 0 otherwise
    score = MissingDict(int)
    if path:
        with open(path, "rt") as f:
            for row in f:
                row = row.split()
                score[row[2].split(".")[0]] = row[4].startswith("PRINCIPAL") + 1
    return score



def reference(path, what, names=(), principal=""):
    TRANSCRIPTS = 0
    EXONS = 1
    CODINGREGIONS = 2
    CODINGEXONS = 3
    ENSEMBL = 0
    REFSEQ = 1
    try:
        what = {"transcripts": TRANSCRIPTS, "exons": EXONS, "codingregions": CODINGREGIONS, "codingexons": CODINGEXONS}[what]
    except KeyError:
        raise ValueError("invalid 'what' for reference()")
    
    needed_genes = set()
    needed_transcripts = MissingDict(str)
    for name in names:
        name = name.split()
        if len(name) == 1:
            needed_genes.add(name[0])
        elif len(name) == 2:
            needed_transcripts[name[0]] = name[1]
    
    # returns 2 if principal transcript, 1 if alternative and 0 otherwise
    score = appris(principal)
    
    transcript = None
    buffer = []
    found = set()
    found_genes = defaultdict(lambda:defaultdict(list))
    found_transcripts = []
    sort_order = {}
    source = None
    with open(path, "rt") as f:
        for row in f:
            row = row.strip("\n ;").split("\t")

            if source is None:
                source = ENSEMBL if row[0].startswith("#!genome-build") else REFSEQ

            if source == REFSEQ:
                gene = row[0]
                transcript = row[1]
                
                if transcript == needed_transcripts[gene]:
                    name = f"{gene} {transcript}"
                    entries = found_transcripts
                elif gene in needed_genes:
                    name = gene
                    entries = found_genes[gene][transcript]
                else:
                    continue
                found.add(name)
                
                chrom = row[2]
                strand = row[3]
                transcript_len = 0
                coding_len = 0
                if what == TRANSCRIPTS:
                    entries.append(Entry(chrom, int(row[4])+1, int(row[5]), name, strand))

                elif what == CODINGREGIONS:
                    start = int(row[6]) + 1
                    stop = int(row[7])
                    if stop >= start:
                        entries.append(Entry(chrom, start, stop, name, strand))

                if what == EXONS or entries is not found_transcripts:
                    exon_numbers = range(1,int(row[8])+1) if (strand=="+") else range(int(row[8]),0,-1)            
                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
                        start = int(start) + 1
                        stop = int(stop)
                        transcript_len += stop - start + 1
                        if what == EXONS:
                            #name = f"{name} e{exon}"
                            entries.append(Entry(chrom, start, stop, name, strand))
                
                if what == CODINGEXONS or entries is not found_transcripts:
                    exon_numbers = range(1,int(row[8])+1) if (strand=="+") else range(int(row[8]),0,-1)
                    codingstart = int(row[6]) + 1
                    codingstop = int(row[7])
                    for start, stop, exon in zip(row[9].rstrip(",").split(","), row[10].rstrip(",").split(","), exon_numbers):
                        start = int(start) + 1
                        stop = int(stop)
                        if stop >= codingstart and start <= codingstop:
                            start = max(codingstart, start)
                            stop = min(codingstop, stop)
                            coding_len += stop - start + 1
                            if what == CODINGEXONS:
                                #name = f"{name} e{exon}"
                                entries.append(Entry(chrom, start, stop, name, strand))
                
                if entries is not found_transcripts:
                    sort_order[transcript] = (score[transcript], coding_len, transcript_len, -int(transcript[3:]))

            else:
                if not row[0].startswith("#"):
                    #chrom, source, feature, start, stop, score, strand, phase, keyvals = row
                    feature = row[2]

                    if feature == "transcript":
                        keyvals = {key: val for key, val in [keyval.split() for keyval in row[8].split("; ")]}
                        gene = keyvals["gene_name"].strip('"')
                        transcript = keyvals["transcript_id"].strip('"')
                        
                        if transcript == needed_transcripts[gene]:
                            name = f"{gene} {transcript}"
                            entries = found_transcripts
                        elif gene in needed_genes:
                            name = gene
                            entries = found_genes[gene][transcript]
                        else:
                            transcript = None
                            continue
                        found.add(name)
                        cds_startstop = None

                        if what == TRANSCRIPTS:
                            entries.append(Entry(row[0], int(row[3]), int(row[4]), name, row[6]))
                        
                        if entries is not found_transcripts:
                            sort_order[transcript] = [row[1] == "ensembl_havana", score[transcript], 0, 0, -int(transcript[4:])]

                    elif transcript:
                        keyvals = {key: val for key, val in [keyval.split() for keyval in row[8].split("; ")]}
                        start = int(row[3])
                        stop = int(row[4])
                        if feature == "exon":
                            if what == EXONS:
                                #"{} e{}".format(name, keyvals["exon_number"].strip('"'))
                                entries.append(Entry(row[0], start, stop, name, row[6]))
                            if entries is not found_transcripts:
                                sort_order[transcript][3] += stop - start + 1

                        elif feature == "CDS":
                            strand = row[6]
                            if what == CODINGEXONS:
                                #"{} e{}".format(name, keyvals["exon_number"].strip('"'))
                                entries.append(Entry(row[0], start, stop, name, strand))
                            elif what == CODINGREGIONS and cds_startstop is None:
                                cds_startstop = start if strand == "+" else stop
                            if entries is not found_transcripts:
                                sort_order[transcript][2] += stop - start + 1

                        elif feature == "stop_codon":
                            strand = row[6]
                            if what == CODINGREGIONS:
                                if strand == "+":
                                    cds_start = cds_startstop
                                    cds_stop = stop
                                else:
                                    cds_start = start
                                    cds_stop = cds_startstop
                                entries.append(Entry(row[0], int(cds_start), int(cds_stop), name, strand))
                            elif what == CODINGEXONS:
                                if strand == "+" and entries[-1].stop + 1 == start:
                                    entries[-1].stop = stop
                                elif strand == "-" and entries[-1].start - 1 == stop:
                                    entries[-1].start = start
                                else:
                                    # Never tested by the looks of it! ? exon+1 correct or strand dependent
                                    #"{} e{}".format(name, entries[-1].exon+1)
                                    entries.append(Entry(row[0], start, stop, name, strand))

    notfound = ", ".join(sorted(set(names) - found))
    if notfound:
        print(f"WARNING: {notfound} not found in reference file", file=sys.stderr)

    for entry in found_transcripts:
        yield entry

    for candidates in found_genes.values():
        for entry in sorted(candidates.items(), key=lambda x: sort_order[x[0]])[-1][1]:
            yield entry



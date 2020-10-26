import csv, pdb
import sys
import gzip
from itertools import islice, repeat, chain, cycle
from collections import defaultdict, namedtuple
import collections.abc
from copy import copy



nucleotide = {"A": "R", "G": "R", "C": "Y", "T": "Y"}# puRine: A, G,  pYrimadine: T, C
standard_chrom = {"1": "chr1",
                  "2": "chr2",
                  "3": "chr3",
                  "4": "chr4",
                  "5": "chr5",
                  "6": "chr6",
                  "7": "chr7",
                  "8": "chr8",
                  "9": "chr9",
                  "10": "chr10",
                  "11": "chr11",
                  "12": "chr12",
                  "13": "chr13",
                  "14": "chr14",
                  "15": "chr15",
                  "16": "chr16",
                  "17": "chr17",
                  "18": "chr18",
                  "19": "chr19",
                  "20": "chr20",
                  "21": "chr21",
                  "22": "chr22",
                  "23": "chrX",
                  "chr23": "chrX",
                  "x": "chrX",
                  "chrx": "chrX",
                  "X": "chrX",
                  "24": "chrY",
                  "chr24": "chrY",
                  "y": "chrY",
                  "chry": "chrY",
                  "Y": "chrY",
                  "25": "chrM",
                  "chr25": "chrM",
                  "m": "chrM",
                  "chrm": "chrM",
                  "M": "chrM"}
                  


__all__ = ["Entry", "Gr", "Variant", "bed", "appris", "reference", "vcf", "DepthAltDepths", "MissingDict"]



def gzopen(fn, *args, **kwargs):
    return (gzip.open if fn.endswith(".gz") else open)(fn, *args, **kwargs)



class Entry(object):
    __slots___ = ("chrom", "start", "stop", "name", "strand")

    def __init__(self, chrom, start, stop, name=".", strand="."):
        self.chrom = standard_chrom.get(chrom, chrom)
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



class Variant(Entry):
    __slots___ = ("ref", "alt", "depth", "alt_depth")

    def __str__(self):
        return "{}:{} {}/{}".format(self.chrom, self.start, self.ref, self.alt)

    def __init__(self, chrom, pos, ref, alt, name=".", depth=None, alt_depth=None):
        if ref == alt:
            raise ValueError("Alt allele cannot be equal to ref allele")
        while ref and alt and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1
        while ref and alt and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]

        self.ref = ref or "-"
        self.alt = alt or "-"
        self.depth = depth
        self.alt_depth = alt_depth
        super().__init__(chrom, pos, pos-1 if self.ref=="-" else pos+len(ref)-1, name)

    @property
    def pos(self):
        return self.start

    @property
    def vartype(self):
        if self.ref == "-":
            return "ins"
        elif self.alt == "-":
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
    def vaf(self):
        return float(self.alt_depth)/self.depth

    @property
    def zygosity(self):
        vaf = self.vaf
        if 0.95 <= vaf:
            zygosity = "hemizygous" if self.chrom in ("chrX", "chrY", "chrM") else "homozygous"
        elif 0.45 <= vaf <= 0.55:
            zygosity = "heterozygous"
        else:
            zygosity = "unknown"
        return zygosity



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



class Deduplicate(dict):
    def __missing__(self, key):
        self[key] = key
        return key



def bed(path):
    deduplicate = Deduplicate()
    with open(path, "rt") as f:
        for row in csv.reader(f, delimiter="\t"):
            row[0] = deduplicate[row[0]]
            row[1] = int(row[1]) + 1
            row[2] = int(row[2])
            yield Entry(*row[:min(5, len(row))])



class MissingDict(dict):
    def __init__(self, func):
        self._func = func
    
    def __repr__(self):
        return "{}({}, {})".format(type(self).__name__, self._func.__name__, super().__repr__())
    
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
        elif len(name) > 1:
            needed_transcripts[name[0]] = name[1]
    
    # returns 2 if principal transcript, 1 if alternative and 0 otherwise
    score = appris(principal)
    
    deduplicate = Deduplicate()
    transcript = None
    buffer = []
    found = set()
    found_genes = defaultdict(lambda:defaultdict(list))
    found_transcripts = []
    sort_order = {}
    source = None
    with gzopen(path, "rt") as f:
        for row in f:
            row = row.rstrip("\n ;").split("\t")

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
                
                chrom = deduplicate[row[2]]
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
                                    entries.append(Entry(deduplicate[row[0]], start, stop, name, strand))

    notfound = ", ".join(sorted(set(names) - found))
    if notfound:
        print(f"WARNING: {notfound} not found in reference file", file=sys.stderr)

    for entry in found_transcripts:
        yield entry

    for candidates in found_genes.values():
        for entry in sorted(candidates.items(), key=lambda x: sort_order[x[0]])[-1][1]:
            yield entry



def vcf(path, name=".", format_index=9):
    depth_alt_depths = DepthAltDepths(format_index)
    with open(path, "rt") as f:
        for row in f:
            row = row.rstrip("\n ;").split("\t")
            depth, alt_depths = depth_alt_depths(row)
            for alt, alt_depth in zip(row[4].split(","), alt_depths):
                if alt_depth and alt not in (ref, "."):
                    yield Variant(row[0], int(row[1]), row[3], alt, name=name, depth=depth, alt_depth=alt_depth)



class DepthAltDepths(object):
    def __init__(self, format_values=9):
        self._format_values = format_values
        
    def __call__(self, row):
        try:
            return self._func(row)
        except AttributeError:
            pass
        self._func = self.choose_algorithm(row)
        return self._func(row)


    def choose_algorithm(self, row):
        infodict = dict(kv.split("=") for kv in row[7].split(";") if "=" in kv)
        
        for self.key1, self.key2, algorithm in (("FRO", "FAO", self.ro_ao_depth), ("RO", "AO", self.ro_ao_depth), ("FDP", "FAO", self.dp_ao_depth), ("DP", "AO", self.dp_ao_depth)):
            if self.key1 in infodict and self.key2 in infodict:
                return algorithm
        
        if len(row) > self._format_values: # has a format section
            formatdict = dict(zip(row[8].split(":"), row[self._format_values].split(":")))
            
            # In illumina vcfs AD stands for allelic depths and contains a comma separated list of ref and all alt depths.
            # In some other vcfs AD stands for alt depth and contains the depth of the alt read only with RD containing the ref depth!!!!!!
            if len(row[4].split(",")) + 1 == len(formatdict.get("AD", "").split(",")): # must have depth for ref and each alt.
                return self.ad_depth
            
            if "AD" in formatdict and "RD" in formatdict:
                return self.ad_rd_depth
            
            if "DP4" in formatdict:
                return self.dp4_format_depth
            
        if "DP4" in infodict:
            return self.dp4_info_depth
        
        if len(row) > self._format_values:
            if all(key in formatdict for key in ["GU", "CU", "AU", "TU"]) or all(key in formatdict for key in ["TAR", "TIR"]):
                return self.strelka_depth
            
        return self.no_depth
            
    
    def ro_ao_depth(self, row):
        info = dict(kv.split("=") for kv in row[7].split(";") if "=" in kv)
        ref_depth = int(info[self.key1])
        alt_depths =  [int(d) for d in info[self.key2].split(",")]
        return (ref_depth + sum(alt_depths), alt_depths)


    def dp_ao_depth(self, row):
        info = dict(kv for kv in row[7] if len(kv) == 2)
        tot_depth = int(info[self.key1])
        alt_depths =  [int(d) for d in info[self.altkey].split(",")]
        return (tot_depth, alt_depths)


    def ad_depth(self, row):
        ad = [int(depth) for depth in row[self._format_values].split(":")[row[8].split(":").index("AD")].split(",")]
        return (sum(ad), ad[1:])


    def ad_rd_depth(self, row):
        keys = row[8].split(":")
        vals = row[self._format_values].split(":")
        ad = int(vals[keys.index("AD")])
        return (ad + int(vals[keys.index("RD")]), [ad])


    def dp4_info_depth(self, row):
        dp4 = [int(depth) for depth in row[7].split("DP4=")[1].split(";")[0].split(",")]
        return (sum(dp4), [dp4[2] + dp4[3]])


    def dp4_format_depth(self, row):
        dp4 = [int(depth) for depth in dict(zip(row[8].split(":"), row[self._format_values].split(":")))["DP4"].split(",")]
        return (sum(dp4), [dp4[2] + dp4[3]])
        

    def strelka_depth(self, row):
        # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
        formatdict = dict(zip(row[8].split(":"), row[self._format_values].split(":")))
        alt = row[4]
        if "," in alt:
            raise RuntimeError("Multiple variants per row in strelka vcf.")
        if "CU" in formatdict:
            depth = sum(int(formatdict[key].split(",")[0]) for key in ("AU", "TU", "CU", "GU"))
            alt_depth = int(formatdict[f"{alt}U"].split(",")[0])
        else:
            alt_depth = int(formatdict["TIR"].split(",")[0])
            depth = int(formatdict["TAR"].split(",")[0]) + alt_depth
        return (depth, [alt_depth])
    

    def no_depth(self, row):
        return (1, [1])




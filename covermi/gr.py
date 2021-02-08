import csv
import pdb
import sys
import gzip
import os
from itertools import islice, chain
from collections import defaultdict
import collections.abc
from copy import copy



__all__ = ["Entry", "Gr", "Variant", "bed", "appris", "gff3", "vcf", "DepthAltDepths"]



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
                  "M": "chrM",
                  "NC_000001": "chr1",
                  "NC_000002": "chr2",
                  "NC_000003": "chr3",
                  "NC_000004": "chr4",
                  "NC_000005": "chr5",
                  "NC_000006": "chr6",
                  "NC_000007": "chr7",
                  "NC_000008": "chr8",
                  "NC_000009": "chr9",
                  "NC_000010": "chr10",
                  "NC_000011": "chr11",
                  "NC_000012": "chr12",
                  "NC_000013": "chr13",
                  "NC_000014": "chr14",
                  "NC_000015": "chr15",
                  "NC_000016": "chr16",
                  "NC_000017": "chr17",
                  "NC_000018": "chr18",
                  "NC_000019": "chr19",
                  "NC_000020": "chr20",
                  "NC_000021": "chr21",
                  "NC_000022": "chr22",
                  "NC_000023": "chrX",
                  "NC_000024": "chrY",
                  "NC_012920": "chrM"}



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



def bed(paths):
    
    try:
        paths = (os.fspath(paths),)
    except TypeError:
        pass

    for path in paths:
        with gzopen(path, "rt") as f:
            for row in csv.reader(f, delimiter="\t"):
                row[1] = int(row[1]) + 1
                row[2] = int(row[2])
                if len(row) < 6:
                    yield Entry(*row)
                else:
                    yield Entry(*row[:5])



def appris(paths):
    # returns 2 if principal transcript, 1 if alternative
    
    try:
        paths = (os.fspath(paths),)
    except TypeError:
        pass

    score = {}
    for path in paths:
        with gzopen(path, "rt") as f:
            for row in f:
                row = row.split()
                score[row[2].split(".")[0]] = row[4].startswith("PRINCIPAL") + 1
    return score



def bases(components):
    return sum(c.stop - c.start + 1 for c in components)



def cannonical(elements):
    return (bases(elements.get("CDS", ())), 
            bases(elements.get("exon", ())), 
            bases(elements.get("transcript", ())))



SEQID = 0
SOURCE = 1
TYPE = 2
START = 3
END = 4
SCORE = 5
STRAND = 6
PHASE = 7
ATTRIBUTES = 8

#def gff3(path, what, names=(), principal=""):
    
    #if what not in ("transcripts", "exons", "codingregions", "codingexons"):
        #raise ValueError(f"Invalid value for what: {what}")
    
    #what = what[:-1]
    #codingregion = what == "codingregion"
    #if what == "codingexon":
        #what = "CDS"
    
    #needed = defaultdict(list)
    #for name in names:
        #splitname = name.split()
        #if splitname:
            #needed[splitname[0]].extend(splitname[1:])
    
    #score = appris(principal)    
    #matches = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    
    #with gzopen(path, "rt") as f_in:
        #match = False
        #reader = csv.reader(f_in, delimiter="\t")
        #for row in reader:
            #if row[0].startswith("#"):
                #continue
            
            #if match:
                #if ";Parent=" in row[ATTRIBUTES]:
                    #feature = row[TYPE]
                    #if feature.endswith("transcript") or feature.endswith("RNA"):
                        #try:
                            #stop = row[ATTRIBUTES].index(";")
                        #except ValueError:
                            #stop = len(row[ATTRIBUTES])
                        #transcript = row[ATTRIBUTES][3:stop].split(".")[0]
                        #feature = "transcript"
                    
                    #entry = Entry(chrom, int(row[START]), int(row[END]), gene, strand)
                    #matches[gene][transcript][feature].append(entry)
                    #continue
                #match = False
        
            #if row[TYPE] == "gene":
                #try:
                    #start = row[ATTRIBUTES].index(gene_name)
                #except NameError:
                    #if "gene_name=" in row[ATTRIBUTES]:
                        #gene_name = "gene_name="
                    #elif "Name=" in row[ATTRIBUTES]:
                        #gene_name = "Name="
                    #else:
                        #raise ValueError("Gene name not found in attributes")
                    #start = row[ATTRIBUTES].index(gene_name)
                        
                #start += len(gene_name)
                #try:
                    #stop = row[ATTRIBUTES].index(";", start)
                #except ValueError:
                    #stop = len(row[ATTRIBUTES])
                #gene = row[ATTRIBUTES][start:stop]
                #if not names or gene in needed:
                    #chrom = row[SEQID].split(".")[0]
                    #chrom = standard_chrom.get(chrom)
                    #if chrom: # Remove patches
                        #strand = row[STRAND]
                        #match = True
    
    
    #for gene in sorted(set(needed) - set(matches)):
        #print(f"WARNING: {gene} not found in reference file", file=sys.stderr)
    
    #for gene, transcripts in sorted(matches.items()):
        #selected = needed[gene]
        #if not selected:
            #if len(transcripts) == 1:
                #selected = transcripts.keys()
            #else:
                #scored = ([], [])
                #for transcript in transcripts:
                    #try:
                        #scored[score[transcript] - 1].append(transcript)
                    #except KeyError:
                        #pass
                
                #candidates = scored[1] or scored[0] or transcripts
                #if len(candidates) == 1:
                    #selected = candidates
                #else:
                    #selected = sorted(candidates, key=lambda t:cannonical(transcripts[t]))[-1:]
        
        #for transcript in selected:
            #try:
                #features = transcripts[transcript]
            #except KeyError:
                #print(f"WARNING: {transcript} not found in reference file", file=sys.stderr)
                #continue
            
            #if codingregion:
                #first = features["CDS"][0]
                #last = features["CDS"][-1]
                #start = min(first.start, last.start)
                #stop = max(first.stop, last.stop)
                #features["codingregion"] = (Entry(first.chrom, start, stop, first.name, first.strand),)
            
            #for entry in features[what]:
                #if len(selected) > 1:
                    #entry.name = f"{entry.name} {transcript}"
                #yield entry


ENSEMBL_REFSEQ_PREFIXES = set(["ENS", "NM_", "NR_", "XM_", "XR_"])


def gff3(paths, what, names=(), principal=()):
    
    if what not in ("transcripts", "exons", "codingregions", "codingexons"):
        raise ValueError(f"Invalid value for what: {what}")
    
    what = what[:-1] # remove trailing s
    codingregion = what == "codingregion"
    if what in ("codingregion", "codingexon"):
        what = "CDS"
    
    try:
        paths = (os.fspath(paths),)
    except TypeError:
        pass
    
    needed = defaultdict(list)
    for name in names:
        splitname = name.split()
        if splitname:
            needed[splitname[0]].extend(splitname[1:])
    
    score = appris(principal)    
    matches = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    
    gene_name = ""
    for path in paths:
        del gene_name
        with gzopen(path, "rt") as f_in:
            transcript = ""
            reader = csv.reader(f_in, delimiter="\t")
            for row in reader:
                if row[0].startswith("#"):
                    continue
                
                feature = row[TYPE]
                
                if feature.endswith("transcript") or feature.endswith("RNA"):
                    transcript = ""
                    attributes = row[ATTRIBUTES]
                    try:
                        start = attributes.index(gene_name)
                    except NameError:
                        if "gene_name=" in attributes:
                            gene_name = "gene_name="
                        elif "Name=" in attributes:
                            gene_name = "gene="
                        else:
                            raise ValueError("Gene name not found in transcript attributes")
                        start = attributes.index(gene_name)
                    except ValueError:
                        continue
                            
                    start += len(gene_name)
                    try:
                        stop = attributes.index(";", start)
                    except ValueError:
                        stop = len(attributes)
                    gene = attributes[start:stop]
                    
                    if not names or gene in needed:
                        chrom = row[SEQID].split(".")[0]
                        chrom = standard_chrom.get(chrom, chrom)
                        strand = row[STRAND]
                        feature = "transcript"
                        start = attributes.index("transcript_id=") + 14
                        try:
                            stop = attributes.index(";", start)
                        except ValueError:
                            stop = len(attributes)
                        transcript = attributes[start:stop]
                        if transcript[:3] in ENSEMBL_REFSEQ_PREFIXES:
                            transcript = transcript.split(".")[0]
                                                
                if transcript and feature in ("transcript", "exon", "CDS"):
                    entry = Entry(chrom, int(row[START]), int(row[END]), gene, strand)
                    matches[gene][transcript][feature].append(entry)
    
    
    for gene in sorted(set(needed) - set(matches)):
        print(f"WARNING: {gene} not found in reference file", file=sys.stderr)
    
    for gene, transcripts in sorted(matches.items()):
        selected = needed[gene]
        if not selected:
            if len(transcripts) == 1:
                selected = transcripts.keys()
            else:
                scored = ([], [])
                for transcript in transcripts:
                    try:
                        scored[score[transcript] - 1].append(transcript)
                    except KeyError:
                        pass
                
                candidates = scored[1] or scored[0] or transcripts
                if len(candidates) == 1:
                    selected = candidates
                else:
                    selected = sorted(candidates, key=lambda t:cannonical(transcripts[t]))[-1:]
        
        for transcript in selected:
            if transcript not in transcripts:
                print(f"WARNING: {transcript} not found in reference file", file=sys.stderr)
                continue
            
            features = transcripts[transcript]
            if codingregion and "CDS" in features:
                first = features["CDS"][0]
                last = features["CDS"][-1]
                yield Entry(first.chrom,
                            min(first.start, last.start),
                            max(first.stop, last.stop),
                            first.name,
                            first.strand)

            else:
                for entry in features.get(what, ()):
                    if len(selected) > 1:
                        entry.name = f"{entry.name} {transcript}"
                    yield entry



def vcf(paths, name=".", format_index=9):
    
    try:
        paths = (os.fspath(paths),)
    except TypeError:
        pass

    for path in paths:
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




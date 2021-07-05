import csv
import pdb
import sys
import gzip
import os
from itertools import islice, chain
from collections import defaultdict
import collections.abc
from copy import copy


__all__ = ["Entry", "Gr", "Variant", "bed", "appris", "gff3", "vcf", "depth_alt_depth_function", "chromosome_order"]



nucleotide = {"A": "R", "G": "R", "C": "Y", "T": "Y"}# puRine: A, G,  pYrimadine: T, C

ENSEMBL_REFSEQ_PREFIXES = set(["ENS", "NM_", "NR_", "XM_", "XR_"])

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



def chromosome_order(key):
    if key.startswith("chr"):
        key = key[3:]
    return ("", int(key)) if key.isnumeric() else (key, 0)



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
        for key in sorted(self._data.keys(), key=chromosome_order):
            yield from self._data[key]


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
        with gzopen(path, "rt") as f_in:
            for row in f_in:
                row = row.strip()
                if row:
                    row = row.split("\t")
                    row[1] = int(row[1]) + 1
                    row[2] = int(row[2])
                    if len(row) < 5:
                        yield Entry(*row)
                    else:
                        yield Entry(*row[:4])



def appris(paths):
    # returns 2 if principal transcript, 1 if alternative
    score = {}
    
    if paths:
        try:
            paths = (os.fspath(paths),)
        except TypeError:
            pass

        for path in paths:
            with gzopen(path, "rt") as f:
                for row in f:
                    row = row.split("\t")
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
        
        # These 3 variabls are used for weeding out refseq duplicate genes
        gene_attr = ""
        refseq = True
        valid_gene = True
        
        with gzopen(path, "rt") as f_in:
            transcript = ""
            reader = csv.reader(f_in, delimiter="\t")
            for row in reader:
                if row[0].startswith("#"):
                    continue
                
                feature = row[TYPE]
                
                if feature == "gene":
                    gene_attr = row[ATTRIBUTES]
                
                if (feature.endswith("transcript") or feature.endswith("RNA")):
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
                        
                        # Remove refseq duplicate genes
                        if refseq:
                            name = ""
                            ident = ""
                            for attr in gene_attr.split(";"):
                                if attr.startswith("ID="):
                                    ident = attr[8:]
                                elif attr.startswith("Name="):
                                    name = attr[5:]
                            if ident and name:
                                valid_gene = (ident == name)
                            else:
                                refseq = False
                        
                        if valid_gene:
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
                        else:
                            print(f"Excluding duplicate gene {ident}", file=sys.stderr)
                                                
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


CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8

def vcf(paths, name="."):
    
    try:
        paths = (os.fspath(paths),)
    except TypeError:
        pass

    depth_alt_depth = None
    for path in paths:
        del depth_alt_depth
        with open(path, "rt") as f:
            for row in f:
                if not row.startswith("#"):
                    row = row.rstrip("\n ;").split("\t")
                    try:
                        dp, ad = depth_alt_depth(row)
                    except NameError:
                        depth_alt_depth = depth_alt_depth_function(row)
                        dp, ad = depth_alt_depth(row)
                    if ad != 0 and row[ALT] not in (row[REF], "."):
                        yield Variant(row[CHROM], int(row[POS]), row[REF], row[ALT], name, dp, ad)
                            


def vd_dp_dad(row):
    info = infodict(row)
    return (int(info["DP"]), int(info["VD"]))



def ao_ro_dad(row):
    info = infodict(row)
    ad = int(info["AO"])
    return (int(info["RO"]) + ad, ad)



def fao_fro_dad(row):
    info = infodict(row)
    ad = int(info["FAO"])
    return (int(info["FRO"]) + ad, ad)



def ao_dp_dad(row):
    info = infodict(row)
    return (int(info["DP"]), int(info["AO"]))



def fao_fdp_dad(row):
    info = infodict(row)
    return (int(info["FDP"]), int(info["FAO"]))



def ad_dad(row):
    ref, alt = formatdict(row)["AD"].split(",")
    ad = int(alt)
    return (int(ref) + ad, ad)



def ad_rd_dad(row):
    fmt = formatdict(row)
    ad = int(fmt["AD"])
    return (int(fmt["RD"]) + ad, ad)



def dp4_format_dad(row):
    fmt = formatdict(row)
    dp4 = [int(num) for num in fmt["DP4"].split(",")]
    ad = dp4[2] + dp4[3]
    return (dp4[0] + dp4[1] + ad, ad)



def dp4_info_dad(row):
    info = infodict(row)
    dp4 = [int(num) for num in info["DP4"].split(",")]
    ad = dp4[2] + dp4[3]
    return (dp4[0] + dp4[1] + ad, ad)

    

def strelka_dad(row):
    # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
    fmt = formatdict(row)
    alt = row[ALT]
    if "," in alt:
        raise RuntimeError("Multiple variants per row")
    ad = int(fmt[f"{alt}U"].split(",")[0])
    dp = sum(int(fmt[key].split(",")[0]) for key in ("AU", "TU", "CU", "GU"))
    return (dp, ad)



def tar_tir_dad(row):
    # https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic
    fmt = formatdict(row)
    alt = row[4]
    if "," in alt:
        raise RuntimeError("Multiple variants per row")
    ad = int(fmt["TIR"].split(",")[0])
    dp = int(fmt["TAR"].split(",")[0]) + ad
    return (dp, ad)



def no_dad(row):
    return (None, None)



def infodict(row):
    infodict = {}
    for token in row[INFO].split(";"):
        try:
            k, v = token.split("=")
        except ValueError:
            k = token
            v = None
        infodict[k] = v
    return infodict



def formatdict(row):
    return dict(zip(row[FORMAT].split(":"), row[FORMAT+1].split(":")))



def depth_alt_depth_function(row):
    ### Will fall down if a vcf contains multiple variants on the same line
    ###
    if "," in row[ALT]:
        raise RuntimeError("Multiple variants per row")
    
    info = infodict(row)
    fmt = formatdict(row) if len(row) > FORMAT + 1 else {}

    # Vardict vcf
    if "VD" in fmt and "DP" in info:
        return vd_dp_dad

    if "FRO" in info and "FAO" in info:
        return fao_fro_dad
    
    if "RO" in info and "AO" in info:
        return ao_ro_dad
    
    if "FDP" in info and "FAO" in info:
        return fao_fdp_dad
    
    if "DP" in info and "AO" in info:
        return ao_dp_dad
    
    # In illumina vcfs AD stands for allelic depths and contains a comma separated list of ref and all alt depths.
    # In some other vcfs AD stands for alt depth and contains the depth of the alt read only with RD containing the ref depth!!!!!!
    if "," in fmt.get("AD", ()): # Will fail with multiple variants on same line
        return ad_dad

    if "AD" in fmt and "RD" in fmt:
        return ad_rd_dad
    
    if "DP4" in fmt:
        return dp4_format_dad
    
    if all(key in fmt for key in ("GU", "CU", "AU", "TU")):
        return strelka_dad
        
    if "TAR" in fmt and "TIR" in fmt:
        return tar_tir_dad
    
    if "DP4" in info:
        return dp4_info_dad
        
    return no_dad
            
    

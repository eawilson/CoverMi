#
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
#s
#    [Targets]
#    0.  Target-ID
#    1.  Target-ID
#    2.  Target number (1 = on target, >1 = off target)
#    3.  Chromosome
#    4.  Start position (including primer, zero based)
#    5.  End position (including primer, zero based)
#    6.  Probe strand
#    7.  Sequence (in direction of probe strand)
# 
#
# Bedfile format 
#    0.  Chromosome
#    1.  Start (zero based)
#    2.  End   (one based)
#    3.  Name
#    4.  Score
#    5.  Strand
# 
#
#  Variant file format
#
# 
#
# 
#
# 
#

import copy, pdb

class CoverMiException(Exception):
    pass

class Entry(object):

    def __repr__(self):
        return "[\"{0}\", {1}, {2}, \"{3}\", \"{4}\"]".format(self.chrom, self.start, self.stop, self.name, self.strand)

    def __lt__(self, other):
        return self.start<other.start

    def __init__(self, chrom, start, stop, name, strand):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        self.strand = strand


class Gr(dict):

    KARYOTYPE = {"chr1":1, "chr2":2, "chr3":3, "chr4":4, "chr5":5, "chr6":6, "chr7":7, "chr8":8, "chr9":9, "chr10":10, "chr11":11, "chr12":12, 
                 "chr13":13, "chr14":14, "chr15":15, "chr16":16, "chr17":17, "chr18":18, "chr19":19, "chr20":20, "chr21":21, "chr22":22, "chrX":23, "chrY":24, "chrM":25}
    @classmethod
    def chromosome_number(genomicrange, chromosome): return genomicrange.KARYOTYPE[chromosome]
    MAX_CHR_LENGTH = 1000000000
    SPLICE_SITE_BUFFER = 5


    def __init__(self, entry=None):
        if entry is not None:
            self.construct(copy.copy(entry))


    @classmethod
    def load_bedfile(genomicrange, path):
        with file(path, "rU") as f:
            gr = genomicrange()
            for line in f:
                splitline = line.strip().split("\t")
                if splitline != [""]:
                    name = splitline[3] if len(splitline)>3 else "."
                    strand = splitline[5] if len(splitline)>5 else "."
                    gr.construct(Entry(splitline[0], 
                                        int(splitline[1])+1,
                                        splitline[2],
                                        name,
                                        strand))
        gr.sort()
        return gr


    @classmethod
    def load_manifest(genomicrange, path, ontarget=True, offtarget=False):
        with file(path, "rU") as f:
            gr = genomicrange()
            section = "Header"
            skip_column_names = False
            probes = {}
            rename_offtarget = {}

            for line in f:
                splitline = line.rstrip("\n").split("\t")

                if skip_column_names:
                    skip_column_names = False

                elif line[0] == "[":
                    section = line[1:].split("]")[0]
                    if section != "Header":
                            skip_column_names = True

                elif section == "Probes":
                    probes[splitline[2]] = (len(splitline[9]), len(splitline[11]))

                elif section == "Targets":
                    if splitline[2] == "1":
                        amplicon_name = "{0}:{1}-{2}".format(splitline[3], splitline[4], splitline[5])
                        rename_offtarget[splitline[0]] = amplicon_name
                    else:
                        amplicon_name = rename_offtarget[splitline[0]]
                    if ((splitline[2] == "1") and ontarget) or ((splitline[2] != "1") and offtarget):
                        gr.construct(Entry(splitline[3],
                                            int(splitline[4])+probes[splitline[0]][splitline[6]=="+"], 
                                            int(splitline[5])-probes[splitline[0]][splitline[6]=="-"], 
                                            amplicon_name, 
                                            splitline[6]))
        gr.sort()
        return gr


    @classmethod
    def load_refflat(genomicrange, path, genes_of_interest=None, canonical_transcripts=None):
        canonicaldict = {}
        if canonical_transcripts is not None:
            for line in canonical_transcripts:
                gene, transcript = line.strip().split("\t")
                canonicaldict[gene] = transcript.split(" ")[-1]

        needed = set([])
        doublecheck = {}
        include_everything = True
        if genes_of_interest is not None:
            for line in genes_of_interest:
                line = line.strip()
                if line != "":
                    include_everything = False
                    if line in canonicaldict:
                        needed.add(canonicaldict[line])
                        doublecheck[canonicaldict[line]] = line
                    else:
                        splitline = line.split()
                        if len(splitline) == 1:
                            needed.add(line)
                        elif len(splitline) == 2:
                            needed.add(splitline[1])
                            doublecheck[splitline[1]] = splitline[0]
                        else:
                            raise CoverMiException("Malformed line in gene list: {0}".format(line))
            
        with file(path, "rU") as f:
            multiple_copies = {}
            found = set([])
            exons = genomicrange()
            transcripts = genomicrange()
            for line in f:
                splitline = line.rstrip("\n").split("\t")
                if splitline[2] in genomicrange.KARYOTYPE:
                    if include_everything:
                        if splitline[0] in canonicaldict and splitline[1] != canonicaldict[splitline[0]]:
                            continue
                    elif splitline[0] in needed:
                        found.add(splitline[0])
                    elif splitline[1] in needed and splitline[0] == doublecheck[splitline[1]]:
                        found.add(splitline[1])
                    else:
                        continue

                    transcript = "{0} {1}".format(splitline[0], splitline[1])
                    name = transcript
                    if name not in multiple_copies:
                        multiple_copies[name] = 1
                    else:
                        if multiple_copies[name] == 1:
                            for gr in (exons, transcripts):
                                for entry in gr.all_entries:
                                    if entry.name == name:
                                        entry.name = "{0} copy 1".format(name)
                        multiple_copies[name] += 1
                        name = "{0} copy {1}".format(name, multiple_copies[name])
              
                    entry = Entry(splitline[2], int(splitline[4])+1-genomicrange.SPLICE_SITE_BUFFER, int(splitline[5])+genomicrange.SPLICE_SITE_BUFFER,  name, splitline[3])
                    entry.transcript = transcript
                    transcripts.construct(entry)

                    exon_numbers = range(1,int(splitline[8])+1) if (splitline[3] == "+") else range(int(splitline[8]),0,-1)            
                    for start, stop, exon in zip(splitline[9].rstrip(",").split(","), splitline[10].rstrip(",").split(","), exon_numbers):
                        entry = Entry(splitline[2], int(start)+1-genomicrange.SPLICE_SITE_BUFFER, int(stop)+genomicrange.SPLICE_SITE_BUFFER, name, splitline[3])
                        entry.exon = exon
                        exons.construct(entry)

        missing = needed - found
        if not include_everything and len(missing) > 0:
            missing = ["{0} {1}".format(doublecheck[trans], trans) if trans in doublecheck else trans for trans in missing]
            print "WARNING - {0} not found in reference file".format(", ".join(missing))        

        for gr in (exons, transcripts):
            gr.sort()

            #simplify gene names
            names = {}
            duplicates = set([])
            for entry in gr.all_entries:
                gene, transcript = entry.name.split(" ")[0:2]
                if gene in names and names[gene] != transcript:
                    duplicates.add(gene)
                else:
                    names[gene] = transcript
            for entry in gr.all_entries:
                splitgene = entry.name.split(" ")
                if splitgene[0] not in duplicates:
                    if len(splitgene) <= 2:
                        entry.name = splitgene[0]
                    else:
                        entry.name = splitgene[0]+" "+" ".join(splitgene[2:])

        return (exons, transcripts)


    @classmethod
    def load_variants(genomicrange, path, catagory, genes_of_interest=None, disease_names=None): 
        if genes_of_interest is not None:
            genes_of_interest = set([name.split()[0] for name in genes_of_interest])

        disease_dict = {}
        if disease_names is not None:
            for line in disease_names:
                if line[0] != "#":
                    line2 = line.rstrip().split("=")
                    if len(line2) == 2:
                        disease_dict[line2[0].strip("\" ")] = line2[1].strip("\" ")

        with file(path, "rU") as f:
            mutations = {}
            for line in f:
                splitline = line.rstrip("\n").split("\t")
                if splitline == [""] or splitline[4] not in genomicrange.KARYOTYPE:
                    continue
                if genes_of_interest is not None and splitline[3] not in genes_of_interest:
                    continue    
                mutation = "{0} {1} {2}:{3}-{4}".format(splitline[3], splitline[8], splitline[4], splitline[5], splitline[6])# gene mutation location
                disease = splitline[1].strip("\" ")
                if disease in disease_dict:
                    if disease_dict[disease] == "":
                        continue
                    else:
                        disease = disease_dict[disease]
                if catagory == "disease":
                    name = disease 
                    mutation = (mutation, name)
                elif catagory == "gene":
                    name = splitline[3]
                elif catagory == "mutation":
                    name = mutation
 
                if mutation not in mutations:
                    mutations[mutation] = Entry(splitline[4], splitline[5], splitline[6], name, splitline[7])
                    if catagory == "mutation": 
                        mutations[mutation].diseases=set([])
                    mutations[mutation].weight = 1
                else:
                    mutations[mutation].weight += 1
                if catagory == "mutation":
                    mutations[mutation].diseases.add(disease)

        gr = genomicrange()
        for entry in mutations.values():
            if catagory == "mutation":
                entry.diseases = "; ".join(sorted(entry.diseases))
            gr.construct(entry)
        gr.sort()
        return gr


    @classmethod
    def load_vcf(genomicrange, path, range_of_interest=None): #range_of_interest = genomic range covering variants to be loaded
        gr = genomicrange()
        if range_of_interest is not None:
            range_of_interest = range_of_interest.merged
        with file(path, "rU") as f:
            for line in f:
                if not line.startswith("#"):
                    splitline = line.split("\t")
                    chrom, untrimmedpos, snpid, untrimmedref, alts, qual, filters, junk = splitline[0:8]
                    if len(splitline) == 10:
                        headings, values = splitline[8:10]
                    assert len(splitline) == 8 or len(splitline) == 10
                    
                    if not chrom.startswith("chr"):
                        chrom = "chr{0}".format(chrom)
                    if chrom == "chr23":
                        chrom = "chrX"
                    elif chrom == "chr24":
                        chrom = "chrY"
                    elif chrom == "chr25":
                        chrom = "chrM"

                    quality = float(qual) if (qual!=".") else "."
                    passfilter = "." if filters == "." else all([code in ("PASS", "LowVariantFreq") for code in filters.split(";")])
#                    if filters != ".":
#                        passfilter = True
#                        for code in filters.split(";"):
#                            if code in ("LowGQX", "LowGQ"):
#                                passfilter = False
#                            elif code not in ("PASS", "LowVariantFreq", "R8"):
#                                passfilter = False
#                                print code
#                                #raise CoverMiException("Unknown filter code {0} in VCF {1}".format(code, path))
                    try:
                        alt_depths = [ float(depth) for depth in values.split(":")[headings.split(":").index("AD")].split(",") ]
                        total_depth = sum(alt_depths)
                    except (ValueError, UnboundLocalError):
                        pass

                    untrimmedpos = int(untrimmedpos)
                    alts = alts.split(",")
                    for alt_number in range(0, len(alts)):
                        alt = alts[alt_number]
                        ref = untrimmedref
                        start = untrimmedpos
                        while len(ref) > 1 and len(alt) > 1:
                            if ref[0:2] == alt[0:2]:
                                ref = ref[1:]
                                alt = alt[1:]
                                start += 1
                            elif ref[-1] == alt[-1]:
                                ref = ref[:-1]
                                alt = alt[:-1]
                            else:
                                break
                        refalt_length_difference = len(ref) - len(alt)
                        if refalt_length_difference == 0:#snp
                            stop = start
                        elif refalt_length_difference < 0:#alt longer than ref therefore insertion
                            stop = start + 1
                        else:#deletion 
                            stop = start + 1 + refalt_length_difference

                        if range_of_interest is not None:
                            if chrom not in range_of_interest:
                                continue
                            covered = False
                            for roi_entry in range_of_interest[chrom]:
                                if roi_entry.start <= start and roi_entry.stop >= stop:
                                    covered = True
                                    break
                                if roi_entry.start > stop:
                                    break
                            if not covered:
                                continue

                        entry = Entry(chrom, start, stop, "{0}:{1} {2}>{3}".format(chrom, start, ref, alt), "+")#?????strand
                        entry.indel = refalt_length_difference != 0
                        if quality != ".":
                            entry.quality = quality
                        if passfilter != ".":
                            entry.passfilter = passfilter
                        try:
                            if total_depth == 0:
                                continue
                            else:
                                entry.vaf = alt_depths[alt_number+1]/total_depth
                        except UnboundLocalError:
                            pass
                        gr.construct(entry)
        gr.sort()
        return gr


    def construct(self, entry):################################# ?add
        if entry.chrom not in self:
            self[entry.chrom] = []
        self[entry.chrom].append(entry)
        return self
            

    def sort(self):
        for chrom in self:
            self[chrom].sort()


    @property
    def all_entries(self):
        for chrom in sorted(self, key=type(self).chromosome_number):
            for entry in self[chrom]:
                yield entry

    @property
    def merged(self): # The name of a range will be the name of the first sub-range that was merged into it
        merged = type(self)()
        for chrom in self:
            merged[chrom] = [copy.copy(self[chrom][0])]
            for entry in self[chrom]:
                if entry.start-1 <= merged[chrom][-1].stop:
                    merged[chrom][-1].stop = max(entry.stop, merged[chrom][-1].stop)
                else:
                    merged.construct(copy.copy(entry))
        merged.sort()
        return merged


    def extended_to_include_touching_amplicons(self, gr):
         extended = type(self)()
         amplicons = type(self)()########################################## = ?Gr
         for chrom in self:
             for a in self[chrom]:
                 new_entry = copy.copy(a)
                 if chrom in gr:
                     for b in gr[chrom]:
                         if b.start > a.stop:
                             break
                         if (a.start <= b.start <= a.stop) or (a.start <= b.stop <= a.stop) or b.stop > a.stop:
                             new_entry.start = min(new_entry.start, b.start)
                             new_entry.stop = max(new_entry.stop, b.stop)
                             amplicons.construct(copy.copy(b))
                 extended.construct(new_entry)
         return (extended, amplicons)


    def overlapped_by(self, gr):# If range gr contains overlapping ranges then we may get multiple copies of the overlapping ranges
         trimmed = type(self)()
         for chrom in self:
            if chrom in gr:
                for a in self[chrom]:
                    for b in gr[chrom]:
                        if b.start > a.stop:
                            break
                        if (a.start <= b.start <= a.stop) or (a.start <= b.stop <= a.stop) or b.stop > a.stop:
                            new_entry = copy.copy(a)
                            new_entry.start = max(a.start, b.start)
                            new_entry.stop = min(a.stop, b.stop)
                            trimmed.construct(new_entry)
         trimmed.sort() # Only needed if gr contains overlapping ranges which should not happen with sane use
         return trimmed


    def not_touched_by(self, gr):
        untouched = type(self)() 
        for chrom in self:
           for a in self[chrom]:
               touching = False
               if chrom in gr:
                   for b in gr[chrom]:
                       if b.start > a.stop:
                           break
                       if (a.start <= b.start <= a.stop) or (a.start <= b.stop <= a.stop) or b.stop > a.stop:
                           touching = True
                           break
               if not touching:
                   untouched.construct(copy.copy(a))
        untouched.sort()
        return untouched


    def touched_by(self, gr):
         touched = type(self)() 
         for chrom in self:
            for a in self[chrom]:
                if chrom in gr:
                    for b in gr[chrom]:
                        if b.start > a.stop:
                            break
                        if (a.start <= b.start <= a.stop) or (a.start <= b.stop <= a.stop) or b.stop > a.stop:
                            touched.construct(copy.copy(a))
                            break
         touched.sort()
         return touched

    
    def subranges_covered_by(self, gr):
         covered = type(self)()
         for chrom in self:
            if chrom in gr:
                for a in self[chrom]:
                    for b in gr[chrom]:
                        if b.start > a.start:
                            break
                        if b.stop >= a.stop:
                            covered.construct(copy.copy(a))
                            break
         covered.sort()
         return covered 


    def combined_with(self, other):
        combined = type(self)()
        for gr in (self, other):
            for chrom in gr:
                for entry in gr[chrom]:
                    combined.construct(copy.copy(entry))
        combined.sort()
        return combined


    def subset(self, names, exclude=False, genenames=False):
        if not hasattr(names, "__iter__"):
            names = [names]
        if genenames:
            names = [name.split()[0] for name in names]
        gr = type(self)()
        for entry in self.all_entries:
            if ((entry.name.split()[0] if genenames else entry.name) in names) ^ exclude:
                gr.construct(copy.copy(entry))
        return gr


    @property
    def names(self):
        names = set([])
        for entry in self.all_entries:
            names.add(entry.name)
        return sorted(names)


    @property
    def number_of_components(self):
        output = 0
        for chrom in self:
            output += len(self[chrom])
        return output


    @property
    def number_of_weighted_components(self):
        output = 0
        for entry in self.all_entries:
            output += entry.weight
        return output


    @property
    def base_count(self):
        output = 0
        for entry in self.all_entries:
            output += entry.stop - entry.start + 1
        return output


    @property
    def locations_as_string(self):
        output = []
        for entry in self.all_entries:
            output.append("{0}:{1}-{2}".format(entry.chrom, entry.start, entry.stop))
        return ", ".join(output)


    @property
    def names_as_string(self):   
        namedict = {}
        for entry in self.all_entries:
            if entry.name not in namedict:
                namedict[entry.name] = []
            try:
                namedict[entry.name].append(entry.exon)
            except AttributeError:
                pass
        namelist = []
        for name, numbers in namedict.items():
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
            namelist.append("{0} {1}".format(name, ",".join(exons)))
            namelist.sort()
        return ", ".join(namelist).strip()####################?strip needed


    @property
    def is_empty(self):
        return len(self) == 0


    def save(self, f): #Save type(self)() object in bedfile format, START POSITION IS ZERO BASED
        try:
            for entry in self.all_entries:
                f.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(entry.chrom, entry.start-1, entry.stop, entry.name, entry.strand))
        except AttributeError:
            with file(f, "wt") as realfile:
                self.save(realfile)
            


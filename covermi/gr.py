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

import re, pdb

KARYOTYPE = {"chr1":1, "chr2":2, "chr3":3, "chr4":4, "chr5":5, "chr6":6, "chr7":7, "chr8":8, "chr9":9, "chr10":10, "chr11":11, "chr12":12, 
             "chr13":13, "chr14":14, "chr15":15, "chr16":16, "chr17":17, "chr18":18, "chr19":19, "chr20":20, "chr21":21, "chr22":22, "chrX":23, "chrY":24, "chrM":25}
def chromosome_number(chromosome): return KARYOTYPE[chromosome]
MAX_CHR_LENGTH = 1000000000
SPLICE_SITE_BUFFER = 5

CHROM = 0
START = 1
STOP = 2
NAME = 3
STRAND = 4
NAME2 = 5
WEIGHT = 6

def load_manifest(path, excluded=[]):
    with file(path, "rU") as f:
        gr1 = GenomicRange()
        #gr2 = GenomicRange()
        section = "Header"
        skip_column_names = False
        probes = {}
        #rename_offtarget = {}

        for line in f:
            splitline = line.rstrip("\n").split("\t")

            if skip_column_names:
                skip_column_names = False

            elif line[0] == "[":
                section = line[1:len(line)-2]
                if section != "Header":
                        skip_column_names = True

            elif section == "Probes":
                probes[splitline[2]] = (len(splitline[9]), len(splitline[11]))

            elif section == "Targets":
                amplicon_name = "{0}:{1}-{2}".format(splitline[3], splitline[4], splitline[5])
                if amplicon_name not in excluded and splitline[2] == "1":
                #    target = gr1
                #    rename_offtarget[splitline[0]] = amplicon_name
                #else:
                #    target = gr2
                #    amplicon_name = splitline[0]
                #target.add( [splitline[3],
                    gr1.construct( [splitline[3],
                        int(splitline[4])+probes[splitline[0]][splitline[6]=="+"],
                        int(splitline[5])-probes[splitline[0]][splitline[6]=="-"],
                        amplicon_name,
                        splitline[6],
                        "",
                        1] )#CHANGED GR
        #for chr_name in gr2:
        #    for entry in gr2[chr_name]:
        #        if entry[NAME] in rename_offtarget:
        #            entry[NAME] = rename_offtarget[entry[NAME]]
    gr1.sort()
    #gr2.sort()
    #return (gr1, gr2)
    return gr1

# 
# 
# 
# 
# 
def _simplify_gene_names(gr):
    names = {}
    duplicates = set([])
    for entry in gr.all_entries:
        gene, transcript = entry[NAME].split(" ")[0:2]
        if gene in names and names[gene] != transcript:
            duplicates.add(gene)
        else:
            names[gene] = transcript
    for entry in gr.all_entries:
        splitgene = entry[NAME].split(" ")
        if splitgene[0] not in duplicates:
            if len(splitgene) <= 2:
                entry[NAME] = splitgene[0]
            else:
                entry[NAME] = splitgene[0]+" "+" ".join(splitgene[2:])


def load_refflat(path, genes_of_interest, canonical_transcripts):
    canonicaldict = {}
    if canonical_transcripts != "": 
        with file(canonical_transcripts, "rU") as f:
            for line in f:
                gene, transcript = line.rstrip().split("\t")
                canonicaldict[gene] = transcript.split()[1]

    needed = set([])
    doublecheck = {}
    include_everything = True
    if genes_of_interest != "":
        with file(genes_of_interest, "rU") as f:
            for line in f:
                line = line.strip()
                splitline = line.split()
                if line != "":
                    include_everything = False
                    if line in canonicaldict:
                        needed.add(canonicaldict[line])
                        doublecheck[canonicaldict[line]] = line 
                    elif len(splitline) == 1:
                        needed.add(line)
                    elif len(splitline) == 2:
                        needed.add(splitline[1])
                        doublecheck[splitline[1]] = splitline[0]
                    else:
                        raise #CoverMiException("Malformed line in gene file: {0}".format(line))
            
    with file(path, "rU") as f:
        multiple_copies = {}
        found = set([])
        exons = GenomicRange()
        transcripts = GenomicRange()
        for line in f:
            splitline = line.rstrip("\n").split("\t")
            if splitline[2] in KARYOTYPE:
                if include_everything:
                    if splitline[0] in canonicaldict and splitline[1] != canonicaldict[splitline[0]]:
                        continue
                elif splitline[0] in needed:
                    found.add(splitline[0])
                elif splitline[1] in needed and splitline[0] == doublecheck[splitline[1]]:
                    found.add(splitline[1])
                else:
                    continue

                name = "{0} {1}".format(splitline[0], splitline[1])
                if name not in multiple_copies:
                    multiple_copies[name] = 1
                else:
                    if multiple_copies[name] == 1:
                        for gr in (exons, transcripts):
                            for entry in gr.all_entries:
                                if entry[NAME] == name:
                                    entry[NAME] = "{0} copy 1".format(name)
                    multiple_copies[name] += 1
                    name = "{0} copy {1}".format(name, multiple_copies[name])
            
                transcripts.construct( [splitline[2], int(splitline[4])+1-SPLICE_SITE_BUFFER, int(splitline[5])+SPLICE_SITE_BUFFER, name, splitline[3], "", 1] )

                exon_numbers = range(1,int(splitline[8])+1) if (splitline[3] == "+") else range(int(splitline[8]),0,-1)            
                for start, stop, exon in zip(splitline[9].rstrip(",").split(","), splitline[10].rstrip(",").split(","), exon_numbers):
                    exons.construct( [splitline[2], int(start)+1-SPLICE_SITE_BUFFER, int(stop)+SPLICE_SITE_BUFFER, name, splitline[3], exon, 1] )

    missing = needed - found
    if not include_everything and len(missing) > 0:
        missing = ["{0} {1}".format(doublecheck[trans], trans) if trans in doublecheck else trans for trans in missing]
        print "WARNING - {0} not found in reference file".format(", ".join(missing))        

    for gr in (exons, transcripts):
        gr.sort()
        _simplify_gene_names(gr)
    return (exons, transcripts)


def load_variants(path, catagory, genes_of_interest=None, disease_names=None): ######################################################################################################
    if genes_of_interest is not None:
        genes_of_interest = set([name.split()[0] for name in genes_of_interest.names])

    disease_dict = {}
    if disease_names is not None:
        with file(disease_names, "rU") as f:
            for line in f:
                if line[0] != "#":
                    line2 = line.rstrip().split("=")
                    if len(line2) == 2:
                        disease_dict[line2[0].strip("\" ")] = line2[1].strip("\" ")

    with file(path, "rU") as f:

        mutations = {}
        for line in f:
            splitline = line.rstrip().split("\t")
            if splitline == [""] or splitline[4] not in KARYOTYPE:
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
                mutations[mutation] = [splitline[4], int(splitline[5]), int(splitline[6]), name, splitline[7], set([]) if (catagory=="mutation") else "", 1]
            else:
                mutations[mutation][WEIGHT] += 1
            if catagory == "mutation":
                mutations[mutation][NAME2].add(disease)

    gr1 = GenomicRange()
    for entry in mutations.values():
        if catagory == "mutation":
            entry[NAME2] = "; ".join(sorted(entry[NAME2]))
        gr1.construct(entry)
    gr1.sort()
    return gr1


def load(path):
    with file(path, "rU") as f:
        gr1 = GenomicRange()
        for line in f:
            splitline = line.rstrip("\n").split("\t")
            if splitline != [""]:
                name1, name2, strand = ("", "", "")
                if len(splitline) >= 4:
                    names = splitline[3].split("\\")
                    name1 = names[0]
                    if len(names) > 1:
                        name2 = names[1]
                    if len(splitline) >= 6:
                        strand = splitline[5]
                gr1.construct( [splitline[0], int(splitline[1])+1, int(splitline[2]), name1, strand, name2, 1] )############################## DOES NOT LOAD/SAVE WEIGHT
    gr1.sort()
    return gr1




def new(entry=None):#############################################?add new constructor
    return GenomicRange() if (entry is None) else GenomicRange().construct(list(entry))




class GenomicRange(dict):


    def construct(self, entry): # [chr_name, start, stop, name1, strand, exon/diseases, weight]
        assert len(entry) == 7, "Genomic Range entry of incorrect length"
        if entry[CHROM] not in self:
            self[entry[CHROM]] = []
        self[entry[CHROM]].append(entry)
        return self
            

    def sort(self):
        for chr_name in self:
            self[chr_name].sort()


    @property
    def all_entries(self):
        for chr_name in sorted(self, key=chromosome_number):
            for entry in self[chr_name]:
                yield entry

    @property
    def merged(self): # The name of a range will be the name of the first sub-range that was merged into it
        merged = GenomicRange()
        for chr_name in self:
            merged[chr_name] = [list(self[chr_name][0])]
            for entry in self[chr_name]:
                if entry[START]-1 <= merged[chr_name][-1][STOP]:
                    merged[chr_name][-1][STOP] = max(entry[STOP], merged[chr_name][-1][STOP])
                else:
                    merged.construct(list(entry))
        merged.sort()
        return merged


    def extended_to_include_touching_amplicons(self, gr):
         extended = GenomicRange() 
         amplicons = GenomicRange()
         for chr_name in self:
             for a in self[chr_name]:
                 new_entry = list(a)
                 if chr_name in gr:
                     for b in gr[chr_name]:
                         if b[START]>a[STOP]:
                             break
                         if (a[START]<=b[START]<=a[STOP]) or (a[START]<=b[STOP]<=a[STOP]) or b[STOP]>a[STOP]:
                             new_entry[START] = min(new_entry[START], b[START])
                             new_entry[STOP] = max(new_entry[STOP], b[STOP])
                             amplicons.construct(list(b))
                 extended.construct(new_entry)
         return (extended, amplicons)


    def overlapped_by(self, gr):# If range gr contains overlapping ranges then we may get multiple copies of the overlapping ranges
         trimmed = GenomicRange()
         for chr_name in self:
            if chr_name in gr:
                for a in self[chr_name]:
                    for b in gr[chr_name]:
                        if b[START]>a[STOP]:
                            break
                        if (a[START]<=b[START]<=a[STOP]) or (a[START]<=b[STOP]<=a[STOP]) or b[STOP]>a[STOP]:
                            new_entry = list(a)
                            new_entry[START] = max(a[START],b[START])
                            new_entry[STOP] = min(a[STOP],b[STOP])
                            trimmed.construct(new_entry)
         trimmed.sort() # Only needed if gr contains overlapping ranges which should not happen with sane use
         return trimmed


    def not_touched_by(self, gr):
        untouched = GenomicRange() 
        for chr_name in self:
           for a in self[chr_name]:
               touching = False
               if chr_name in gr:
                   for b in gr[chr_name]:
                       if b[START]>a[STOP]:
                           break
                       if (a[START]<=b[START]<=a[STOP]) or (a[START]<=b[STOP]<=a[STOP]) or b[STOP]>a[STOP]:
                           touching = True
                           break
               if not touching:
                   untouched.construct(list(a))
        untouched.sort()
        return untouched


    def touched_by(self, gr):
         touched = GenomicRange() 
         for chr_name in self:
            for a in self[chr_name]:
                if chr_name in gr:
                    for b in gr[chr_name]:
                        if b[START]>a[STOP]:
                            break
                        if (a[START]<=b[START]<=a[STOP]) or (a[START]<=b[STOP]<=a[STOP]) or b[STOP]>a[STOP]:
                            touched.construct(list(a))
                            break
         touched.sort()
         return touched

    
    def subranges_covered_by(self, gr):
         covered = GenomicRange()
         for chr_name in self:
            if chr_name in gr:
                for a in self[chr_name]:
                    for b in gr[chr_name]:
                        if b[START]>a[START]:
                            break
                        if b[STOP]>=a[STOP]:
                            covered.construct(list(a))
                            break
         covered.sort()
         return covered 


    def combined_with(self, gr):
        combined = GenomicRange()
        for chr_name in self:
            for entry in self[chr_name]:
                combined.construct(list(entry))
        for chr_name in gr:
            for entry in gr[chr_name]:
                combined.construct(list(entry))
        combined.sort()
        return combined


    def subset2(self, names, exclude=False, genenames=False):
        if type(names) != set:
            names = set(names) if (type(names)==list) else set([names])
        if genenames:
            names = set([name.split()[0] for name in names])
        gr = GenomicRange()
        for entry in self.all_entries:
            name = entry[NAME].split()[0] if genenames else entry[NAME]
            if (name in names) ^ exclude:
                gr.construct(list(entry))
        return gr


    @property
    def names(self):
        names = set([])
        for entry in self.all_entries:
            names.add(entry[NAME])
        return sorted(names)


    @property
    def number_of_components(self):
        output = 0
        for chr_name in self:
            output += len(self[chr_name])
        return output


    @property
    def number_of_weighted_components(self):
        output = 0
        for entry in self.all_entries:
            output += entry[WEIGHT]
        return output


    @property
    def base_count(self):
        output = 0
        for entry in self.all_entries:
            output += entry[STOP] - entry[START] + 1
        return output


    @property
    def locations_as_string(self):
        output = []
        for entry in self.all_entries:
            output.append("{0}:{1}-{2}".format(entry[CHROM], entry[START], entry[STOP]))
        return ", ".join(output)


    @property
    def names_as_string(self):
        namedict = {}
        for entry in self.all_entries:
            if entry[NAME] not in namedict:
                namedict[entry[NAME]] = []
            if type(entry[NAME2]) == int:
                namedict[entry[NAME]].append(entry[NAME2])
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
        return ", ".join(namelist).strip()


    @property
    def is_empty(self):
        return len(self) == 0


    def save(self, f): #Save GenomicRange object in bedfile format, START POSITION IS ZERO BASED
        if type(f) == str:
            with file(f, "wt") as f2:
                self._save(f2)
        else:
            self._save(f)

    def _save(self, f): #Save GenomicRange object in bedfile format, START POSITION IS ZERO BASED
        for chr_name in self:
            for entry in self[chr_name]:
                f.write("{0}\t{1}\t{2}\t{3}\t.\t{4}\n".format(entry[CHROM], entry[START]-1, entry[STOP], "{0}\\{1}".format(entry[NAME], entry[NAME2]), entry[STRAND]))


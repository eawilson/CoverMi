# { chr_name : [ [start, stop, depth], [start, stop, depth] ] }

import re, subprocess, os, tempfile, copy, pdb
from gr import Gr


class CoverMiException(Exception):
    pass


CSTART = 0
CSTOP = 1
CDEPTH = 2
FORWARD = True
REVERSE = False

NAME = 0
DEPTH = 1
BASES = 2
RANGE = 3
COMPONENTS = 4
DISEASES = 5
WEIGHTED_COMPONENTS = 6
COVERED = True
UNCOVERED = False
        
AM_NAME = 0
AM_FDEPTH = 1
AM_RDEPTH = 2
AM_GR = 3


class Cov(dict):


    @staticmethod
    def _amplicon_number(name):
        chrom, start = re.match("(chr[0-9XYM]+):([0-9]+)-", name[0]).groups()[0:2]
        return (Gr.KARYOTYPE[chrom]*Gr.MAX_CHR_LENGTH)+int(start)


    @classmethod
    def load_bam(coverage, bam, gr, amplicons=False):
        return coverage._load_bam_amplicons(bam, gr) if amplicons else  coverage._load_bam_exons(bam, gr)
        

    @classmethod
    def _load_bam_exons(coverage, bamfile_name, exons):
        cov = coverage()
        bedfile_temp = tempfile.TemporaryFile()
        exons.merged.save(bedfile_temp)
        bedfile_temp.seek(0, 0)

        bedtools = coverage._create_bedtools(bamfile_name, bedfile_temp, stranded=False)
        
        first_pass = True
        last_chr_name = ""
        while True:
            splitline = bedtools.stdout.readline().rstrip("\n").split("\t")
            if splitline == [""]:
                if not first_pass:
                    cov[last_chr_name].append([first_pos, last_pos, last_depth])
                break
            pos = int(splitline[1])+int(splitline[6])
            depth = int(splitline[7])
            chr_name = splitline[0]                 
            if chr_name != last_chr_name or pos != last_pos + 1 or depth != last_depth:
                if not first_pass:
                    if last_chr_name not in cov:
                        cov[last_chr_name] = []
                    cov[last_chr_name].append([first_pos, last_pos, last_depth])
                first_pass = False   
                first_pos = pos
            last_pos = pos
            last_depth = depth
            last_chr_name = chr_name                 
            
        bedfile_temp.close()
        
        for chr_name in cov:
            current = cov[chr_name]
            current.sort()
            current.append([current[-1][CSTOP]+1, Gr.MAX_CHR_LENGTH, 0])
            current.append([0, current[0][CSTART]-1, 0])
            for index in range(0, len(current)-3):
                if current[index][CSTOP] != current[index+1][CSTART]-1:
                    current.append([current[index][CSTOP]+1, current[index+1][CSTART]-1, 0])
            current.sort()
        return cov
        

#    @classmethod
#    def _load_bam_exons(coverage, bamfile_name, exons):
#        cov = coverage()
#        bedfile_temp = tempfile.TemporaryFile()
#        exons.merged.save(bedfile_temp)
#        bedfile_temp.seek(0, 0)

#        bedtools = coverage._create_bedtools(bamfile_name, bedfile_temp, stranded=False)

#        name = ""
#        while True:
#            splitline = bedtools.stdout.readline().rstrip("\n").split("\t")
#            if splitline == [""] or splitline[3] != name:
#                if name != "":
#                    cov[chr_name].append([runstartpos, pos, depth])
#                if splitline == [""]:
#                    break
#                chr_name = splitline[0]
#                if chr_name not in cov:
#                    cov[chr_name] = []
#                runstartpos = int(splitline[1])+int(splitline[6])
#                runstartdepth = int(splitline[7])
#                cov[chr_name].append([runstartpos-1, runstartpos-1, 0])
#        
#            pos = int(splitline[1])+int(splitline[6])
#            depth = int(splitline[7])
#            if runstartdepth != depth:
#                cov[chr_name].append([runstartpos, pos-1, runstartdepth])
#                runstartpos = pos
#                runstartdepth = depth
#            name = splitline[3]
#        bedfile_temp.close()
#            
#        for chr_name in cov:
#            cov[chr_name].sort()
#    
#        for chr_name in cov:
#            cov[chr_name][0][CSTART] = 0
#            for index in range(1, len(cov[chr_name])):
#                #if cov[chr_name][index][CSTART] != (cov[chr_name][index-1][CSTOP]+1) and (cov[chr_name][index][CSTART] != cov[chr_name][index][CSTOP] or cov[chr_name][index][CDEPTH] != 0):
#                #    print cov[chr_name][index]###########################################????????????????????????????????????????????????????????????????????????????
#                cov[chr_name][index][CSTART] = cov[chr_name][index-1][CSTOP]+1
#            cov[chr_name].append([cov[chr_name][-1][CSTOP]+1, Gr.MAX_CHR_LENGTH, 0])
#        return cov


    @classmethod
    def _load_bam_amplicons(coverage, bamfile_name, amplicons):
        cov = coverage()

        amplicon_metrics = {}
        for chr_name in amplicons:
            for entry in amplicons[chr_name]:
                amplicon_metrics[entry.name] = Amplicon_info([entry.name, "NA", "NA", Gr(entry)])

        bedfile_temp = tempfile.TemporaryFile()
        for chr_name in sorted(amplicons):
            if len(amplicons[chr_name]) == 1:
                bedfile_temp.write("{0}\t{2}\t{3}\t{1}\t.\t+\n{0}\t{2}\t{3}\t{1}\t.\t-\n".format(chr_name, amplicons[chr_name][0].name,
                    amplicons[chr_name][0].start-1, 
                    amplicons[chr_name][0].stop))

            else:
                bedfile_temp.write("{0}\t{2}\t{3}\t{1}\t.\t+\n{0}\t{4}\t{5}\t{1}\t.\t-\n".format(chr_name, amplicons[chr_name][0].name,
                    amplicons[chr_name][0].start-1, 
                    min(amplicons[chr_name][0].stop, amplicons[chr_name][1].start-1),
                    amplicons[chr_name][0].start-1, 
                    amplicons[chr_name][0].stop))

                for index in range(1, len(amplicons[chr_name])-1):
                    bedfile_temp.write("{0}\t{2}\t{3}\t{1}\t.\t+\n{0}\t{4}\t{5}\t{1}\t.\t-\n".format(chr_name, amplicons[chr_name][index].name,
                        amplicons[chr_name][index].start-1,
                        min(amplicons[chr_name][index].stop, amplicons[chr_name][index+1].start-1),
                        max(amplicons[chr_name][index].start, amplicons[chr_name][index-1].stop+1)-1, 
                        amplicons[chr_name][index].stop))

                bedfile_temp.write("{0}\t{2}\t{3}\t{1}\t.\t+\n{0}\t{4}\t{5}\t{1}\t.\t-\n".format(chr_name, amplicons[chr_name][-1].name,
                    amplicons[chr_name][-1].start-1, 
                    amplicons[chr_name][-1].stop,
                    max(amplicons[chr_name][-1].start, amplicons[chr_name][-2].stop+1)-1, 
                    amplicons[chr_name][-1].stop))
        bedfile_temp.seek(0, 0)

        bedtools = coverage._create_bedtools(bamfile_name, bedfile_temp, stranded=True)

        # loop one through bedtools output and transcribe positions and depths into "firstpass" 
        firstpass = {}# { "chr1" : [ [pos, depth], [pos, depth], [pos, depth], ...] }
        name = "" 
        while True:
            splitline = bedtools.stdout.readline().rstrip("\n").split("\t")
            if splitline == [""] or splitline[3] != name or splitline[5] != strand:
                if name != "":
                    if chr_name not in firstpass:
                        firstpass[chr_name] = []
                    firstpass[chr_name] += tempdata
                    if amplicon_metrics[name][(strand=="-")+1] != "NA":  
                        raise CoverMiException("Duplicate amplicon {0}".format(name))
                    amplicon_metrics[name][(strand=="-")+1] = tempdata[len(tempdata)/2][1] if (len(tempdata)>0) else 0
                        
                if splitline == [""]:
                    break
                tempdata = []
                chr_name = splitline[0]
        
            depth = int(splitline[7])
            strand = splitline[5]
            name = splitline[3]
            if depth > 0: 
                tempdata.append([int(splitline[1])+int(splitline[6]), depth])
        bedfile_temp.close()

        cov = coverage._firstpass_into_coverage(firstpass)
        cov.amplicon_info = amplicon_metrics.values()
        try:
            cov.amplicon_info.sort(key=coverage._amplicon_number)
        except AttributeError:
            pass
        return cov


    @staticmethod
    def _create_bedtools(bamfile_name, bedfile_temp, stranded=False):
        if os.name == "nt":
            bamfile_name = os.path.abspath(bamfile_name).replace("C:", "/cygdrive/c").replace("\\", "/")
            if os.path.exists("C:\\cygwin64\\bin\\"):
                cygwin_path = "C:\\cygwin64\\bin\\"
            elif os.path.exists("C:\\cygwin\\bin\\"):
                cygwin_path = "C:\\cygwin\\bin\\"
            else:
                raise CoverMiException("Unable to find Cygwin path")
            def decorate(command): return "START /B {0}bash.exe -c \"PATH=/usr/bin:/usr/local/bin:$PATH; {1}; exit\"".format(cygwin_path, command)
        else:
            def decorate(command): return command
        version = subprocess.Popen(decorate("bedtools --version"), stdout=subprocess.PIPE, universal_newlines=True, shell=True).stdout.readline().strip()
        if not version.startswith("bedtools v"):
            raise CoverMiException("Unable to start bedtools")
        major, minor, patch =  [int(num) for num in version[10:].split(".")]
        if major > 2 or (major == 2 and minor >=24):
            command = "bedtools coverage -d {0} -a '{2}' -b '{1}'"
        else:
            command = "bedtools coverage -d {0} -abam '{1}' -b '{2}'"
            
        command = decorate(command.format("-s" if stranded else "", bamfile_name, "stdin"))
        return subprocess.Popen(command, stdout=subprocess.PIPE, stdin=bedfile_temp, universal_newlines=True, shell=True)


    @classmethod
    def _firstpass_into_coverage(coverage, fp):
        cov = coverage()

        for chr_name in fp:
            fp[chr_name].sort()
            cov[chr_name] =[[0, 0, 0]]
            lastpos = 0

            for index in range(0,len(fp[chr_name])):
                if index < len(fp[chr_name])-2 and fp[chr_name][index][0] == fp[chr_name][index+1][0]:#if depth exists for adjacent + and - reads then combine
                    index += 1
                    fp[chr_name][index][1] += fp[chr_name][index-1][1]

                if fp[chr_name][index][0] > lastpos+1 and cov[chr_name][-1][CDEPTH] > 0:
                    cov[chr_name].append([lastpos+1, 0, 0])
                    
                if fp[chr_name][index][1] <> cov[chr_name][-1][CDEPTH]:
                    cov[chr_name].append([fp[chr_name][index][0], 0, fp[chr_name][index][1]])
                lastpos=fp[chr_name][index][0]

            if cov[chr_name][-1][CDEPTH] > 0:
                cov[chr_name].append([lastpos+1, 0, 0])
            cov[chr_name].append([Gr.MAX_CHR_LENGTH, Gr.MAX_CHR_LENGTH, 0])

        for chr_name in cov: #Loop over coverage object start positions and add in stop positions
            for index in range(0, len(cov[chr_name])-1):
                cov[chr_name][index][CSTOP] = cov[chr_name][index+1][CSTART]-1
        return cov


    @classmethod
    def perfect_coverage(coverage, gr1, perfect_depth=1000):
        firstpass = {}
        for chr_name in gr1:
            firstpass[chr_name] = []
            for entry in gr1[chr_name]:
                read_length = (entry.stop - entry.start + 1) * 2 / 3
                for pos in range(entry.start, entry.start+read_length+1):
                    firstpass[chr_name].append([pos, perfect_depth])
                for pos in range(entry.stop-read_length, entry.stop+1):
                    firstpass[chr_name].append([pos, perfect_depth])
        return coverage._firstpass_into_coverage(firstpass)


    @classmethod
    def load(coverage, path):
        cov = coverage()
        with file(path, "rU") as f:
            for row in f:
                chr_name, cstart, cstop, depth = row.rstrip("\n").split("\t")
                if chr_name not in cov:
                    cov[chr_name] = []
                cov[chr_name].append([int(cstart), int(cstop) ,int(depth)])
        return cov


    def __getitem__(self, chr_name):
        try:
            return super(type(self), self).__getitem__(chr_name)
        except KeyError:
            return [[0, Gr.MAX_CHR_LENGTH-1, 0], [Gr.MAX_CHR_LENGTH, Gr.MAX_CHR_LENGTH, 0]]


    def calculate(self, gr1, min_depth, exons=False, total=False):
        results = []
        resultkey = {}
        for entry in gr1.all_entries:
            name = entry.name if (not total) else ""
            name = name if (not exons) else "{0} e{1}".format(name, entry.exon)
            if name not in resultkey:
                resultkey[name] = len(results)
                results.append(Coverage_info([ name, [0, 0], [0, 0], [Gr(), Gr()], [0, 0], None, [0, 0] ]))
                try:
                    results[resultkey[name]][DISEASES] = entry.diseases
                except AttributeError:
                    pass
            line = results[resultkey[name]]

            allcovered = True
            for cstart, cstop, cdepth in self[entry.chrom]:
                if cstart > entry.stop:
                    break
                elif cstop >= entry.start:
                    bases = min(entry.stop, cstop) - max(entry.start, cstart) + 1
                    line[DEPTH][cdepth>=min_depth] += (bases*cdepth)
                    line[BASES][cdepth>=min_depth] += bases
                    allcovered = allcovered and (cdepth>=min_depth)
                    if entry.chrom in (line[RANGE][cdepth>=min_depth]) and cstart == line[RANGE][cdepth>=min_depth][entry.chrom][-1].stop+1:
                            line[RANGE][cdepth>=min_depth][entry.chrom][-1].stop = min(entry.stop, cstop)
                    else:
                        line[RANGE][cdepth>=min_depth].construct(copy.copy(entry)) 
                        line[RANGE][cdepth>=min_depth][entry.chrom][-1].start = max(entry.start, cstart) 
                        line[RANGE][cdepth>=min_depth][entry.chrom][-1].stop = min(entry.stop, cstop)
            line[COMPONENTS][COVERED] += int(allcovered)
            line[COMPONENTS][UNCOVERED] += int(not(allcovered))
            try:
                line[WEIGHTED_COMPONENTS][COVERED] += int(allcovered) * entry.weight
                line[WEIGHTED_COMPONENTS][UNCOVERED] += int(not(allcovered)) * entry.weight
            except AttributeError:
                pass

        for result in results:
            for index in [0,1]:
                result[DEPTH][index] /= max(result[BASES][index], 1)
        results.sort()
        return results if (not total) else results[0]


    def save_plot(self, gr1, filename): #Replaced by decorated plots but may be useful for debugging plotting
        with file(filename, "wt") as f:
            f.write("chrom\tpos\tdepth\tname\n")
            for chr_name in gr1:
                for entry in gr1[chr_name]:
                    f.write("{0}\t{1}\t{2}\t{3}\n".format(entry.chrom, entry.start-1, 0, entry.name))
                    for cstart, cstop, cdepth in self[chr_name]:
                        if cstop >= entry.start:
                            f.write("{0}\t{1}\t{2}\t{3}\n".format(entry.chrom, max(entry.start, cstart), cdepth, entry.name))
                        if cstop >= entry.stop:
                            if cstart < entry.stop:
                                f.write("{0}\t{1}\t{2}\t{3}\n".format(entry.chrom, entry.stop, cdepth, entry.name))
                            break
                    f.write("{0}\t{1}\t{2}\t{3}\n".format(entry.chrom, entry.stop+1, 0, entry.name))


    def filter_by_germline(self, germ, germmin):
        newcov = type(self)()

        for chrom in self:

            newcov[chrom] = []
            selfindex = 0
            germindex = 0
            start = 0
            current = 0
            depth = 0

            selflist = self[chrom]
            germlist = germ[chrom]
            while True:
  
                newdepth = 0 if (germlist[germindex][CDEPTH] < germmin) else selflist[selfindex][CDEPTH]
                if depth != newdepth:
                    newcov[chrom].append([start, current-1, depth])
                    start = current
                    depth = newdepth

                if germlist[germindex][CSTOP] < selflist[selfindex][CSTOP]:
                    germindex += 1
                    current = germlist[germindex][CSTART]
                elif germlist[germindex][CSTOP] > selflist[selfindex][CSTOP]:
                    selfindex += 1
                    current = selflist[selfindex][CSTART]
                elif len(selflist) > selfindex+1:
                    selfindex += 1
                    germindex += 1
                    current = selflist[selfindex][CSTART]
                else:
                    newcov[chrom].append([start, selflist[selfindex][CSTOP], depth])
                    break
        return newcov


    def save(self, path):
        with file(path, "wt") as f:
            for chr_name in self:
                for entry in self[chr_name]:
                    f.write("{0}\t{1}\t{2}\t{3}\n".format(chr_name, entry[CSTART], entry[CSTOP], entry[CDEPTH]))


class Amplicon_info(list):
    @property
    def name(self):
        return self[AM_NAME]
    @property
    def chrom(self):
        return self[AM_NAME].split(":")[0]
    @property
    def f_depth(self):
        return self[AM_FDEPTH]
    @property
    def r_depth(self):
        return self[AM_RDEPTH]
    @property
    def minimum_depth(self):
        return self[AM_FDEPTH] if self[AM_FDEPTH]<self[AM_RDEPTH] else self[AM_RDEPTH]
    @property
    def maximum_depth(self):
        return self[AM_RDEPTH] if self[AM_FDEPTH]<self[AM_RDEPTH] else self[AM_FDEPTH]
    @property
    def mean_depth(self):
        return (self[AM_FDEPTH]+self[AM_RDEPTH])/2
    @property
    def ratio(self):
        return float(min(self[AM_FDEPTH], self[AM_RDEPTH]))/max(self[AM_FDEPTH], self[AM_RDEPTH], 1)
    @property
    def gr(self):
        return self[AM_GR]


class Coverage_info(list):
    @property
    def name(self):
        return self[NAME]
    @property
    def diseases(self):
        return self[DISEASES]
    @property
    def percent_covered(self):
        return float(self[BASES][COVERED]*100)/sum(self[BASES])
    @property
    def percent_uncovered(self):
        return float(self[BASES][UNCOVERED]*100)/sum(self[BASES])
    @property
    def depth_covered(self):
        return self[DEPTH][COVERED]
    @property
    def depth_uncovered(self):
        return self[DEPTH][UNCOVERED]
    @property
    def bases_covered(self):
        return self[BASES][COVERED]
    @property
    def bases_uncovered(self):
        return self[BASES][UNCOVERED]
    @property
    def range_covered(self):
        return self[RANGE][COVERED]
    @property
    def range_uncovered(self):
        return self[RANGE][UNCOVERED]
    @property
    def range_combined(self):
        return self[RANGE][COVERED].combined_with(self[RANGE][UNCOVERED]).merged
    @property
    def bases(self):
        return sum(self[BASES])
    @property
    def components_covered(self):
        return self[COMPONENTS][COVERED]
    @property
    def components_uncovered(self):
        return self[COMPONENTS][UNCOVERED]
    @property
    def components(self):
        return sum(self[COMPONENTS])
    @property
    def percent_components_covered(self):
        return float(self[COMPONENTS][COVERED]*100)/sum(self[COMPONENTS])
    @property
    def percent_components_uncovered(self):
        return float(self[COMPONENTS][UNCOVERED]*100)/sum(self[COMPONENTS])
    @property
    def weighted_components_covered(self):
        return self[WEIGHTED_COMPONENTS][COVERED]
    @property
    def weighted_components_uncovered(self):
        return self[WEIGHTED_COMPONENTS][UNCOVERED]
    @property
    def weighted_components(self):
        return sum(self[WEIGHTED_COMPONENTS])
    @property
    def percent_weighted_components_covered(self):
        return float(self[WEIGHTED_COMPONENTS][COVERED]*100)/sum(self[WEIGHTED_COMPONENTS])
    @property
    def percent_weighted_components_uncovered(self):
        return float(self[WEIGHTED_COMPONENTS][UNCOVERED]*100)/sum(self[WEIGHTED_COMPONENTS])
    @property
    def completely_covered(self):
        return self[BASES][UNCOVERED]==0
    @property
    def incompletely_covered(self):
        return self[BASES][UNCOVERED]>0



from __future__ import print_function, absolute_import, division

import sys, os, time, re, getopt, pdb
from collections import namedtuple

from .cov import Cov, fake_paired_end_reads
from .panel import Panel
from . import technicalreport, clinicalreport, designreport, covermiplot, tablereport
from .include import CoverMiException, __version__



illumina_suffix = re.compile("(.+)_S[0-9]+$")
Sample = namedtuple("Sample", "name run path")

def walk_samples(path, ext=".bam", strip_trailing_underscore=True):
    found = []

    if os.path.isfile(path):
        if not os.path.splitext(path)[1] == ext:
            raise CoverMiException("file "+path+" is not of type "+ext)
        path_noext, extension = os.path.splitext(path)
        runpath, sample = os.path.split(path_noext)
        run = os.path.basename(runpath)
        yield Sample(strip_underscore(sample), run, path)

    elif os.path.isdir(path):
        for fn in os.listdir(path):
            if os.path.isfile(os.path.join(path, fn)) and os.path.splitext(fn)[1] == ext: # Single run format
                for sample in walk_run(path, ext):
                    yield sample
                return

        found = [] # Muliptle run format
        for fn in os.listdir(path):
            if os.path.isdir(os.path.join(path, fn)):
               for sample in walk_run(os.path.join(path, fn), ext):
                    yield sample

    else:
        raise CoverMiException(path+" is not a file or directory")



def walk_run(path, ext):
    run = os.path.basename(path)
    for root, dirnames, filenames in os.walk(path):
        for fn in filenames:
            sample, extension = os.path.splitext(fn)
            if os.path.isfile(os.path.join(root, fn)) and extension == ext:
                yield Sample(strip_underscore(sample), run, os.path.join(root, fn))
        


def strip_underscore(sample):
    matchobj = illumina_suffix.search(sample)
    return matchobj.group(1) if matchobj else sample 



def create_output_dir(output_path, bam_path):
    output_path = os.path.join(output_path, "{0}_covermi_output".format(os.path.splitext(os.path.basename(bam_path))[0]))
    try:
        os.mkdir(output_path)
    except OSError:
        if os.path.isdir(output_path):
            raise CoverMiException("{0} folder already exists".format(output_path))
        else:
            raise CoverMiException("Unable to create folder {0}".format(output_path))
    return output_path


def covermimain(panel_path, bam_path, output_path, depth=None):
    print("CoverMi v{} (Python {}.{}.{})".format(__version__, *sys.version_info[:3]))

    panel = Panel(panel_path, verbose=True, splice_site_buffer=5)
    if depth is not None:
        panel.properties["depth"] = depth
    output_path = create_output_dir(output_path, bam_path if bam_path!="" else panel_path)
    print("Processing...")

    if bam_path != "":
        bam_file_list = list(sample for sample in walk_samples(bam_path))
        if len(bam_file_list) == 1:
            clinical_report_path = output_path
            technical_report_path = output_path
        else:
            clinical_report_path = os.path.join(output_path, "clinical")
            technical_report_path = os.path.join(output_path, "technical")
            os.mkdir(clinical_report_path)
            os.mkdir(technical_report_path)

        output_stems = set([])
        for sample in bam_file_list:
            start_time = time.time()
            print("{0}/{1}".format(sample.run, sample.name))

            output_stem = sample.name
            dup_num = 1
            while output_stem in output_stems:
                output_stem = "{0}({1})".format(sample.name, dup_num)
                dup_num += 1
            output_stems.add(output_stem)

            extra = {"amplicons": panel.amplicons} if "amplicons" in panel else {"locations": panel.transcripts, "print_progress": True}
            cov = Cov(sample.path, **extra)

            if cov.amplicon_info: technicalreport.create(cov.amplicon_info, panel, sample, os.path.join(technical_report_path, output_stem))
            clinicalreport.create(cov, panel, sample, os.path.join(clinical_report_path, output_stem))
            #tablereport.create(cov, panel, sample, os.path.join(clinical_report_path, output_stem))
            covermiplot.plot(cov, panel, sample.name, os.path.join(clinical_report_path, output_stem+"_plot.pdf"))

            seconds = int(time.time() - start_time)
            time_string = "{0} sec".format(seconds) if (seconds<60) else "{0} min {01} sec".format(seconds//60, seconds%60)
            print("file {0} of {1} completed in {2}".format(len(output_stems), len(bam_file_list), time_string))
    else:
        cov = Cov(fake_paired_end_reads(panel.amplicons, depth=panel.depth*2))
        designreport.create(cov, panel, None, os.path.join(output_path, panel.name))
        covermiplot.plot(cov, panel, panel.name, os.path.join(output_path, panel.name+"_plot.pdf"))


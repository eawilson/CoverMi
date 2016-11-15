#!/usr/bin env python
import sys, os, time, re, getopt, pdb
from cov import Cov
from panel import Panel
from files import Files
import covermiconf, technicalreport, clinicalreport, designreport, covermiplot
#import testreport as clinicalreport

class CoverMiException(Exception):
    pass


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


def main(panel_path, bam_path, output_path, depth=None):

    if not os.path.exists(panel_path):
        panel_path = os.path.join(covermiconf.load_conf["panel_root"], panel_path)

    panel = Panel(panel_path).load(bam_path=="")
    if depth is not None:
        panel["Options"]["Depth"] = int(depth)
    output_path = create_output_dir(output_path, bam_path if bam_path!="" else panel_path)
    print "Processing..."

    if bam_path != "":
        bam_file_list = Files(bam_path, ".bam")
        if len(bam_file_list) == 1:
            clinical_report_path = output_path
            technical_report_path = output_path
        else:
            clinical_report_path = os.path.join(output_path, "clinical")
            technical_report_path = os.path.join(output_path, "technical")
            os.mkdir(clinical_report_path)
            os.mkdir(technical_report_path)

        output_stems = set([])
        for panel["Filenames"]["Run"], panel["Filenames"]["Sample"], path in bam_file_list:
            path += ".bam"
            start_time = time.time()
            print "{0}/{1}".format(panel["Filenames"]["Run"], panel["Filenames"]["Sample"])

            output_stem = panel["Filenames"]["Sample"]
            dup_num = 1
            while output_stem in output_stems:
                output_stem = "{0}({1})".format(panel["Filenames"]["Sample"], dup_num)
                dup_num += 1
            output_stems.add(output_stem)
            
            if "Amplicons" in panel:
                cov = Cov.load_bam(path, panel["Amplicons"], amplicons=True)
                technicalreport.create(cov.amplicon_info, panel, os.path.join(technical_report_path, output_stem))
            else:
                cov = Cov.load_bam(path, panel["Exons"], amplicons=False)
            
            clinicalreport.create(cov, panel, os.path.join(clinical_report_path, output_stem))
            covermiplot.plot(cov, panel, os.path.join(clinical_report_path, output_stem))
            seconds = int(time.time() - start_time)
            time_string = "{0} sec".format(seconds) if (seconds<60) else "{0} min {01} sec".format(seconds/60, seconds%60)
            print"file {0} of {1} completed in {2}".format(len(output_stems), len(bam_file_list), time_string)
    else:
        cov = Cov.perfect_coverage(panel["Amplicons"])
        designreport.create(cov, panel, os.path.join(output_path, panel["Filenames"]["Panel"]))
        covermiplot.plot(cov, panel, os.path.join(output_path, panel["Filenames"]["Panel"]))


if __name__ == "__main__":

    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:b:o:", ["panel=", "bams=", "output="])
    except getopt.GetoptError as err:
        raise CoverMiException(str(err))

    output = None
    bams = None
    panel = None
    depth = None
    for o, a in opts:
        if o in ("-p", "--panel"):
            panel = a
        elif o in ("-b", "--bams"):
            bams = a
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-d", "--depth"):
            depth = a
        else:
            raise CoverMiException("Unrecognised option {0}".format(o))

    main(panel, bams, output, depth)









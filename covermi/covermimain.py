#!/usr/bin env python
import sys, os, time, re, getopt, pdb
import covermipanel, coverage, technicalreport, clinicalreport, designreport, covermiplot


class CoverMiException(Exception):
    pass


def identify_bam_files(bam_path):
    bam_file_list = []
    if os.path.isfile(bam_path) and bam_path.endswith(".bam"):
        bam_file_list.append(bam_path)
    elif not os.path.isdir(bam_path):
        raise CoverMiException("'{0}' is not a bam file or a directory".format(bam_path))
    else:
        for root, dirnames, filenames in os.walk(bam_path):
            for filename in filenames:
                if filename.endswith(".bam"):
                    bam_file_list.append(os.path.join(root, filename))
    if len(bam_file_list) == 0:
        raise CoverMiException("No bam files found at '{0}'".format(bam_path))
    return bam_file_list


def create_output_dir(output_path):
    output_path = os.path.join(output_path, "covermi_output")
    try:
        os.mkdir(output_path)
    except OSError:
        if os.path.isdir(output_path):
            raise CoverMiException("{0} folder already exists".format(output_path))
        else:
            raise CoverMiException("Unable to create folder {0}".format(output_path))
    return output_path


def covermi_main(panel_path, bam_path, output_path, depth=None):

    panel = covermipanel.load(panel_path, bam_path=="")
    if depth is not None:
        panel["Depth"] = int(depth)
    output_path = create_output_dir(output_path)
    print "Processing..."

    if bam_path != "":
        bam_file_list = identify_bam_files(bam_path)
        if len(bam_file_list) == 1:
            clinical_report_path = output_path
            technical_report_path = output_path
        else:
            clinical_report_path = os.path.join(output_path, "clinical")
            technical_report_path = os.path.join(output_path, "technical")
            os.mkdir(clinical_report_path)
            os.mkdir(technical_report_path)

        output_stems = set([])
        for bam_file in bam_file_list:
            start_time = time.time()
            print bam_file
            panel["Filenames"]["Sample"] = os.path.basename(bam_file).split("_")[0]
            panel["Filenames"]["Run"] = os.path.basename(os.path.dirname(bam_file))

            output_stem = panel["Filenames"]["Sample"]
            dup_num = 1
            while output_stem in output_stems:
                output_stem = "{0}({1})".format(panel["Filenames"]["Sample"], dup_num)
                dup_num += 1
            output_stems.add(output_stem)

            if "Amplicons" in panel:
                cov = coverage.load_bam_amplicons(bam_file, panel["Amplicons"])
                technicalreport.create(cov.amplicon_data, panel, os.path.join(technical_report_path, output_stem))
            else:
                cov = coverage.load_bam_exons(bam_file, panel["Exons"])

            clinicalreport.create(cov, panel, os.path.join(clinical_report_path, output_stem))
            covermiplot.plot(cov, panel, os.path.join(clinical_report_path, output_stem))
            seconds = int(time.time() - start_time)
            time_string = "{0} sec".format(seconds) if (seconds<60) else "{0} min {01} sec".format(seconds/60, seconds%60)
            print"file {0} of {1} completed in {2}".format(len(output_stems), len(bam_file_list), time_string)

    else:
        cov = coverage.perfect_coverage(panel["Amplicons"])
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

    covermimain(panel, bams, output, depth)









#!/usr/bin env python
import sys, os, time, re, getopt, pdb
from cov import Cov
from panel import Panel
from files import Files
import clinicalreport, covermiplot

class CoverMiException(Exception):
    pass


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


def main(panel_path, bamlist, output_path):
    panel = Panel(panel_path).load()
    output_path = create_output_dir(output_path)
    print "Processing..."
    if len(set([os.path.splitext(os.path.basename(fn))[0].split("_")[0] for fn, depths in bamlist])) > len(bamlist):
        print "Duplicate sample numbers. Quiting"
        exit()

    starttime = time.time()
    for bam_path, depths in bamlist:
        panel["Options"]["Depths"] = depths
        panel["Filenames"]["Sample"] = os.path.splitext(os.path.basename(bam_path))[0].split("_")[0]
        cov = Cov.load_bam(bam_path, panel["Exons"], amplicons=False)
        for depth in depths:
            panel["Options"]["Depth"] = depth
            clinicalreport.create(cov, panel, os.path.join(output_path, "{0}_{1}x_".format(panel["Filenames"]["Sample"], depth)))
            seconds = int(time.time() - starttime)
            time_string = "{0} sec".format(seconds) if (seconds<60) else "{0} min {01} sec".format(seconds/60, seconds%60)
            print "File {}, Depth {}, {}".format(panel["Filenames"]["Sample"], depth, time_string)
        panel["Options"]["Depth"] = max(depths)
        covermiplot.plot(cov, panel, os.path.join(output_path, panel["Filenames"]["Sample"]))


#    try:
#        opts, args = getopt.getopt(sys.argv[1:], "p:b:o:", ["panel=", "bams=", "output="])
#    except getopt.GetoptError as err:
#        raise CoverMiException(str(err))

#    output = None
#    bams = None
#    panel = None
#    depth = None
#    for o, a in opts:
#        if o in ("-p", "--panel"):
#            panel = a
#        elif o in ("-b", "--bams"):
#            bams = a
#        elif o in ("-o", "--output"):
#            output = a
#        elif o in ("-d", "--depth"):
#            depth = a
#        else:
#            raise CoverMiException("Unrecognised option {0}".format(o))

#    main(panel, bams, output, depth)









#!/usr/bin env python
import sys, os, tkFileDialog, Tkinter, pdb, traceback, pkg_resources

from cov import Cov
from panel import Panel
from files import Files
from itertools import izip
from covermimain import CoverMiException
import covermiconf, technicalreport, clinicalreport


def main():
    try:
        version = pkg_resources.require("CoverMi")[0].version

        rootwindow = Tkinter.Tk()
        rootwindow.withdraw()

        print "CoverMi version {}".format(version)
        print("Please select the verification directory. This must contain a single panel and one or more bams run with that panel.")
        path = tkFileDialog.askdirectory(parent=rootwindow, title='Please select the verification directory')
        if not bool(path):
            sys.exit()
        path = os.path.abspath(path)
        print("{0} directory selected".format(os.path.basename(path)))

        for root, dirnames, filenames in os.walk(path):
            break
        
        panel = None
        for dirname in dirnames:
            try:
                panel = Panel(os.path.join(root, dirname)).load()
                break
            except CoverMiException:
                pass
        if panel is None:
            raise CoverMiException("No panel within verification directory.")
        if "Amplicons" not in panel:
            raise CoverMiException("No amplicon information within panel.")

        bams = Files(path, ".bam")
        suffixes = ("_amp_covermi_clinical_report.txt", "_noamp_covermi_clinical_report.txt", "_covermi_technical_report.txt", "_amplicon_data.tsv")

        matched = 0
        for run, sample, path in bams:
            matched += sum([os.path.exists(os.path.join(root, sample+suffix)) for suffix in suffixes])

        if matched == 0:
            verify = ""
            print "Creating baseline verification files for CoverMi version {}".format(version)
        elif matched == len(suffixes) * len(bams):
            verify = "_"
            print "Verifying CoverMi version {}".format(version)
        else:
            raise CoverMiException("Incomplete number of verification files.")

        for run, sample, path in bams:
            panel["Filenames"]["Run"] = "verification"
            panel["Filenames"]["Sample"] = sample
            print "Processing test bam {}.".format(sample)
            path += ".bam"
            cov = Cov.load_bam(path, panel["Amplicons"], amplicons=True)
            technicalreport.create(cov.amplicon_info, panel, os.path.join(root, verify+sample))
            clinicalreport.create(cov, panel, os.path.join(root, verify+sample+"_amp"))

            cov = Cov.load_bam(path, panel["Exons"], amplicons=False)
            clinicalreport.create(cov, panel, os.path.join(root, verify+sample+"_noamp"))

            if verify:
                for fn in [sample+suffix for suffix in suffixes]:
                    with file(os.path.join(root, fn), "rt") as old:
                        with file(os.path.join(root, verify+fn), "rt") as new:
                            for new_row, old_row in izip(new, old):
                                if new_row != old_row:
                                    if not (new_row.startswith("Date of CoverMi analysis:") or new_row.startswith("CoverMi version:")):
                                        raise CoverMiException("VERIFICATIO FAILED in file {}\nThe following lines do not match.\n{}{}".format(fn, new_row, old_row))
                    print "{} verified.".format(fn)
                    os.unlink(os.path.join(root, verify+fn))



    except Exception as e:
        if type(e).__name__ == "CoverMiException":
            print(e.message)
        else:
            traceback.print_exc()
            print("UNEXPECTED ERROR. QUITTING.")
    finally:
        raw_input("Press any key to continue...")

if __name__ == "__main__":
    main()

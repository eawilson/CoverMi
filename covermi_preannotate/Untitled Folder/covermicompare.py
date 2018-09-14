#!/usr/bin env python
#!/usr/bin env python
import sys, os, tkFileDialog, Tkinter, pdb, traceback
import covermiconf
from panel import Panel


class DepthDialog(object):
    def __init__(self, parent, default_depth):
        self.parent = parent
        self.parent.depth = ""
        self.window = Tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        Tkinter.Label(self.window, text="Please enter minimum depth").grid(column=0, row=0)
        self.entry = Tkinter.Entry(self.window)
        self.entry.grid(column=1, row=0)
        self.entry.insert(0, str(default_depth))
        self.entry.bind('<Return>', self.return_pressed)
        self.entry.focus_set()

    def return_pressed(self, event):
        if self.entry.get().isdigit():
            self.parent.depth = self.entry.get()
            self.window.destroy()


class M_S_D_Dialog(object):
    def __init__(self, parent):
        self.parent = parent
        self.parent.msd = ""
        self.window = Tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        Tkinter.Label(self.window, text="Do you want to check a single bam, multiple bams or review the panel design?").grid(column=0, row=0, columnspan=3, padx=10, pady=5)
        Tkinter.Button(self.window, text="Folder of bam files", width=15, command=self.multiple_pressed).grid(column=0, row=1, pady=10)
        Tkinter.Button(self.window, text="Single bam file", width=15, command=self.single_pressed).grid(column=1, row=1, pady=10)
        Tkinter.Button(self.window, text="Review panel design", width=15, command=self.design_pressed).grid(column=2, row=1, pady=10)

    def single_pressed(self):
        self.parent.msd = "single"
        self.window.destroy()

    def multiple_pressed(self):
        self.parent.msd = "multiple"
        self.window.destroy()

    def design_pressed(self):
        self.parent.msd = "design"
        self.window.destroy()


class SomCon_Dialog(object):
    def __init__(self, parent):
        self.parent = parent
        self.parent.somcon = ""
        self.window = Tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        Tkinter.Label(self.window, text="Is this a somatic or constitutional panel?").grid(column=0, row=0, columnspan=2, padx=10, pady=5)
        Tkinter.Button(self.window, text="Somatic", width=15, command=self.som_pressed).grid(column=0, row=1, pady=10, padx=10)
        Tkinter.Button(self.window, text="Constitutional", width=15, command=self.con_pressed).grid(column=1, row=1, pady=10, padx=10)

    def som_pressed(self):
        self.parent.somcon = "Somatic"
        self.window.destroy()

    def con_pressed(self):
        self.parent.somcon = "Constitutional"
        self.window.destroy()


def guimain():
    try:
        conf = covermiconf.load_conf()

        rootwindow = Tkinter.Tk()
        rootwindow.withdraw()

        print("Please select a panel")
        panelpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=conf["panel_root"], title='Please select a panel')
        if not bool(panelpath):
            sys.exit()
        panelpath = os.path.abspath(panelpath)
        print("{0} panel selected".format(os.path.basename(panelpath)))

        panel = Panel(panelpath)
        options = panel.read_options()
        rootwindow.wait_window(DepthDialog(rootwindow, str(options["Depth"]) if ("Depth" in options) else "").window)
        if rootwindow.depth == "":
            sys.exit()
        options["Depth"] = int(rootwindow.depth)
        print("Depth {0} selected".format(options["Depth"]))

        if "ReportType" not in options:
            print("Is this a somatic or constitutional panel?")
            rootwindow.wait_window(SomCon_Dialog(rootwindow).window)    
            if rootwindow.somcon == "":
                sys.exit()
            options["ReportType"] = rootwindow.somcon

        panel.write_options(options)

        print("Do you wish to coverage check multiple bams, a single bam or review the panel design?")
        rootwindow.wait_window(M_S_D_Dialog(rootwindow).window)    
        mode = str(rootwindow.msd)
        if mode == "":
            sys.exit()
        elif mode == "design":
            bampath = ""
            outputpath = conf["panel_root"]  
            print("Design review selected")
        else:
            if mode == "multiple":
                print("Please select the folder containing the bam files")
                bampath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=conf["bam_root"], title='Please select a folder')
            elif mode == "single":
                print("Please select a bam file")
                bampath = tkFileDialog.askopenfilename(parent=rootwindow, initialdir=conf["bam_root"], filetypes=[("bamfile", "*.bam")], title='Please select a bam file')
            if bampath == "":
                sys.exit()
            bampath = os.path.abspath(bampath)
            outputpath = os.path.dirname(bampath) if (mode=="single") else bampath

            print("{0} selected".format(bampath))

        print("Please select a location for the output")
        outputpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=outputpath, title='Please select a location for the output')
        if outputpath == "":
            sys.exit()
        outputpath = os.path.abspath(outputpath)
        print("Output folder {0} selected".format(outputpath))

        main(panelpath, bampath, outputpath, "old")
        main(panelpath, bampath, outputpath, "python")

        print("Finished")
    except Exception as e:
        if type(e).__name__ == "CoverMiException":
            print(e.message)
        else:
            traceback.print_exc()
            print("UNEXPECTED ERROR. QUITTING.")
    finally:
        raw_input("Press any key to continue...")


import sys, os, time, re, getopt, pdb
from cov import Cov
from panel import Panel
from files import Files
import covermiconf, technicalreport, clinicalreport, designreport, covermiplot
#import testreport as clinicalreport
from covermiexception import CoverMiException


def create_output_dir(output_path, bam_path, suffix):
    output_path = os.path.join(output_path, "{0}_covermi_output_{1}".format(os.path.splitext(os.path.basename(bam_path))[0], suffix))
    try:
        os.mkdir(output_path)
    except OSError:
        if os.path.isdir(output_path):
            raise CoverMiException("{0} folder already exists".format(output_path))
        else:
            raise CoverMiException("Unable to create folder {0}".format(output_path))
    return output_path


def main(panel_path, bam_path, output_path, test, depth=None):

    if not os.path.exists(panel_path):
        panel_path = os.path.join(covermiconf.load_conf["panel_root"], panel_path)

    panel = Panel(panel_path).load(bam_path=="")
    if depth is not None:
        panel["Options"]["Depth"] = int(depth)
    output_path = create_output_dir(output_path, bam_path if bam_path!="" else panel_path, test)
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

            if test == "python":
                cov = Cov.load_bam_python(path, panel["Amplicons"])            
            elif test == "old":
                if "Amplicons" in panel:
                    cov = Cov.load_bam(path, panel["Amplicons"], amplicons=True)
                    technicalreport.create(cov.amplicon_info, panel, os.path.join(technical_report_path, output_stem))
                else:
                    cov = Cov.load_bam(path, panel["Exons"], amplicons=False)
            else:
                raise CoverMiException("Hard coding error!")
            
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

    main(panel, bams, output, "python", depth)
    main(panel, bams, output, "old", depth)








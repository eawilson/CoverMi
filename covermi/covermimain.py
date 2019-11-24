import sys
import os
import argparse
import platform
try:
    import tkinter
    from tkinter import filedialog
except ImportError:
    pass
import time
import re
import shutil
from collections import namedtuple
import pdb

from .panel import Panel
from .gr import variants
from .include import CoverMiException, __version__
from . import technicalreport, clinicalreport, designreport, covermiplot, tablereport
from .cov import Cov, fake_paired_end_reads



class DepthDialog(object):
    def __init__(self, parent, default_depth):
        self.parent = parent
        self.parent.depth = ""
        self.window = tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        tkinter.Label(self.window, text="Please enter minimum depth").grid(column=0, row=0)
        self.entry = tkinter.Entry(self.window)
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
        self.window = tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        tkinter.Label(self.window, text="Do you want to check a single bam, multiple bams or review the panel design?").grid(column=0, row=0, columnspan=3, padx=10, pady=5)
        tkinter.Button(self.window, text="Folder of bam files", width=15, command=self.multiple_pressed).grid(column=0, row=1, pady=10)
        tkinter.Button(self.window, text="Single bam file", width=15, command=self.single_pressed).grid(column=1, row=1, pady=10)
        tkinter.Button(self.window, text="Review panel design", width=15, command=self.design_pressed).grid(column=2, row=1, pady=10)

    def single_pressed(self):
        self.parent.msd = "single"
        self.window.destroy()

    def multiple_pressed(self):
        self.parent.msd = "multiple"
        self.window.destroy()

    def design_pressed(self):
        self.parent.msd = "design"
        self.window.destroy()



def main():
    
    try:
        if len(sys.argv) > 1: # Command line arguments
            parser = argparse.ArgumentParser()
            parser.add_argument("-p", "--panel", help="Path of panel directory.", required=True)
            parser.add_argument("-i", "--input", help="Path of input bam or run directory.")
            parser.add_argument("-o", "--output", help="Path of output directory.")
            parser.add_argument("-d", "--depth", type=int, help="Panel depth (optional).")
            parser.add_argument("--profile", help="File in which to save profiling information (development use only).")
            parser.add_argument("--overwrite", help="Overwrite previous output.", action="store_true")
            args = parser.parse_args()
            
            if args.input is None and args.output is None: # Initialise panel
                panel = Panel(args.panel)
                print("Panel = {}.".format(panel.name))
                
                if "properties" in panel.files and not args.overwrite:
                    print("Properties file {} already exists.".format(panel.files["properties"]))
                else:
                    fn = "covermi_{}.txt".format(panel.name)
                    path = os.path.join(args.panel, fn)
                    if os.path.exists(path):
                        raise RuntimeError("Properties file {} already exists but is of incorrect format.".format(fn))
                    print("Writing properties file {}.".format(fn))
                    
                    with open(path, "wt") as f:
                        f.write("#covermi\n")
                        f.write("#assembly=GRCH37\n")
                        f.write("#transcript_source=refseq/ensembl\n")
                        f.write("#reporttype=somatic/constitutional\n")
                        f.write("#depth=100\n")
                        f.write("#filters=filters == 'PASS'\n")
                                    
                if "variants" in panel.files and not args.overwrite:
                    if "diseases" in panel.files:
                        print("Diseases file {} already exists.".format(panel.files["properties"]))
                    else:
                        fn = "diseass_{}.txt".format(panel.name)
                        path = os.path.join(args.panel, fn)
                        if os.path.exists(path):
                            raise RuntimeError("Diseases file {} already exists but is of incorrect format.".format(fn))
                        print("Writing diseases file {}.".format(fn))
                        
                        with open(path, "wt") as f:
                            f.write("#diseases\n")
                            f.write("\n".join(sorted(set(variant.name for variant in variants(panel.files["variants"], "disease")))))
                            f.write("\n")
            
            elif args.output is None:
                parser.error("the following arguments are required: -o/--output")

            else:
                if args.profile:
                    from cProfile import Profile
                    from pstats import Stats
                    profiler = Profile()
                    profiler.runctx('covermimain(args.panel, args.output, bam_path=args.input, depth=args.depth, overwrite=args.overwrite)', globals(), locals())
                    with open(args.profile + ".txt", "wt") as f:
                        stats = Stats(profiler, stream=f).strip_dirs().sort_stats("cumulative")
                        stats.print_stats()
                    stats.dump_stats(args.profile + ".stats")

                else:
                    covermimain(args.panel, args.output, bam_path=args.input, depth=args.depth, overwrite=args.overwrite)
                
        else: # No arguments therefore use gui
            rootwindow = tkinter.Tk()
            rootwindow.withdraw()

            print("Please select a panel.")
            panelpath = filedialog.askdirectory(parent=rootwindow, initialdir="~", title='Please select a panel')
            if not bool(panelpath):
                sys.exit()
            panelpath = os.path.abspath(panelpath)
            print("{0} panel selected.".format(os.path.basename(panelpath)))

            print("Please select minimum coverage depth.")
            rootwindow.wait_window(DepthDialog(rootwindow, Panel(panelpath).properties.get("depth", "")).window)
            if rootwindow.depth == "":
                sys.exit()
            depth = rootwindow.depth
            print("Depth {} selected.".format(depth))

            print("Do you wish to coverage check multiple bams, a single bam or review the panel design?")
            rootwindow.wait_window(M_S_D_Dialog(rootwindow).window)    
            mode = str(rootwindow.msd)
            if mode == "":
                sys.exit()
            elif mode == "design":
                bampath = ""
                print("Design review selected.")
            else:
                if mode == "multiple":
                    print("Please select the folder containing the bam files.")
                    bampath = filedialog.askdirectory(parent=rootwindow, initialdir="~", title='Please select a folder')
                elif mode == "single":
                    print("Please select a bam file.")
                    bampath = filedialog.askopenfilename(parent=rootwindow, initialdir="~", filetypes=[("bamfile", "*.bam")], title='Please select a bam file')
                if bampath == ():
                    sys.exit()
                bampath = os.path.abspath(bampath)
                print("{} selected.".format(bampath))

            print("Please select a location for the output.")
            outputpath = filedialog.askdirectory(parent=rootwindow, initialdir=bampath, title='Please select a location for the output')
            if outputpath == ():
                sys.exit()
            outputpath = os.path.abspath(outputpath)
            print("Output location {0} selected.".format(outputpath))

            print('covermi -p "{}" -o "{}" -i "{}" -d {}'.format(panelpath, outputpath, bampath, depth))
            covermimain(panelpath, outputpath, bam_path=bampath, depth=depth)
            
    except CoverMiException as e:
        print(e)
        
    finally:
        # In windows if started from the gui the terminal will close too quickly to allow for warnings and errors to be seen
        if platform.system == "Windows" and len(sys.argv) == 1: 
            print("Finished.")
            input("Press enter to continue...")



def walk_samples(path, ext=".bam", strip_trailing_underscore=True):
    
    Sample = namedtuple("Sample", "name run path")

    illumina_suffix = re.compile("(.+)_S[0-9]+$")
    def strip_underscore(sample):
        matchobj = illumina_suffix.search(sample)
        return matchobj.group(1) if matchobj else sample 

    def walk_run(path, ext):
        run = os.path.basename(path)
        for root, dirnames, filenames in os.walk(path):
            for fn in filenames:
                sample, extension = os.path.splitext(fn)
                if os.path.isfile(os.path.join(root, fn)) and extension == ext:
                    yield Sample(strip_underscore(sample), run, os.path.join(root, fn))
            
    found = []
    if os.path.isfile(path):
        if not os.path.splitext(path)[1] == ext:
            raise CoverMiException("File "+path+" is not of type {}.".format(ext))
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
        raise CoverMiException("Invalid bam path {}.".format(path))



def covermimain(panel_path, output_path="", bam_path=None, depth=None, overwrite=False):
    print("CoverMi v{} (Python {}.{}.{})".format(__version__, *sys.version_info[:3]))
    print("Processing...")
    
    output_path = os.path.join(output_path, "{0}_covermi_output".format(os.path.splitext(os.path.basename(bam_path if bam_path!="" else panel_path))[0]))
    if overwrite and os.path.isdir(output_path):
        for root, dns, fns in os.walk(output_path):
            for fn in fns:
                if os.path.splitext(fn)[1] not in (".pdf", ".txt", ".tsv"):
                    raise CoverMiException("Unable to overwrite {0} as contains unexpected contents.".format(output_path))
        try:
            shutil.rmtree(output_path)
        except Exception:
            raise CoverMiException("Unable to overwrite {0}.".format(output_path))
        
    try:
        os.mkdir(output_path)
    except Exception:
        if os.path.isdir(output_path):
            raise CoverMiException("{0} already exists.".format(output_path))
        else:
            raise CoverMiException("Unable to create {0}.".format(output_path))


    panel = Panel(panel_path, verbose=True, splice_site_buffer=5)
    if depth is not None:
        panel.properties["depth"] = depth

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

            if cov.amplicon_info:
                technicalreport.create(cov.amplicon_info, panel, sample, os.path.join(technical_report_path, output_stem))
            clinicalreport.create(cov, panel, sample, os.path.join(clinical_report_path, output_stem))
            covermiplot.plot(cov, panel, sample.name, os.path.join(clinical_report_path, output_stem+"_plot.pdf"))

            seconds = int(time.time() - start_time)
            time_string = "{0} sec".format(seconds) if (seconds<60) else "{0} min {01} sec".format(seconds//60, seconds%60)
            print("File {0} of {1} completed in {2}.".format(len(output_stems), len(bam_file_list), time_string))

    else:
        cov = Cov(fake_paired_end_reads(panel.amplicons, depth=panel.depth*2))
        designreport.create(cov, panel, None, os.path.join(output_path, panel.name))
        covermiplot.plot(cov, panel, panel.name, os.path.join(output_path, panel.name+"_plot.pdf"))

    return output_path
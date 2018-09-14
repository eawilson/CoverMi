#!/usr/bin env python
from __future__ import print_function, absolute_import, division

import sys, os, pdb, traceback, argparse

try: # python2
    import Tkinter as tkinter
    import tkFileDialog as filedialog
except ImportError: #python3
    import tkinter
    from tkinter import filedialog
    
from .covermimain import covermimain
from .panel import Panel
from .gr import variants
from .include import CoverMiException

try: # python2
    input = raw_input
except NameError:
    pass


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


class SomCon_Dialog(object):
    def __init__(self, parent):
        self.parent = parent
        self.parent.somcon = ""
        self.window = tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        tkinter.Label(self.window, text="Is this a somatic or constitutional panel?").grid(column=0, row=0, columnspan=2, padx=10, pady=5)
        tkinter.Button(self.window, text="Somatic", width=15, command=self.som_pressed).grid(column=0, row=1, pady=10, padx=10)
        tkinter.Button(self.window, text="Constitutional", width=15, command=self.con_pressed).grid(column=1, row=1, pady=10, padx=10)

    def som_pressed(self):
        self.parent.somcon = "somatic"
        self.window.destroy()

    def con_pressed(self):
        self.parent.somcon = "constitutional"
        self.window.destroy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path of input bam or run directory")
    parser.add_argument("-o", "--output", help="path of output directory")
    parser.add_argument("-p", "--panel", help="path of panel directory")
    parser.add_argument("-d", "--depth", type=int, help="panel depth (optional)")
    args = parser.parse_args()
    if None not in (args.input, args.output, args.panel):
        return covermimain(args.panel, args.input, args.output, depth=args.depth)

    try:
        rootwindow = tkinter.Tk()
        rootwindow.withdraw()

        print("Please select a panel")
        panelpath = filedialog.askdirectory(parent=rootwindow, initialdir="~", title='Please select a panel')
        if not bool(panelpath):
            sys.exit()
        panelpath = os.path.abspath(panelpath)
        print("{0} panel selected".format(os.path.basename(panelpath)))
        panel = Panel(panelpath)

        if "variants" in panel.files and "diseases" not in panel.files:
            path = os.path.join(panel.path, panel.name.lower().replace(" ", "_")+"_diseases.txt")
            if os.path.exists(path):
                raise CoverMiException("Incorrect format of diseases file {}".format(os.path.basename(path)))
            try:
                with open(path, "wt") as f:
                    f.write("#diseases\n")
                    for disease in sorted(set(variant.name for variant in variants(panel.files["variants"], "disease"))):
                        f.write(disease+"\n")
            except IOError:
                pass

        update = {}
        rootwindow.wait_window(DepthDialog(rootwindow, panel.properties.get("depth", "")).window)
        if rootwindow.depth == "":
            sys.exit()
        depth = rootwindow.depth
        print("Depth {} selected".format(depth))

        if "depth" not in panel.properties:
            update["depth"] = depth

        if "reporttype" not in panel.properties:
            print("Is this a somatic or constitutional panel?")
            rootwindow.wait_window(SomCon_Dialog(rootwindow).window)    
            if rootwindow.somcon == "":
                sys.exit()
            update["reporttype"] = rootwindow.somcon

        if update:
            if "properties" in panel.files:
                path = panel.files["properties"]
            else:
                path = os.path.join(panel.path, panel.name+"_properties.txt")
                if os.path.exists(path):
                    raise CoverMiException("Incorrect format of properties file {}".format(os.path.basename(path)))
            try:
                with open(path, "at") as f:
                    for key, val in update.items():
                        f.write("{}={}\n".format(key, val))
            except IOError:
                pass

        print("Do you wish to coverage check multiple bams, a single bam or review the panel design?")
        rootwindow.wait_window(M_S_D_Dialog(rootwindow).window)    
        mode = str(rootwindow.msd)
        if mode == "":
            sys.exit()
        elif mode == "design":
            bampath = ""
            print("Design review selected")
        else:
            if mode == "multiple":
                print("Please select the folder containing the bam files")
                bampath = filedialog.askdirectory(parent=rootwindow, initialdir="~", title='Please select a folder')
            elif mode == "single":
                print("Please select a bam file")
                bampath = filedialog.askopenfilename(parent=rootwindow, initialdir="~", filetypes=[("bamfile", "*.bam")], title='Please select a bam file')
            if bampath == ():
                sys.exit()
            bampath = os.path.abspath(bampath)
            print("{} selected".format(bampath))

        print("Please select a location for the output")
        outputpath = filedialog.askdirectory(parent=rootwindow, initialdir=bampath, title='Please select a location for the output')
        if outputpath == ():
            sys.exit()
        outputpath = os.path.abspath(outputpath)
        print("Output folder {0} selected".format(outputpath))

        covermimain(panelpath, bampath, outputpath, depth=depth)

        print("Finished")
    except Exception as e:
        if type(e).__name__ == "CoverMiException":
            print(e.message)
        else:
            traceback.print_exc()
            print("UNEXPECTED ERROR. QUITTING.")
    finally:
        input("Press enter to continue...")

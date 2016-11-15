#!/usr/bin env python
import sys, os, tkFileDialog, Tkinter, pdb, traceback
import covermimain, covermiconf
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


def main():
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

        covermimain.main(panelpath, bampath, outputpath)

        print("Finished")
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

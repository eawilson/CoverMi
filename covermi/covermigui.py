#!/usr/bin env python
import sys, os, tkFileDialog, Tkinter, pdb, traceback
import covermimain

try:
    input = raw_input
except NameError:
    pass

class CoverMiException(Exception):
    pass
        

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


if __name__ == "__main__":
    try:
        
        root_dir = os.path.dirname(os.path.abspath(__file__))
        root_root_dir = os.path.dirname(root_dir)
        if os.path.isdir(os.path.join(root_dir, "panels")):
            root_dir = os.path.join(root_dir, "panels")
        elif os.path.isdir(os.path.join(root_root_dir, "panels")):
            root_dir = os.path.join(root_root_dir, "panels")

        rootwindow = Tkinter.Tk()
        rootwindow.withdraw()

        print("Please select a panel")
        panelpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=root_dir, title='Please select a panel')
        if not bool(panelpath):
            sys.exit()
        print("{0} panel selected".format(os.path.basename(panelpath)))

        path = covermimain.covermipanel.identify(panelpath)
        if "Depth" in path:
            depthpath = path["Depth"]
            with file(depthpath, "rU") as f:
                default_depth = f.readline().strip()
        else:
            depthpath = os.path.join(panelpath, "depth")
            if os.path.lexists(depthpath):
                raise CoverMiException("File '{0}' exists but is of incorrect format".format(depthpath))
            default_depth = ""
        rootwindow.wait_window(DepthDialog(rootwindow, default_depth).window)
        depth = rootwindow.depth
        if depth == "":
            sys.exit()
        if depth != default_depth:
            with file(depthpath, "wt") as f:
                f.write("{0}\n".format(depth))
        print("Depth {0} selected".format(depth))
        depth = int(depth)

        print("Do you wish to coverage check multiple bams, a single bam or review the panel design?")
        rootwindow.wait_window(M_S_D_Dialog(rootwindow).window)    
        mode = str(rootwindow.msd)
        if mode == "":
            sys.exit()
        elif mode == "design":
            bampath = ""
            outputpath = root_dir    
            print("Design review selected")
        else:
            if mode == "multiple":
                print("Please select the folder containing the bam files")
                bampath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=root_dir, title='Please select a folder')
            elif mode == "single":
                print("Please select a bam file")
                bampath = tkFileDialog.askopenfilename(parent=rootwindow, initialdir=root_dir, filetypes=[("bamfile", "*.bam")], title='Please select a bam file')
            if bampath == "":
                sys.exit()
            outputpath = os.path.dirname(bampath) if (mode=="single") else bampath

            print("{0} selected".format(bampath))

        print("Please select a folder for the output")   
        outputpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=outputpath, title='Please select a folder for the output')
        if outputpath == "":
            sys.exit()
        print("Output folder {0} selected".format(outputpath))

        covermimain.covermi_main(panelpath, bampath, outputpath)

        print("Finished")
    except Exception as e:
        if type(e).__name__ == "CoverMiException":
            print(e.message)
        else:
            traceback.print_exc()
            print("UNEXPECTED ERROR. QUITTING.")
    finally:
        raw_input("Press any key to continue...")


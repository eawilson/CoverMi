#!/usr/bin env python
import sys, os, tkFileDialog, Tkinter, pdb, traceback
import covermiwgsmain
from panel import Panel


class DepthDialog(object):
    def __init__(self, parent, default_depth):
        self.parent = parent
        self.parent.depth = ""
        self.window = Tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        Tkinter.Label(self.window, text="Please enter a comma separated list of depths").grid(column=0, row=0)
        self.entry = Tkinter.Entry(self.window)
        self.entry.grid(column=1, row=0)
        self.entry.insert(0, str(default_depth))
        self.entry.bind('<Return>', self.return_pressed)
        self.entry.focus_set()

    def return_pressed(self, event):
        depths = [depth.strip() for depth in self.entry.get().split(",")]
        if all([depth.isdigit() for depth in depths]):
            self.parent.depths = [int(depth) for depth in depths]
            self.window.destroy()


class YesNo_Dialog(object):
    def __init__(self, parent):
        self.parent = parent
        self.parent.yesno = ""
        self.window = Tkinter.Toplevel(self.parent)
        self.window.title("CoverMi")
        Tkinter.Label(self.window, text="Do you wish to coverage check any more bams?").grid(column=0, row=0, columnspan=2, padx=10, pady=5)
        Tkinter.Button(self.window, text="Yes", width=15, command=self.yes_pressed).grid(column=0, row=1, pady=10, padx=10)
        Tkinter.Button(self.window, text="No", width=15, command=self.no_pressed).grid(column=1, row=1, pady=10, padx=10)

    def yes_pressed(self):
        self.parent.yesno = "yes"
        self.window.destroy()

    def no_pressed(self):
        self.parent.yesno = "no"
        self.window.destroy()


def main():

    try:        
        root_dir = os.path.dirname(os.path.expanduser("~"))
        rootwindow = Tkinter.Tk()
        rootwindow.withdraw()

        print("Please select a panel")
        panelpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=root_dir, title='Please select a panel')
        if not bool(panelpath):
            sys.exit()
        print("{0} panel selected".format(os.path.basename(panelpath)))

        bamlist = []
        while True:
                print("Please select a bam file")
                bampath = tkFileDialog.askopenfilename(parent=rootwindow, initialdir=root_dir, filetypes=[("bamfile", "*.bam")], title='Please select a bam file')
                if bampath == "":
                    sys.exit()
                root_dir = os.path.dirname(bampath)
                print("{0} selected".format(bampath))

                rootwindow.wait_window(DepthDialog(rootwindow, "").window)
                if rootwindow.depths == "":
                    sys.exit()
                depths = rootwindow.depths
                print("Depth {0} selected".format(", ".join([str(depth) for depth in depths])))

                bamlist += [(bampath, depths)]

                rootwindow.wait_window(YesNo_Dialog(rootwindow).window)    
                if rootwindow.yesno == "":
                    sys.exit()
                elif rootwindow.yesno == "no":
                    break

        print("Please select a location for the output")   
        outputpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=root_dir, title='Please select a location for the output')
        if outputpath == "":
            sys.exit()
        print("Output location {0} selected".format(outputpath))

        covermiwgsmain.main(panelpath, bamlist, outputpath)

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

#!/usr/bin env python
import os, sys, tkFileDialog, Tkinter
from .gr import Gr, Entry
from .files import Files


def shared(proband=None, test=None, chromosome=None):
    if chromosome:
        proband = proband.touched_by(chromosome)
        test = test.touched_by(chromosome)
    p = set([entry.name for entry in proband.all_entries if entry.passfilter])
    t = set([entry.name for entry in test.all_entries if entry.passfilter])
    return "{:.1f}%".format(float(len(p & t)*100)/len(t)) if len(t)>0 else "NA"


def main():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    probandpath = tkFileDialog.askopenfilename(parent=rootwindow, initialdir=root_dir, title='Please select the proband vcf')
    if not bool(probandpath):
        sys.exit()
    print "proband: {}".format(os.path.basename(probandpath))

    testpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=root_dir, title='Please select the folder of the vcfs to be checked')
    if not bool(testpath):
        sys.exit()

    chrY = Gr(Entry("chrY", 1, Gr.MAX_CHR_LENGTH, "", "+"))
    chrM = Gr(Entry("chrM", 1, Gr.MAX_CHR_LENGTH, "", "+"))
    proband_variants = Gr.load_vcf(probandpath)

    for sample, run, testpath in Files(testpath, ".vcf"): 
        test_variants = Gr.load_vcf(testpath+".vcf")
        print "{} {}, common variants: {}, chrY: {}, chrM {}".format(
              os.path.basename(sample), 
              "male" if len(test_variants.touched_by(chrY))>0 else "female",
              shared(proband=proband_variants, test=test_variants), 
              shared(proband=proband_variants, test=test_variants, chromosome=chrY),
              shared(proband=proband_variants, test=test_variants, chromosome=chrM))

if __name__ == "__main__":
    main()

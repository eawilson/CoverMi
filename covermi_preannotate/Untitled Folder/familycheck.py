#!/usr/bin env python
import os, tkFileDialog, Tkinter
from .gr import Gr, Entry


def notshared(proband=None, test=None, chromosome=None):
    if chromosome:
        proband = proband.touched_by(chromosome)
        test = test.touched_by(chromosome)
    p = set([entry.name for entry in proband.all_entries if entry.passfilter])
    t = set([entry.name for entry in test.all_entries if entry.passfilter])
    return "{:.1f}%".format(float(len(p - t)*100)/len(p)) if len(p)>0 else "NA"


def main():
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    probandpath = tkFileDialog.askopenfilename(parent=rootwindow, title='Please select the probands vcf')
    if not bool(probandpath):
        sys.exit()
    print "proband: {}".format(os.path.basename(probandpath))

    motherpath = tkFileDialog.askopenfilename(parent=rootwindow, title='Please select the mothers vcf')
    if not bool(motherpath):
        sys.exit()
    print "mother: {}".format(os.path.basename(motherpath))

    fatherpath = tkFileDialog.askopenfilename(parent=rootwindow, title='Please select the fathers vcf')
    if not bool(fatherpath):
        sys.exit()
    print "father: {}".format(os.path.basename(fatherpath))

    chrY = Gr().add(Entry("chrY", 1, Gr.MAX_CHR_LENGTH, "", "+"))
    chrM = Gr().add(Entry("chrM", 1, Gr.MAX_CHR_LENGTH, "", "+"))
    proband_variants = Gr.load_vcf(probandpath)
    mother_variants = Gr.load_vcf(motherpath)
    father_variants = Gr.load_vcf(fatherpath)

    print "{} of the probands variants are not present in the parents".format(
        notshared(proband=proband_variants, test=mother_variants.combined_with(father_variants)))
    print "{} of the probands variants on chromosome M are not present in the mother".format(
        notshared(proband_variants, mother_variants, chrM))
    print "{} of the probands variants on chromosome Y are not present in the father".format(
        notshared(proband_variants, father_variants, chrY))

if __name__ == "__main__":
    main()

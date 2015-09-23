import sys, os, tkFileDialog, Tkinter, pdb


def main():
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print("Please select knownCanonical file downloaded from UCSC")
    kcpath = tkFileDialog.askopenfilename(parent=rootwindow, filetypes=[("textfile", "*.txt")], title='Please select knownCanonical file downloaded from UCSC')
    if not bool(kcpath):
        sys.exit()
    print("{0} file selected".format(os.path.basename(kcpath)))
    path = os.path.dirname(kcpath)
    print("Please select kgXref file downloaded from UCSC")
    kgpath = tkFileDialog.askopenfilename(parent=rootwindow, filetypes=[("textfile", "*.txt")], initialdir=path, 
                                          title='Please select kgXref file downloaded from UCSC')
    if not bool(kgpath):
        sys.exit()
    print("{0} file selected".format(os.path.basename(kgpath)))


    canonicals = set([])
    with file(kcpath) as f:
        for line in f:
            chrom, start, stop, cluster, kgid, protein = line.split("\t")
            canonicals.add(kgid)

    path = os.path.join(path, "covermi_"+os.path.basename(kcpath))
    if os.path.exists(path):
        raise RuntimeError("File {} already exists".format(path))
    with file(path, "wt") as g:
        with file(kgpath) as f:
            for line in f:
                kgid, mrna, upan, upid, gene, transcript, protein, des, rfam, rna = line.split("\t")
                if kgid in canonicals and transcript != "":
                    g.write("{0}\t{1}\n".format(gene, transcript))

    print "Written file {}".format(path)
    print "This file needs to be placed in all panel directories"


if __name__ == "__main__":
    main()

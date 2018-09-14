import sys, os, tkFileDialog, Tkinter, pdb


def main():
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print "BringMiUp"
    print "Please select a folder"
    folderpath = tkFileDialog.askdirectory(parent=rootwindow, title='BringMiUp. Please select a folder')
    if folderpath == "":
        sys.exit()

    skiplist = []
    print "{0} folder selected".format(os.path.basename(folderpath))
    for root, dirnames, filenames in os.walk(folderpath):
        if root != folderpath:
            for filename in filenames:
                if filename.endswith(".bam") or filename.endswith(".bai") or filename.endswith(".vcf"):

                    oldfile = os.path.join(root, filename)
                    newfile = os.path.join(folderpath, filename)
                    if os.path.exists(newfile):
                        skiplist.append(oldfile)
                        print "File \""+oldfile+"\" already exists in \""+folderpath+"\". Skipping"
                    else:
                        print "Moving "+filename
                        os.rename(oldfile, newfile)

    if len(skiplist) > 0:
        print "WARNING. The following files were not moved"
        print ", ".join(skiplist)


if __name__ == "__main__":
    main()

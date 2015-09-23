#!/usr/bin env python
from panel import Panel
import sys, os, tkFileDialog, Tkinter, re
import pdb




def main():
    chr_translate = { "23":"X", "24":"Y", "25":"M" }

    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print "Please select the COSMIC mutation tsv file"
    cosmic_path = tkFileDialog.askopenfilename(parent=rootwindow, title='Please select the COSMIC mutation file')
    if cosmic_path == "":
        sys.exit()
    if  len(file(cosmic_path, "rU").readline().split("\t")) != 29:
        print "File {0} is of incorrect format".format(os.path.basename(cosmic_path))
        sys.exit()

    print "Please select a panel"
    panel_path = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a panel')
    if panel_path == "":
        sys.exit()

    panel = Panel(panel_path)
    if "Variants" in panel:
        print "Variants file {0} already exists within panel".format(os.path.basename(panel["Variants"]))
        sys.exit()
    if "Targets" not in panel:
        print "Unable to identify target gene file. Exiting"
        sys.exit()
    if "Disease_Names" in panel:
        print "Disease names file {0} will be deleted and updated".format(os.path.basename(panel["Disease_Names"]))
        if raw_input("Press y to continue, any other key to abort") != "y":
            sys.exit()

    genes = set([])
    with file(panel["Targets"] ,"rU") as f:
        for line in f:
            gene = line.strip().split()[0]
            if gene != "":
                genes.add(gene)

    outputfilepath = os.path.join(panel_path, 
                              "{0}_{1}.txt".format(os.path.basename(panel_path.rstrip(os.pathsep)), os.path.splitext(os.path.basename(cosmic_path.rstrip(os.pathsep)))[0]))
    if os.path.exists(outputfilepath):
        print "File {0} already exists".format(outputfilepath)
        sys.exit()

    print "Writing file {0}".format(os.path.basename(outputfilepath))
    with file(outputfilepath, "wt") as of:
        of.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format( 
                 "HGMD ID", "Disease", "Variant Class", "Gene Symbol", "chromosome", "coordinate start", "coordinate end", "strand", "hgvs"))

        diseases = {}
        with file(cosmic_path ,"rU") as f:
            for line in f:
                splitline = line.split("\t")
                cosmic_gene = splitline[0].split("_")[0]
                if cosmic_gene in genes and splitline[7] == "haematopoietic_and_lymphoid_tissue":
                    for disease in reversed(splitline[7:11]):
                        if disease != "NS":
                            break
                    diseases[disease] = ""
                    chrom, startstop = splitline[18].split(":")
                    if chrom in chr_translate:
                        chrom = chr_translate[chrom]
                    start, stop = startstop.split("-")
                    of.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                             ".", disease, ".", cosmic_gene, "chr"+chrom, start, stop, splitline[19], splitline[14]))

        with file(cosmic_path ,"rU") as f:
            for line in f:
                splitline = line.split("\t")
                cosmic_gene = splitline[0].split("_")[0]
                if cosmic_gene not in genes and splitline[7] == "haematopoietic_and_lymphoid_tissue":
                    for disease in reversed(splitline[7:11]):
                        if disease != "NS":
                            break
                    if disease in diseases:
                        chrom, startstop = splitline[18].split(":")
                        if chrom in chr_translate:
                            chrom = chr_translate[chrom]
                        start, stop = startstop.split("-")
                        of.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                                 ".", disease, ".", cosmic_gene, "chr"+chrom, start, stop, splitline[19], splitline[14]))

    if "Disease_Names" in panel:
        with file(panel["Disease_Names"]) as f:
            for l in f:
                if l[0] != "#":
                    old_name, new_name = l.split("=").strip()
                    if old_name in diseases:
                        diseases[old_name] = "="+new_name
        os.unlink(panel["Disease_Names"])

    diseasefilepath = "{0}_Disease_Names.txt".format(os.path.splitext(outputfilepath)[0])
    if os.path.exists(diseasefilepath):
        print "File {0} already exists".format(os.path.basename(diseasefilepath))
        sys.exit()

    print "Writing file {0}".format(os.path.basename(diseasefilepath))
    with file(diseasefilepath, "wt") as f:
        f.write("#Variants Disease Name Translation\n")
        for t in sorted(diseases.items()):
            f.write("".join(t)+"\n")    

    options = panel.read_options()
    options["VariantFrequency"] = "True"
    panel.write_options(options)

if __name__ == "__main__":
    main()



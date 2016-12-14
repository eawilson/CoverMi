#!/usr/bin env python
from panel import Panel
import sys, os, tkFileDialog, Tkinter, re, csv
import pdb



def translate(row):
    if row["Primary site"] == "haematopoietic_and_lymphoid_tissue":
        try:
            return row["Histology subtype 1"]
        except KeyError:
            return row["Histology subtype"] 
    return ""



def main():
    chr_translate = {str(x): "chr{}".format(y) for x, y in zip(range(1, 26), range(1, 23)+["X", "Y", "M"])}

    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print "Please select the COSMIC mutation tsv file"
    cosmic_path = tkFileDialog.askopenfilename(parent=rootwindow, title='Please select the COSMIC mutation file')
    if cosmic_path == "":
        sys.exit()
    with file(cosmic_path, "rU") as f:
        if  "Gene name" not in f.readline().split("\t"):
            abort("File {0} is of incorrect format".format(os.path.basename(cosmic_path)))

    print "Please select a panel"
    panel_path = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a panel')
    if panel_path == "":
        sys.exit()

    panel = Panel(panel_path)
    if "Variants" in panel:
        abort("Variants file {0} already exists within panel".format(os.path.basename(panel["Variants"])))
    if "Targets" not in panel:
        abort("Unable to identify target gene file")
    if "Disease_Names" in panel:
        print "Disease names file {0} will be deleted and updated".format(os.path.basename(panel["Disease_Names"]))
        if raw_input("Press y to continue, any other key to abort") != "y":
            sys.exit()
    target_genes = set(panel.load()["Transcripts"].names)

    outputfilepath = os.path.join(panel_path, 
                              "{0}_{1}.txt".format(os.path.basename(panel_path.rstrip(os.pathsep)), os.path.splitext(os.path.basename(cosmic_path.rstrip(os.pathsep)))[0]))
    if os.path.exists(outputfilepath):
        abort("File {0} already exists".format(outputfilepath))

    print "Writing file {0}".format(os.path.basename(outputfilepath))
    with file(outputfilepath, "wb") as of:
        writer = csv.writer(of, delimiter="\t")
        writer.writerow(["HGMD ID", "Disease", "Variant Class", "Gene Symbol", "chromosome", "coordinate start", "coordinate end", "strand", "hgvs"])

        diseases = {}
        with file(cosmic_path ,"rU") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                gene = row["Gene name"].split("_")[0]
                if gene in target_genes:
                    disease = translate(row)
                    if disease and row["Mutation genome position"]:
                        diseases[disease] = ""
                        chrom, startstop = row["Mutation genome position"].split(":")
                        start, stop = startstop.split("-")
                        writer.writerow([ ".", disease, ".", gene, chr_translate[chrom], start, stop, row["Mutation strand"], row["Mutation AA"]])

        with file(cosmic_path ,"rU") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                gene = row["Gene name"].split("_")[0]
                if gene not in target_genes:
                    disease = translate(row)
                    if disease and row["Mutation genome position"]:
                        if disease in diseases:
                            chrom, startstop = row["Mutation genome position"].split(":")
                            start, stop = startstop.split("-")
                            writer.writerow([ ".", disease, ".", gene, chr_translate[chrom], start, stop, row["Mutation strand"], row["Mutation AA"]])

    if "Disease_Names" in panel:
        with file(panel["Disease_Names"]) as f:
            for l in f:
                if l[0] != "#":
                    old_name, new_name = l.strip().split("=")
                    if old_name in diseases:
                        diseases[old_name] = "="+new_name
        os.unlink(panel["Disease_Names"])

    diseasefilepath = "{0}_Disease_Names.txt".format(os.path.splitext(outputfilepath)[0])
    if os.path.exists(diseasefilepath):
        abort("File {0} already exists".format(os.path.basename(diseasefilepath)))

    print "Writing file {0}".format(os.path.basename(diseasefilepath))
    with file(diseasefilepath, "wt") as f:
        f.write("#Variants Disease Name Translation\n")
        for t in sorted(diseases.items()):
            f.write("".join(t)+"\n")    

    options = panel.read_options()
    options["VariantFrequency"] = "True"
    panel.write_options(options)



def abort(message):
    print message
    raw_input("Press any key to quit")
    sys.exit()



if __name__ == "__main__":
    main()




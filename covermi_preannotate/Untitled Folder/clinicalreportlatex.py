from reportfunctions import TextTable, header
from gr import Gr
import pdb
import os
import time


def textable(body):
    return "\\begin{tabular}{ " + " ".join(["l"]*len(body[0]))+"}\n" + \
           "\\\\\n".join([" & ".join([str(item)for item in row]) for row in body]) + \
           "\\\\\n\\end{tabular}\n"

def texbigtable(headings, body):
    return 





def create(coverage, panel, outputstem):

    tex = [r"""
\documentclass[12pt]{article}
\usepackage{underscore}

\title{CoverMi Report}

\begin{document}
\maketitle
"""]


    if "Exons" not in panel or "Transcripts" not in panel or "Depth" not in panel["Options"]:
        return

    frequency = panel["Options"]["VariantFrequency"]=="True" if ("VariantFrequency" in panel["Options"]) else False

    minimum_depth = panel["Options"]["Depth"]
    if "Amplicons" in panel:
        targeted_range = panel["Amplicons"].merged
        targeted_exons = panel["Exons"].overlapped_by(targeted_range)
    else:
        targeted_range = panel["Transcripts"]
        targeted_exons = panel["Exons"]


    body = []
    for descriptor, item in (("Sample: ", "Sample"), ("Run: ", "Run"), ("Panel: ", "Panel"), ("Manifest: ", "Manifest"), ("Design Studio Bedfile:", "DesignStudio"),
                             ("Known variants file: ", "Variants"), ("Panel type:", "ReportType"), ("Minimum Depth: ", "Depth")):
        for catagory in ("Filenames", "Options"):
            if item in panel[catagory]:
                body += [[descriptor, str(panel[catagory][item])]]
                break
    body += [["CoverMi version: ", "XXX"]]
    tex += [textable(body)]
        

    # Total Coverage
    i = coverage.calculate(targeted_exons, minimum_depth, total=True)
    tex += ["\\subsection*{Overall Coverage}\n"] + \
           ["{0} of {1} bases ({2:0.1f}\\%) covered with a mean depth of {3}\n".format(i.bases_covered, i.bases, i.percent_covered, i.depth_covered)]


 #    # Coverage by Gene
    tex += ["\\subsection*{Coverage by Gene}\n"]
    headers = ["Gene", "Coverage of Targeted Region", "Coverage of Whole Gene", "Mean Depth"]
    body = []
    for i in coverage.calculate(targeted_exons, minimum_depth):
        body += [[i.name, "{}".format(i.percent_covered), "{}".format(i.bases_covered*100.0/panel["Exons"].subset(i.name).base_count), i.depth_covered]]
    if body:
        tex += [textable(body)]




















#    # Summary Table for WGS
#    if "Depths" in panel["Options"]:
#        table = TextTable()
#        table.headers.append(["Depth", "Proportion of genes with"])
#        table.headers.append(["",      "at least 90% coverage"])

#        for depth in sorted(panel["Options"]["Depths"], reverse=True):
#            info = coverage.calculate(targeted_exons, depth)
#            table.rows.append([depth, (float(sum([int(i.percent_covered>=90) for i in info]))*100/len(info), "{:.0f}%")])
#        if len(table.rows) > 0:
#            report += ["\n\n"] + table.formated(sep="    ")

#    # Coverage by Variant per Gene
#    if "Variants_Gene" in panel:
#        table = TextTable()
#        table.headers.append(["Gene", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
#        table.headers.append(["", "in Targeted Region", "in Whole Gene     ", "Sensitivity" if frequency else ""])
#        targeted_variants = panel["Variants_Gene"].subranges_covered_by(targeted_range)
#        for i in coverage.calculate(panel["Variants_Gene"], minimum_depth):
#            detectable = targeted_variants.subset(i.name).number_of_components
#            table.rows.append([i.name, 
#                              [i.components_covered, "/", detectable, "(", (float(i.components_covered)*100/max(detectable,1), "{:.0f}%)")],
#                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
#                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
#        if len(table.rows) > 0:
#           report += ["\n\n"] + table.formated(sep="    ")
#        

#    # Coverage by Variant per Disease
#    if "Variants_Disease" in panel:
#        table = TextTable()
#        table.headers.append(["Disease", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
#        table.headers.append(["", "in Targeted Region", "in Whole Geneome  ", "Sensitivity" if frequency else ""])
#        targeted_variants = panel["Variants_Disease"].subranges_covered_by(targeted_range)
#        for i in coverage.calculate(panel["Variants_Disease"], minimum_depth):
#            detectable = targeted_variants.subset(i.name).number_of_components
#            table.rows.append([i.name,
#                              [i.components_covered, "/", detectable , "(", (float(i.components_covered*100)/max(detectable,1), "{:.0f}%)")],
#                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
#                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
#        if len(table.rows) > 0:
#            report += ["\n\n"] + table.formated(sep="    ")


#    # Coverage by Individual Variant
#    if "Variants_Mutation" in panel:
#        table = TextTable()
#        table.headers.append(["Gene", "Mutation", "Location", "Depth", "Proportion of" if frequency else "", "Disease"])
#        if frequency:
#            table.headers.append(["", "", "", "", "Mutations in" if frequency else "", ""])
#            table.headers.append(["", "", "", "", "Gene" if frequency else "", ""])
#        weighted_mutations_per_gene = {}
#        for entry in panel["Variants_Mutation"].all_entries:
#            gene = entry.name.split()[0]
#            if gene not in weighted_mutations_per_gene:
#                weighted_mutations_per_gene[gene] = 0
#            weighted_mutations_per_gene[gene] += entry.weight
#        for i in coverage.calculate(panel["Variants_Mutation"], minimum_depth):
#            if i.incompletely_covered:
#                table.rows.append([i.name.split()[0],
#                                   i.name.split()[1],
#                                   i.range_combined.locations_as_string,
#                                   i.depth_uncovered,
#                                   (float(i.weighted_components_uncovered)*100/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%") if frequency else "",
#                                   i.diseases])
#        if len(table.rows) > 0:
#            report += ["\n\n"] + ["Inadequately covered targeted variants\n"] 
#            report += table.formated(sep="  ", sortedby=4, reverse=True, trim_columns=((5, 20), (1, 20))) if frequency else table.formated(sep="  ", trim_columns=((4, 20),))


#    # Coverage by Exon
#    table = TextTable()
#    table.headers.append(["Exon", "Coverage of", "Location of", "Mean Depth of"])
#    table.headers.append(["", "Targeted Region", "Uncovered Region", "Uncovered Region"])
#    for i in coverage.calculate(targeted_exons, minimum_depth, exons=True):
#        if i.bases_uncovered > 0:
#            table.rows.append([i.name, 
#                               (i.percent_covered, "{:.0f}%"),
#                               i.range_uncovered.locations_as_string, 
#                               i.depth_uncovered])
#    if len(table.rows) > 0:
#        report += ["\n\n"] + ["Inadequately covered targeted exons\n"] + table.formated(sep="  ")



    tex += [r"""\end{document}
"""]



    
    with file(outputstem+".tex", "wt") as f:
        f.writelines(tex)


import reader
from panel import Panel
from cov import Cov
import pickle

if __name__ == "__main__":

    panel = Panel("/home/ed/Desktop/bampy/Iron Panel v2 - 13.04.15_HD").load()
    #cov = Cov(data=reader.Bam("/home/ed/Desktop/bampy/G165527L_S3.bam").coverage())
    #with open("cov.pickle", "wb") as f:
    #    pickle.dump(cov, f)
    with open("cov.pickle", "rb") as f:
        cov = pickle.load(f)
    create(cov, panel, os.path.join("./textest"))






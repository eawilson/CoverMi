from __future__ import print_function, absolute_import, division
from collections import Counter
from .reportfunctions import TextTable, header
import pdb


def create(coverage, panel, identifier, outputstem):

    if "exons" not in panel or "transcripts" not in panel or "depth" not in panel:
        return

    if "amplicons" in panel:
        targeted_range = panel.amplicons.merged
        targeted_exons = panel.exons.overlapped_by(targeted_range)
    else:
        targeted_range = panel.transcripts
        targeted_exons = panel.exons

    report = header(panel, identifier)
        

    # Total Coverage
    i = coverage.calculate(targeted_exons, panel.depth, total=True)
    report += ["\n\n", "{0} of {1} bases ({2:0.1f}%) covered with a mean depth of {3}\n".format(i.bases_covered, i.bases, i.percent_covered, i.depth_covered)]


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


    # Coverage by Gene
    table = TextTable()
    table.headers.append(["Gene", "Coverage of", "Coverage of", "Mean Depth"])
    table.headers.append(["", "Targeted Region", "Whole Gene     ", ""])
    for i in coverage.calculate(targeted_exons, panel.depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%"), (float(i.bases_covered*100)/panel.exons.subset(lambda e: e.name==i.name).base_count, \
                           "{:.0f}%"), i.depth_covered])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Variant per Gene
    if "variants_gene" in panel:
        frequency = hasattr(next(panel.variants_gene.__iter__()), "weight")
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
        table.headers.append(["", "in Targeted Region", "in Whole Gene     ", "Sensitivity" if frequency else ""])
        targeted_variants = panel.variants_gene.subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel.variants_gene, panel.depth):
            detectable = targeted_variants.subset(lambda e: e.name==i.name).components
            table.rows.append([i.name, 
                              [i.components_covered, "/", detectable, "(", (float(i.components_covered)*100/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")
        

    # Coverage by Variant per Disease
    if "variants_disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
        table.headers.append(["", "in Targeted Region", "in Whole Geneome  ", "Sensitivity" if frequency else ""])
        targeted_variants = panel.variants_disease.subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel.variants_disease, panel.depth):
            detectable = targeted_variants.subset(lambda e: e.name==i.name).components
            table.rows.append([i.name,
                              [i.components_covered, "/", detectable , "(", (float(i.components_covered*100)/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Individual Variant
    if "variants_mutation" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Depth", "Proportion of" if frequency else "", "Disease"])
        if frequency:
            table.headers.append(["", "", "", "", "Mutations in" if frequency else "", ""])
            table.headers.append(["", "", "", "", "Gene" if frequency else "", ""])
            weighted_mutations_per_gene = Counter()
            for entry in panel.variants_mutation.all_entries:
                weighted_mutations_per_gene[entry.name.split()[0]] += entry.weight
        for i in coverage.calculate(panel.variants_mutation, panel.depth):
            if i.incompletely_covered:
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range.locations_as_string,
                                   i.depth_uncovered,
                                   (float(i.weighted_components_uncovered)*100/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%") if frequency else "",
                                   i.diseases])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Inadequately covered targeted variants\n"] 
            report += table.formated(sep="  ", sortedby=4, reverse=True, trim_columns=((5, 20), (1, 20))) if frequency else table.formated(sep="  ", trim_columns=((4, 20),))


    # Coverage by Exon
    table = TextTable()
    table.headers.append(["Exon", "Coverage of", "Location of", "Mean Depth of"])
    table.headers.append(["", "Targeted Region", "Uncovered Region", "Uncovered Region"])
    for i in coverage.calculate(targeted_exons, panel.depth, exons=True):
        if i.bases_uncovered > 0:
            table.rows.append([i.name, 
                               (i.percent_covered, "{:.0f}%"),
                               i.range_uncovered.locations_as_string, 
                               i.depth_uncovered])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Inadequately covered targeted exons\n"] + table.formated(sep="  ")

    
    with open(outputstem+"_covermi_clinical_report.txt", "wt") as f:
        f.writelines(report)


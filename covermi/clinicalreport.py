from reportfunctions import TextTable, header
from gr import Gr
import pdb


def create(coverage, panel, outputstem):

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

    report = header(panel)
        

    # Total Coverage
    i = coverage.calculate(targeted_exons, minimum_depth, total=True)
    report += ["\n\n", "{0} of {1} bases ({2:0.1f}%) covered with a mean depth of {3}\n".format(i.bases_covered, i.bases, i.percent_covered, i.depth_covered)]


    # Summary Table for WGS
    if "Depths" in panel["Options"]:
        table = TextTable()
        table.headers.append(["Depth", "Proportion of genes with"])
        table.headers.append(["",      "at least 90% coverage"])

        for depth in sorted(panel["Options"]["Depths"], reverse=True):
            info = coverage.calculate(targeted_exons, depth)
            table.rows.append([depth, (float(sum([int(i.percent_covered>=90) for i in info]))*100/len(info), "{:.0f}%")])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Gene
    table = TextTable()
    table.headers.append(["Gene", "Coverage of", "Coverage of", "Mean Depth"])
    table.headers.append(["", "Targeted Region", "Whole Gene     ", ""])
    for i in coverage.calculate(targeted_exons, minimum_depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%"), (float(i.bases_covered*100)/panel["Exons"].subset(i.name).base_count, "{:.0f}%"), i.depth_covered])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Variant per Gene
    if "Variants_Gene" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
        table.headers.append(["", "in Targeted Region", "in Whole Gene     ", "Sensitivity" if frequency else ""])
        targeted_variants = panel["Variants_Gene"].subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel["Variants_Gene"], minimum_depth):
            detectable = targeted_variants.subset(i.name).number_of_components
            table.rows.append([i.name, 
                              [i.components_covered, "/", detectable, "(", (float(i.components_covered)*100/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")
        

    # Coverage by Variant per Disease
    if "Variants_Disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Variants Covered", "Clinical" if frequency else ""])
        table.headers.append(["", "in Targeted Region", "in Whole Geneome  ", "Sensitivity" if frequency else ""])
        targeted_variants = panel["Variants_Disease"].subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel["Variants_Disease"], minimum_depth):
            detectable = targeted_variants.subset(i.name).number_of_components
            table.rows.append([i.name,
                              [i.components_covered, "/", detectable , "(", (float(i.components_covered*100)/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Individual Variant
    if "Variants_Mutation" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Depth", "Proportion of" if frequency else "", "Disease"])
        if frequency:
            table.headers.append(["", "", "", "", "Mutations in" if frequency else "", ""])
            table.headers.append(["", "", "", "", "Gene" if frequency else "", ""])
        weighted_mutations_per_gene = {}
        for entry in panel["Variants_Mutation"].all_entries:
            gene = entry.name.split()[0]
            if gene not in weighted_mutations_per_gene:
                weighted_mutations_per_gene[gene] = 0
            weighted_mutations_per_gene[gene] += entry.weight
        for i in coverage.calculate(panel["Variants_Mutation"], minimum_depth):
            if i.incompletely_covered:
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range_combined.locations_as_string,
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
    for i in coverage.calculate(targeted_exons, minimum_depth, exons=True):
        if i.bases_uncovered > 0:
            table.rows.append([i.name, 
                               (i.percent_covered, "{:.0f}%"),
                               i.range_uncovered.locations_as_string, 
                               i.depth_uncovered])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Inadequately covered targeted exons\n"] + table.formated(sep="  ")

    
    with file(outputstem+"_covermi_clinical_report.txt", "wt") as f:
        f.writelines(report)


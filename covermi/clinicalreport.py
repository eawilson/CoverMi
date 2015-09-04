from reportfunctions import TextTable, header
import gr


def create(coverage, panel, outputstem):
    #pdb.set_trace()

    if "Exons" not in panel or "Transcripts" not in panel or "Depth" not in panel:
        return

    minimum_depth = panel["Depth"]
    if "Amplicons" in panel:
        targeted_range = panel["Amplicons"].merged
        targeted_exons = panel["Exons"].overlapped_by(targeted_range)
    else:
        targeted_range = panel["Transcripts"]
        targeted_exons = panel["Exons"]

    report = header(panel)
        
    i = coverage.calculate(targeted_exons, minimum_depth, total=True)
    report += ["\n", "{0} of {1} bases ({2:0.1f}%) covered with a mean depth of {3}\n".format(i.bases_covered, i.bases, i.percent_covered, i.depth_covered)]

    table = TextTable()
    table.headers.append(["Gene", "Coverage of", "Coverage of", "Mean Depth"])
    table.headers.append(["", "Targeted Region", "Whole Gene     ", ""])
    for i in coverage.calculate(targeted_exons, minimum_depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%"), (float(i.bases_covered*100)/panel["Exons"].subset2(i.name).base_count, "{:.0f}%"), i.depth_covered])
    report += ["\n\n"] + table.formated(sep="    ")
    
    if "Variants_Gene" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Variants Covered", "Clinical"])
        table.headers.append(["", "in Targeted Region", "in Whole Gene     ", "Sensitivity"])
        targeted_variants = panel["Variants_Gene"].subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel["Variants_Gene"], minimum_depth):
            detectable = targeted_variants.subset2(i.name).number_of_components
            table.rows.append([i.name, 
                              [i.components_covered, "/", detectable, "(", (float(i.components_covered)*100/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")], 
                              (i.percent_weighted_components_covered, "{:.2f}%")])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")
        
    if "Variants_Disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Variants Covered", "Clinical"])
        table.headers.append(["", "in Targeted Region", "in Whole Geneome  ", "Sensitivity"])
        targeted_variants = panel["Variants_Disease"].subranges_covered_by(targeted_range)
        for i in coverage.calculate(panel["Variants_Disease"], minimum_depth):
            detectable = targeted_variants.subset2(i.name).number_of_components
            table.rows.append([i.name,
                              [i.components_covered, "/", detectable , "(", (float(i.components_covered*100)/max(detectable,1), "{:.0f}%)")],
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%")])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")

    if "Variants_Mutation" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Depth", "Proportion of", "Disease"])
        table.headers.append(["", "", "", "", "Mutations in Gene", ""])
        weighted_mutations_per_gene = {}
        for entry in panel["Variants_Mutation"].all_entries:
            gene = entry[gr.NAME].split()[0]
            if gene not in weighted_mutations_per_gene:
                weighted_mutations_per_gene[gene] = 0
            weighted_mutations_per_gene[gene] += entry[gr.WEIGHT]
        for i in coverage.calculate(panel["Variants_Mutation"].subranges_covered_by(targeted_range), minimum_depth):
            if i.incompletely_covered:
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range_uncovered.locations_as_string,
                                   i.depth_uncovered,
                                   (float(i.weighted_components_uncovered*100)/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%"),
                                   i.diseases])                  
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="  ", sortedby=4, reverse=True)

    table = TextTable()
    table.headers.append(["Gene", "Location", "Depth"])
    for i in coverage.calculate(targeted_exons, minimum_depth, exons=True):
        if i.incompletely_covered:
            table.rows.append([i.name, i.range_uncovered.locations_as_string, i.depth_uncovered])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="  ")

    
    with file(outputstem+"_covermi_clinical_report.txt", "wt") as f:
        f.writelines(report)


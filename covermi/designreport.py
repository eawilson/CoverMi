from reportfunctions import TextTable, header, location
from gr import Gr
import pdb


def create(coverage, panel, outputstem):

    if "Exons" not in panel or "Transcripts" not in panel or "Depth" not in panel["Options"] or "Amplicons" not in panel:
        return

    somatic = panel["Options"]["ReportType"]=="Somatic" if ("ReportType" in panel["Options"]) else False
    frequency = panel["Options"]["VariantFrequency"]=="True" if ("VariantFrequency" in panel["Options"]) else False

    minimum_depth = panel["Options"]["Depth"]

    report = header(panel)


    # Coverage by Gene
    table = TextTable()
    table.headers.append(["Gene", "Coverage"])
    for i in coverage.calculate(panel["Exons"], minimum_depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Gene (for genes not in panel)
    table = TextTable()
    table.headers.append(["Gene", "Coverage"])
    othergenes = panel["AllExons"].touched_by(panel["Amplicons"]).subset(panel["Transcripts"].names, exclude=True, genenames=True)
    if not othergenes.is_empty:
        for i in coverage.calculate(othergenes, minimum_depth):
            table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Coverage of genes not in the targeted panel\n"]  + table.formated(sep="    ")
  

    # Coverage by Variant per Gene
    if "Variants_Gene" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Clinical Sensitivity" if frequency else ""])
        for i in coverage.calculate(panel["Variants_Gene"], minimum_depth):
            table.rows.append([i.name, 
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Variant per Disease
    if "Variants_Disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Clinical Sensitivity" if frequency else ""])
        for i in coverage.calculate(panel["Variants_Disease"], minimum_depth):
            table.rows.append([i.name,
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Individual Variant
    if "Variants_Mutation" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Proportion of" if frequency else "", "Disease"])
        table.headers.append(["", "", "", "Mutations in" if frequency else "", ""])
        table.headers.append(["", "", "", "Gene" if frequency else "", ""])
        weighted_mutations_per_gene = {}
        for entry in panel["Variants_Mutation"].all_entries:
            gene = entry.name.split()[0]
            if gene not in weighted_mutations_per_gene:
                weighted_mutations_per_gene[gene] = 0
            weighted_mutations_per_gene[gene] += entry.weight
        for i in coverage.calculate(panel["Variants_Mutation"], minimum_depth):
            if (somatic and i.completely_covered) or (not somatic and i.incompletely_covered):
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range_combined.locations_as_string,
                                   (float(i.weighted_components_covered if somatic else i.weighted_components_uncovered)*100/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%") if frequency else "",
                                   i.diseases])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Variants "+("" if somatic else "not ")+"covered by panel\n"]
            report += table.formated(sep="  ", sortedby=3, reverse=True, trim_columns=((4, 20), (1,20))) if frequency else table.formated(sep="  ", trim_columns=((3, 20),))


    # Coverage by Exon
    table = TextTable()
    table.headers.append(["Exon", "Coverage", "Covered Region" if somatic else "Uncovered Region"])
    for i in coverage.calculate(panel["Exons"], minimum_depth, exons=True):
        if (somatic and i.bases_covered>0) or (not somatic and i.bases_uncovered>0):
            table.rows.append([i.name, 
                               (i.percent_covered, "{:.0f}%"),
                               i.range_covered.locations_as_string if somatic else i.range_uncovered.locations_as_string])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Exons "+("" if somatic else "not fully ")+"covered by panel\n"] + table.formated(sep="  ")


    rogue_amplicons = panel["Amplicons"].not_touched_by(panel["Transcripts"])
    if not rogue_amplicons.is_empty:
        table = TextTable()
        table.headers.append(["Amplicon", "Location"])
        for chr_name in Gr.KARYOTYPE:
            if chr_name in rogue_amplicons:
                for entry in rogue_amplicons[chr_name]:
                    table.rows.append([entry.name, location(Gr(entry), panel)])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Amplicons not covering a targeted gene\n"] + table.formated(sep="    ")


    if "ExcludedAmplicons" in panel:
        table = TextTable()
        table.headers.append(["Amplicon", "Location"])
        for entry in panel("ExcludedAmplicons").all_entries:
            table.rows.append([entry.name, location(Gr(entry), panel)])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Amplicons in manifest file that have been excluded from analysis\n"] + table.formated(sep="    ")


    table = TextTable()
    table.headers.append(["Exon", "Upstream Padding", "Downstream Padding"])
    for chr_name in panel["Transcripts"]:
        for transcript in panel["Transcripts"][chr_name]:
            exons = panel["Exons"].touched_by(Gr(transcript)).subset(transcript.name)
            amplicons = panel["Amplicons"].touched_by(Gr(transcript)).merged
            loopover = range(0, len(exons[chr_name])) if (transcript.strand=="+") else range(len(exons[chr_name])-1, -1, -1)
            for index in loopover:
                touching_amplicons = amplicons.touched_by(Gr(exons[chr_name][index]))
                if len(touching_amplicons) == 0:
                    continue
                prev_stop = exons[chr_name][index-1].stop if (index>0) else 0
                next_start = exons[chr_name][index+1].start if (index<len(exons[chr_name])-1) else Gr.MAX_CHR_LENGTH
                prev_amp = touching_amplicons[chr_name][0].start
                next_amp = touching_amplicons[chr_name][-1].stop
                padding = [ exons[chr_name][index].start-max(prev_amp, prev_stop)+Gr.SPLICE_SITE_BUFFER,
                            min(next_amp, next_start)-exons[chr_name][index].stop+Gr.SPLICE_SITE_BUFFER ]
                
                table.rows.append([Gr(exons[chr_name][index]).names_as_string,
                                   padding[transcript.strand == "-"] if (padding[transcript.strand == "-"]>0) else "",
                                   padding[transcript.strand == "+"] if (padding[transcript.strand == "+"]>0) else ""])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="  ")

    with file(outputstem+"_covermi_design_report.txt", "wt") as f:
        f.writelines(report)


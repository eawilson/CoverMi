from reportfunctions import TextTable, header, location
import pdb
from gr import Gr

def create(coverage, panel, outputstem):

    if "Exons" not in panel or "Transcripts" not in panel or "Depth" not in panel or "Amplicons" not in panel:
        return

    minimum_depth = panel["Depth"]
    if "Manifest" in panel["Filenames"]:
        manifest = "{0:25} {1}".format("Manifest:", panel["Filenames"]["Manifest"])
    elif "DesignStudio" in panel["Filenames"]:
        manifest = "{0:25} {1}".format("Design Studio bedfile:", panel["Filenames"]["DesignStudio"])
    else:
        raise CovermiException("No manifest or design studio bedfile in panel")

    report = header(panel)

    table = TextTable()
    table.headers.append(["Genes in Panel", "Coverage"])
    for i in coverage.calculate(panel["Exons"], minimum_depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    report += ["\n\n"] + table.formated(sep="    ")

    table = TextTable()
    table.headers.append(["Genes not in Panel", "Coverage"])
    othergenes = panel["AllExons"].touched_by(panel["Amplicons"]).subset2(panel["Transcripts"].names, exclude=True, genenames=True)
    if not othergenes.is_empty:
        for i in coverage.calculate(othergenes, minimum_depth):
            table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    report += ["\n\n"] + table.formated(sep="    ")

    table = TextTable()
    table.headers.append(["Exons in Panel", "Coverage", "Uncovered Region"])
    for i in coverage.calculate(panel["Exons"], minimum_depth, exons=True):
        if i.incompletely_covered:
            table.rows.append([i.name, (i.percent_covered, "{:.0f}%"), i.range_uncovered.locations_as_string])
    report += ["\n\n"] + table.formated(sep="    ")
  
    if "Variants_Gene" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Clinical Sensitivity"])
        for i in coverage.calculate(panel["Variants_Gene"], minimum_depth):
            table.rows.append([i.name, 
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%")])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")

    if "Variants_Disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Clinical Sensitivity"])
        for i in coverage.calculate(panel["Variants_Disease"], minimum_depth):
            table.rows.append([i.name,
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%")])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")

    if "Variants_Mutation" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Proportion of", "Disease"])
        table.headers.append(["", "", "", "Mutations in Gene", ""])
        for i in coverage.calculate(panel["Variants_Mutation"], minimum_depth):
            if i.incompletely_covered:
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range_uncovered.locations_as_string,
                                   (i.percent_weighted_components_covered, "{:.2f}%"),
                                   i.diseases])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="  ", sortedby=4, reverse=True)

    rogue_amplicons = panel["Amplicons"].not_touched_by(panel["Transcripts"])
    if not rogue_amplicons.is_empty:
        table = TextTable()
        table.headers.append(["Amplicons Not Covering", "Location"])
        table.headers.append(["a Panel Gene", ""])
        for chr_name in Gr.KARYOTYPE:
            if chr_name in rogue_amplicons:
                for entry in rogue_amplicons[chr_name]:
                    table.rows.append([entry[Gr.NAME], location(Gr(entry), panel)])
        report += ["\n\n"] + table.formated(sep="  ")

    if "Excluded" in panel:
        table = TextTable()
        table.headers.append(["Excluded Amplicons", "Location"])
        excluded_amplicons = panel["AllAmplicons"].subset2(panel["Excluded"])
        for chr_name in Gr.KARYOTYPE:
            if chr_name in excluded_amplicons:
                for entry in excluded_amplicons[chr_name]:
                    table.rows.append([entry[Gr.NAME], location(Gr(entry), panel)])
        report += ["\n\n"] + table.formated(sep="  ")

    table = TextTable()
    table.headers.append(["Exon", "Upstream Padding", "Downstream Padding"])
    for chr_name in panel["Transcripts"]:
        for transcript in panel["Transcripts"][chr_name]:
            exons = panel["Exons"].touched_by(Gr(transcript)).subset2(transcript[Gr.NAME])
            amplicons = panel["Amplicons"].touched_by(Gr(transcript)).merged
            loopover = range(0, len(exons[chr_name])) if (transcript[Gr.STRAND]=="+") else range(len(exons[chr_name])-1, -1, -1)
            for index in loopover:
                touching_amplicons = amplicons.touched_by(Gr(exons[chr_name][index]))
                if len(touching_amplicons) == 0:
                    continue
                prev_stop = exons[chr_name][index-1][Gr.STOP] if (index>0) else 0
                next_start = exons[chr_name][index+1][Gr.START] if (index<len(exons[chr_name])-1) else Gr.MAX_CHR_LENGTH
                prev_amp = touching_amplicons[chr_name][0][Gr.START]
                next_amp = touching_amplicons[chr_name][-1][Gr.STOP]
                padding = [ exons[chr_name][index][Gr.START]-max(prev_amp, prev_stop)+Gr.SPLICE_SITE_BUFFER,
                            min(next_amp, next_start)-exons[chr_name][index][Gr.STOP]+Gr.SPLICE_SITE_BUFFER ]
                
                table.rows.append([Gr(exons[chr_name][index]).names_as_string,
                                   padding[transcript[Gr.STRAND] == "-"] if (padding[transcript[Gr.STRAND] == "-"]>0) else "",
                                   padding[transcript[Gr.STRAND] == "+"] if (padding[transcript[Gr.STRAND] == "+"]>0) else ""])
    report += ["\n\n"] + table.formated(sep="  ")

    with file(outputstem+"_covermi_design_report.txt", "wt") as f:
        f.writelines(report)


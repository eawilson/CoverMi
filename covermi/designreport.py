from __future__ import print_function, absolute_import, division

from .reportfunctions import TextTable, header, location
from .gr import Gr, SPLICE_SITE_BUFFER, MAX_LEN
import pdb


def create(coverage, panel, identifier, outputstem):

    if "exons" not in panel or "transcripts" not in panel or "depth" not in panel or "amplicons" not in panel:
        return

    somatic = panel.properties.get("reporttype", "constitutional").lower() == "somatic"
    report = header(panel, identifier)

    # Coverage by Gene
    table = TextTable()
    table.headers.append(["Gene", "Coverage"])
    for i in coverage.calculate(panel.exons,panel.depth):
        table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    if len(table.rows) > 0:
        report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Gene (for genes not in panel)
    table = TextTable()
    table.headers.append(["Gene", "Coverage"])
    targeted_names = panel.transcripts.names
    othergenes = panel.allexons.touched_by(panel.amplicons).subset(lambda e: e.name not in targeted_names)
    if not othergenes.is_empty:
        for i in coverage.calculate(othergenes, panel.depth):
            table.rows.append([i.name, (i.percent_covered, "{:.0f}%")])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Coverage of genes not in the targeted panel\n"]  + table.formated(sep="    ")
  

    # Coverage by Variant per Gene
    if "variants_gene" in panel:
        frequency = hasattr(next(panel.variants_gene.__iter__()), "weight")
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered", "Clinical Sensitivity" if frequency else ""])
        for i in coverage.calculate(panel.variants_gene,panel.depth):
            table.rows.append([i.name, 
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
           report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Variant per Disease
    if "variants_disease" in panel:
        table = TextTable()
        table.headers.append(["Disease", "Variants Covered", "Clinical Sensitivity" if frequency else ""])
        for i in coverage.calculate(panel.variants_disease,panel.depth):
            table.rows.append([i.name,
                              [i.components_covered, "/", i.components, "(", (i.percent_components_covered,  "{:.0f}%)")],
                              (i.percent_weighted_components_covered, "{:.0f}%") if frequency else ""])
        if len(table.rows) > 0:
            report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Individual Variant
    if "variants_mutation" in panel:
        frequency = hasattr(next(panel.variants_gene.__iter__()), "weight")
        table = TextTable()
        table.headers.append(["Gene", "Mutation", "Location", "Proportion of" if frequency else "", "Disease"])
        table.headers.append(["", "", "", "Mutations in" if frequency else "", ""])
        table.headers.append(["", "", "", "Gene" if frequency else "", ""])
        if frequency:
            weighted_mutations_per_gene = {}
            for entry in panel.variants_mutation.all_entries:
                gene = entry.name.split()[0]
                if gene not in weighted_mutations_per_gene:
                    weighted_mutations_per_gene[gene] = 0
                weighted_mutations_per_gene[gene] += entry.weight
        for i in coverage.calculate(panel.variants_mutation,panel.depth):
            if (somatic and i.completely_covered) or (not somatic and i.incompletely_covered):
                table.rows.append([i.name.split()[0],
                                   i.name.split()[1],
                                   i.range.locations_as_string,
                                   (float(i.weighted_components_covered if somatic else i.weighted_components_uncovered)*100/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%") if frequency else "",
                                   i.diseases])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Variants "+("" if somatic else "not ")+"covered by panel\n"]
            report += table.formated(sep="  ", sortedby=3, reverse=True, trim_columns=((4, 20), (1,20))) if frequency else table.formated(sep="  ", trim_columns=((3, 20),))


    # Coverage by Exon
    table = TextTable()
    table.headers.append(["Exon", "Coverage", "Covered Region" if somatic else "Uncovered Region"])

    for i in coverage.calculate(panel.exons, panel.depth, exons=True):
        if (somatic and i.bases_covered>0) or (not somatic and i.bases_uncovered>0):
            table.rows.append([i.name, 
                               (i.percent_covered, "{:.0f}%"),
                               i.range_covered.locations_as_string if somatic else i.range_uncovered.locations_as_string])
    if len(table.rows) > 0:
        report += ["\n\n"] + ["Exons "+("" if somatic else "not fully ")+"covered by panel\n"] + table.formated(sep="  ")


    rogue_amplicons = panel.amplicons.not_touched_by(panel.transcripts)
    if not rogue_amplicons.is_empty:
        table = TextTable()
        table.headers.append(["Amplicon", "Location"])
        for entry in rogue_amplicons.all_entries:
            table.rows.append([entry.name, location(Gr().add(entry), panel)])
        if len(table.rows) > 0:
            report += ["\n\n"] + ["Amplicons not covering a targeted gene\n"] + table.formated(sep="    ")


#    if "ExcludedAmplicons" in panel:
#        table = TextTable()
#        table.headers.append(["Amplicon", "Location"])
#        for entry in panel("ExcludedAmplicons").all_entries:
#            table.rows.append([entry.name, location(Gr().add(entry), panel)])
#        if len(table.rows) > 0:
#            report += ["\n\n"] + ["Amplicons in manifest file that have been excluded from analysis\n"] + table.formated(sep="    ")


    table = TextTable()
    table.headers.append(["Exon", "Upstream Padding", "Downstream Padding"])
    for chr_name in panel.transcripts.data:
        for transcript in panel.transcripts.data[chr_name]:
            exons = panel.exons.touched_by(Gr().add(transcript)).subset(lambda e: e.name==transcript.name)
            amplicons = panel.amplicons.touched_by(Gr().add(transcript)).merged
            loopover = range(0, len(exons.data[chr_name])) if (transcript.strand=="+") else range(len(exons.data[chr_name])-1, -1, -1)
            for index in loopover:
                touching_amplicons = amplicons.touched_by(Gr().add(exons.data[chr_name][index]))
                if touching_amplicons.is_empty:
                    continue
                prev_stop = exons.data[chr_name][index-1].stop if (index>0) else 0
                next_start = exons.data[chr_name][index+1].start if (index<len(exons.data[chr_name])-1) else MAX_LEN
                prev_amp = touching_amplicons.data[chr_name][0].start
                next_amp = touching_amplicons.data[chr_name][-1].stop
                padding = [ exons.data[chr_name][index].start-max(prev_amp, prev_stop)+SPLICE_SITE_BUFFER,
                            min(next_amp, next_start)-exons.data[chr_name][index].stop+SPLICE_SITE_BUFFER ]
                
                table.rows.append([Gr().add(exons.data[chr_name][index]).names_as_string,
                                   padding[transcript.strand == "-"] if (padding[transcript.strand == "-"]>0) else "",
                                   padding[transcript.strand == "+"] if (padding[transcript.strand == "+"]>0) else ""])
    if len(table.rows) > 0:
        table.rows.sort()
        report += ["\n\n"] + table.formated(sep="  ")

    with open(outputstem+"_covermi_design_report.txt", "wt") as f:
        f.writelines(report)


import time
import pdb
import os
import argparse
from itertools import chain, repeat
from copy import copy

from covermi import Panel, Gr, Cov, Plot, __version__



class TextTable(object):

    def __init__(self):
        self.rows = []
        self.headers = []

    @classmethod
    def _aligned(cls, table, delete_empty_cols=True):
        columns = len(table[0])
        rows = len(table)
        newtab = [[] for n in range(0, rows)]

        for col in range(0, columns):
            if type(table[0][col]) == list:
                replacement = cls._aligned([table[row][col] for row in range(0, rows)])
                for row in range(0, rows):
                    newtab[row].append("".join(replacement[row]))
            else:
                if type(table[0][col]) == tuple:
                    fstring = table[0][col][1]
                    for row in range(0, rows):
                        table[row][col] = table[row][col][0]
                else:
                    fstring = "{}"

                biggest = max([len(fstring.format(table[row][col])) for row in range(0, rows)])
                if biggest >0 or not delete_empty_cols:
                    for row in range(0, rows):
                        newtab[row].append(fstring.format(table[row][col]).ljust(biggest) if (type(table[row][col])==str) else fstring.format(table[row][col]).rjust(biggest))
        return newtab


    def formated(self, sep="", sortedby=None, reverse=False, maxcolwidth=40):
        if sortedby is not None:
            def bycolumn(line): return [line[sortedby]] + line
            self.rows.sort(key=bycolumn, reverse=reverse)

        if maxcolwidth is not None: 
            for row in range(0, len(self.rows)):
                self.rows[row] = [item[0:maxcolwidth-2]+".." if (type(item)==str and len(item)>maxcolwidth) else item for item in self.rows[row]]

        self.rows = type(self)._aligned(self.rows, delete_empty_cols=(len(self.headers)==0))
        if len(self.headers) > 0:
            table = [sep.join(row)+"\n" for row in type(self)._aligned(self.headers + self.rows)]
            table = [("-"*(len(table[0])-1))+"\n"] + table[0:len(self.headers)] + [("-"*(len(table[0])-1))+"\n"] + table[len(self.headers):]
        else:
            table = [sep.join(row)+"\n" for row in self.rows]
        return table



#def header(panel):
    #table = TextTable()
    #for descriptor, item in (("Sample: ", "Sample"), ("Run: ", "Run"), ("Panel: ", "Panel"), ("Manifest: ", "Manifest"), ("Design Studio Bedfile:", "DesignStudio"),
                             #("Known variants file: ", "Variants"), ("Panel type:", "ReportType"), ("Minimum Depth: ", "Depth")):
        #for catagory in ("Filenames", "Options"):
            #if item in panel[catagory]:
                #table.rows.append([descriptor, str(panel[catagory][item])])
                #break
    #table.rows += [["Date of CoverMi analysis: ", time.strftime("%d/%m/%Y")], ["CoverMi version: ", __version__]]
    #return table.formated()



def designreport():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--panel", help="Path to covermi panel which must contain targets bedfile.", required=True)
    parser.add_argument("-o", "--output", help="Directory to write output files to.", default=".")
    args = parser.parse_args()
    
    title = os.path.basename(args.panel)
    panel = Panel(args.panel)
    output=os.path.join(args.output, f"{title}_design")
    
    if "exons" not in panel or "transcripts" not in panel or "targets" not in panel:
        sys.exit("Panel insufficient to create a design report.")
    
    cov = Cov(reads=chain(*repeat(panel["targets"].merged, 500)))
    
    Plot(panel=panel, coverage=cov, title=title, output=f"{output}.pdf")
    
    min_depth = 400
    report = []
    
    # Coverage by Gene
    table = TextTable()
    table.headers.append(["Gene", "Transcript", "Coverage"])
    for i in cov.calculate(panel["exons"], min_depth):
        gene, transcript = i.name.split()
        table.rows.append([gene, transcript, (i.percent_covered, "{:.0f}%")])
    report += ["\n\n"] + table.formated(sep="    ")
    
    
    # Coverage by Variant by Gene
    targeted_gene = set(name.split()[0] for name in panel["names"])
    if "variants" in panel:
        other_components_covered = 0
        other_components = 0
        table = TextTable()
        table.headers.append(["Gene", "Variants Covered"])
        for i in cov.calculate(panel["variants"], min_depth):
            if i.name in targeted_gene:
                table.rows.append([str(i.name), [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")]])
            else:
                other_components_covered += i.components_covered
                other_components += i.components
        table.rows.append(["Other", [other_components_covered, "/", other_components, "(", (other_components_covered * 100.0 / other_components, "{:.0f}%)")]])
        report += ["\n\n"] + table.formated(sep="    ")
    
    
    # Coverage by Variant by Site
    if "variants" in panel:
        table = TextTable()
        table.headers.append(["Site", "Variants Covered"])
        variants = Gr()
        for entry in panel["variants"]:
            entry = copy(entry)
            try:
                entry.name = entry.name.site
            except AttributeError:
                continue
            variants.add(entry)
        
        for i in cov.calculate(variants, min_depth):
            table.rows.append([str(i.name), [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")]])
        report += ["\n\n"] + table.formated(sep="    ")


    # Coverage by Variant by Histology
    if "variants" in panel:
        table = TextTable()
        table.headers.append(["Histology", "Variants Covered"])
        variants = Gr()
        for entry in panel["variants"]:
            entry = copy(entry)
            try:
                entry.name = entry.name.histology
            except AttributeError:
                continue
            variants.add(entry)
        
        for i in cov.calculate(variants, min_depth):
            table.rows.append([str(i.name), [i.components_covered, "/", i.components, "(", (i.percent_components_covered, "{:.0f}%)")]])
        report += ["\n\n"] + table.formated(sep="    ")
    
    
    
    ## Coverage by Exon
    #table = TextTable()
    #table.headers.append(["Gene", "Transcript", "Exon", "Coverage"])
    #for exon in panel["exons"]:
        
    #for i in cov.calculate(, min_depth):
        #gene, transcript = i.name.split()
        #table.rows.append([gene, transcript, (i.percent_covered, "{:.0f}%")])
    #report += ["\n\n"] + table.formated(sep="    ")
    
    
    
    # Coverage by Individual Variant
    if "variants" in panel:
        table = TextTable()
        table.headers.append(["Gene", "Transcript", "Covered Variants"])
        variants = Gr()
        dedup = set()
        for entry in panel["variants"].covered_by(panel["targets"]):
            try:
                cds = entry.name.cds
                transcript = entry.name.transcript
            except AttributeError:
                continue
            if (transcript, cds) not in dedup:
                dedup.add((transcript, cds))
                table.rows.append([str(entry.name), transcript, cds])
        report += ["\n\n"] + table.formated(sep="    ")
    
    
    
    ## Coverage by Individual Variant
    #if "Variants_Mutation" in panel:
        #table = TextTable()
        #table.headers.append(["Gene", "Mutation", "Location", "Proportion of" if frequency else "", "Disease"])
        #table.headers.append(["", "", "", "Mutations in Gene" if frequency else "", ""])
        #weighted_mutations_per_gene = {}
        #for entry in panel["Variants_Mutation"].all_entries:
            #gene = entry.name.split()[0]
            #if gene not in weighted_mutations_per_gene:
                #weighted_mutations_per_gene[gene] = 0
            #weighted_mutations_per_gene[gene] += entry.weight
        #for i in coverage.calculate(panel["Variants_Mutation"], min_depth):
            #if (somatic and i.bases_covered>0) or (not somatic and i.bases_uncovered>0):
                #table.rows.append([i.name.split()[0],
                                   #i.name.split()[1],
                                   #i.range_covered.locations_as_string if somatic else i.range_uncovered.locations_as_string,
                                   #(float(i.weighted_components_covered if somatic else i.weighted_components_uncovered)*100/weighted_mutations_per_gene[i.name.split()[0]], "{:.2f}%") if frequency else "",
                                   #i.diseases])
        #if len(table.rows) > 0:
            #report += ["\n\n"] + ["Variants "+("" if somatic else "not ")+"covered by panel\n"]
            #report += table.formated(sep="  ", sortedby=3, reverse=True) if frequency else table.formated(sep="  ")


    ## Coverage by Exon
    #table = TextTable()
    #table.headers.append(["Exon", "Coverage", "Covered Region" if somatic else "Uncovered Region"])
    #for i in coverage.calculate(panel["exons"], min_depth, exons=True):
        #if (somatic and i.bases_covered>0) or (not somatic and i.bases_uncovered>0):
            #table.rows.append([i.name, 
                               #(i.percent_covered, "{:.0f}%"),
                               #i.range_covered.locations_as_string if somatic else i.range_uncovered.locations_as_string])
    #if len(table.rows) > 0:
        #report += ["\n\n"] + ["Exons "+("" if somatic else "not fully ")+"covered by panel\n"] + table.formated(sep="  ")


    #rogue_amplicons = panel["targets"].not_touched_by(panel["transcripts"])
    #if not rogue_amplicons.is_empty:
        #table = TextTable()
        #table.headers.append(["Amplicon", "Location"])
        #for chr_name in Gr.KARYOTYPE:
            #if chr_name in rogue_amplicons:
                #for entry in rogue_amplicons[chr_name]:
                    #table.rows.append([entry.name, location(Gr(entry), panel)])
        #if len(table.rows) > 0:
            #report += ["\n\n"] + ["Targets not covering a targeted gene\n"] + table.formated(sep="    ")


    #if "ExcludedAmplicons" in panel:
        #table = TextTable()
        #table.headers.append(["Amplicon", "Location"])
        #for entry in panel("ExcludedAmplicons").all_entries:
            #table.rows.append([entry.name, location(Gr(entry), panel)])
        #if len(table.rows) > 0:
            #report += ["\n\n"] + ["Targets in manifest file that have been excluded from analysis\n"] + table.formated(sep="    ")


    #table = TextTable()
    #table.headers.append(["Exon", "Upstream Padding", "Downstream Padding"])
    #for chr_name in panel["transcripts"]:
        #for transcript in panel["transcripts"][chr_name]:
            #exons = panel["exons"].touched_by(Gr(transcript)).subset2(transcript.name)
            #amplicons = panel["targets"].touched_by(Gr(transcript)).merged
            #loopover = range(0, len(exons[chr_name])) if (transcript.strand=="+") else range(len(exons[chr_name])-1, -1, -1)
            #for index in loopover:
                #touching_amplicons = amplicons.touched_by(Gr(exons[chr_name][index]))
                #if len(touching_amplicons) == 0:
                    #continue
                #prev_stop = exons[chr_name][index-1].stop if (index>0) else 0
                #next_start = exons[chr_name][index+1].start if (index<len(exons[chr_name])-1) else Gr.MAX_CHR_LENGTH
                #prev_amp = touching_amplicons[chr_name][0].start
                #next_amp = touching_amplicons[chr_name][-1].stop
                #padding = [ exons[chr_name][index].start-max(prev_amp, prev_stop)+Gr.SPLICE_SITE_BUFFER,
                            #min(next_amp, next_start)-exons[chr_name][index].stop+Gr.SPLICE_SITE_BUFFER ]
                
                #table.rows.append([Gr(exons[chr_name][index]).names_as_string,
                                   #padding[transcript.strand == "-"] if (padding[transcript.strand == "-"]>0) else "",
                                   #padding[transcript.strand == "+"] if (padding[transcript.strand == "+"]>0) else ""])
    #if len(table.rows) > 0:
        #report += ["\n\n"] + table.formated(sep="  ")

    with open(f"{output}.txt", "wt") as f_out:
        f_out.writelines(report)



def main():
    try:
        designreport()
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

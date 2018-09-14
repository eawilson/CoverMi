from __future__ import print_function, absolute_import, division


from .reportfunctions import TextTable, header, location


def create(info, panel, identifier, outputstem):

    if "depth" not in panel:
        return

    with open(outputstem+"_amplicon_data.tsv", "wt") as f:
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
        f.write(template.format("identifier", "run", "amplicon", "chr", "fdepth", "rdepth"))
        for i in info:
            f.write(template.format(identifier.name, identifier.run, i.name, i.chrom, i.f_depth, i.r_depth))

    report = header(panel, identifier)

    passed = 0
    mean_depth = 0 
    for i in info:
        if i.min_depth >= panel.depth:
            passed += 1
            mean_depth += i.mean_depth
    mean_depth = ", (mean depth {0})".format(mean_depth//passed) if (passed>0) else ""
    report += ["\n", "{0} of {1} amplicons ({2:0.1f}%) provided adequate coverage{3}\n".format(passed, len(info), float(passed)/len(info)*100, mean_depth)]
  
    table = TextTable()
    table.headers.append(["Amplicon", "Gene", "Location", "Forward", "Reverse"])
    table.headers.append(["", "", "", "Depth", "Depth"])
    for i in info:
        if i.min_depth < panel.depth:
            table.rows.append([i.name, location(i.gr, panel), i.gr.locations_as_string, i.f_depth, i.r_depth])
    if len(table.rows) > 0:
        report += ["\n\n", "The following amplicons failed to provide adequate coverage\n"] + table.formated(sep="   ")
 
    table = TextTable()
    table.headers.append(["Amplicon", "Gene", "Location", "Forward", "Reverse", "Ratio"])
    table.headers.append(["", "", "", "Depth", "Depth", ""])
    for i in info:
        if  i.ratio < 0.7 and i.max_depth > panel.depth:
            table.rows.append([i.name, location(i.gr, panel), i.gr.locations_as_string, i.f_depth, i.r_depth, (i.ratio, "{:.2f}")])
    if len(table.rows) > 0:
        report += ["\n\n", "The following amplicons had an abnormal forward:reverse depth ratio (<0.7)\n"] + table.formated(sep="   ")
    
    with open(outputstem+"_covermi_technical_report.txt", "wt") as f:
        f.writelines(report)


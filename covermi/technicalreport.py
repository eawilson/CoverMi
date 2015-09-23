from reportfunctions import TextTable, header, location


def create(info, panel, outputstem):

    if "Depth" not in panel["Options"]:
        return
    minimum_depth = panel["Options"]["Depth"]

    with file(outputstem+"_amplicon_data.tsv", "wt") as f:
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
        f.write(template.format("sample", "run", "amplicon", "chr", "fdepth", "rdepth"))
        for i in info:
            f.write(template.format(panel["Filenames"]["Sample"], panel["Filenames"]["Run"], i.name, i.chrom, i.f_depth, i.r_depth))

    report = header(panel)

    passed = 0
    mean_depth = 0 
    for i in info:
        if i.minimum_depth >= minimum_depth:
            passed += 1
            mean_depth += i.mean_depth     
    mean_depth = ", (mean depth {0})".format(mean_depth/passed) if (passed>0) else ""
    report += ["\n", "{0} of {1} amplicons ({2:0.1f}%) provided adequate coverage{3}\n".format(passed, len(info), float(passed)/len(info)*100, mean_depth)]
  
    table = TextTable()
    table.headers.append(["Amplicon", "Gene", "Location", "Forward", "Reverse"])
    table.headers.append(["", "", "", "Depth", "Depth"])
    for i in info:
        if i.minimum_depth < minimum_depth:
            table.rows.append([i.name, location(i.gr, panel), i.gr.locations_as_string, i.f_depth, i.r_depth])
    if len(table.rows) > 0:
        report += ["\n\n", "The following amplicons failed to provide adequate coverage\n"] + table.formated(sep="   ")
 
    table = TextTable()
    table.headers.append(["Amplicon", "Gene", "Location", "Forward", "Reverse", "Ratio"])
    table.headers.append(["", "", "", "Depth", "Depth", ""])
    for i in info:
        if  i.ratio < 0.7 and i.maximum_depth > minimum_depth:
            table.rows.append([i.name, location(i.gr, panel), i.gr.locations_as_string, i.f_depth, i.r_depth, (i.ratio, "{:.2f}")])
    if len(table.rows) > 0:
        report += ["\n\n", "The following amplicons had an abnormal forward:reverse depth ratio (<0.7)\n"] + table.formated(sep="   ")
    
    with file(outputstem+"_covermi_technical_report.txt", "wt") as f:
        f.writelines(report)


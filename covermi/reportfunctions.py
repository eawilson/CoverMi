import time, pdb, pkg_resources


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
                if not (biggest == 0 and delete_empty_cols):
                    for row in range(0, rows):
                        newtab[row].append(fstring.format(table[row][col]).ljust(biggest) if (type(table[row][col])==str) else fstring.format(table[row][col]).rjust(biggest))
        return newtab


    def formated(self, sep="", sortedby=None, reverse=False, maxwidth=99, trimcolumn=None):
        if sortedby is not None:
            def bycolumn(line): return [line[sortedby]] + line
            self.rows.sort(key=bycolumn, reverse=reverse)

        rows = type(self)._aligned(self.rows, delete_empty_cols=(len(self.headers)==0))
        if len(self.headers) > 0:
            rows = type(self)._aligned(self.headers + rows)

        if trimcolumn is not None:
            excess_width = len(sep.join(rows[0])) - maxwidth
            if excess_width > 0:
                trim_len = len(rows[0][trimcolumn]) - (excess_width + 2)
                if trim_len > 0:
                    for row in rows:
                        row[trimcolumn] = row[trimcolumn][0:trim_len] + ".."

        if len(self.headers) > 0:
            table = [sep.join(row)+"\n" for row in rows]
            table = [("-"*(len(table[0])-1))+"\n"] + table[0:len(self.headers)] + [("-"*(len(table[0])-1))+"\n"] + table[len(self.headers):]
        else:
            table = [sep.join(row)+"\n" for row in rows]
        return table


def header(panel):
    table = TextTable()
    for descriptor, item in (("Sample: ", "Sample"), ("Run: ", "Run"), ("Panel: ", "Panel"), ("Manifest: ", "Manifest"), ("Design Studio Bedfile:", "DesignStudio"),
                             ("Known variants file: ", "Variants"), ("Panel type:", "ReportType"), ("Minimum Depth: ", "Depth")):
        for catagory in ("Filenames", "Options"):
            if item in panel[catagory]:
                table.rows.append([descriptor, str(panel[catagory][item])])
                break
    table.rows += [["Date of CoverMi analysis: ", time.strftime("%d/%m/%Y")], ["CoverMi version: ", pkg_resources.require("CoverMi")[0].version]]
    return table.formated()


def location(gr1, panel):
    if "Exons" and "Transcripts" in panel:
        exons = panel["AllExons"] if ("AllExons" in panel) else panel["Exons"]
        transcripts = panel["AllTranscripts"] if ("AllTranscripts" in panel) else panel["Transcripts"]
        loc = exons.overlapped_by(gr1).names_as_string
        if loc =="":
            loc = "{0} intron".format(transcripts.overlapped_by(gr1).names_as_string)
            if loc == " intron":
                loc = "not covering a "+("" if ("AllTranscripts" in panel) else"targeted ")+"gene"
    else:
        loc = ""
    return loc


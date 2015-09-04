import time, pdb

def _rjust(s, l):
    return s.rjust(l, " ")

def _ljust(s, l):
    return s.ljust(l, " ")



class Table(list):
    rows = [] 
    headers = []

    def align(self):
        justify = []
        columns = len(self.rows[0])
        for col in range(0, columns):
            justify.append(_ljust if (type(self.rows[0][col])==str) else _rjust)
        for col in range(0, columns):
            biggest = 0
            for row in range(0, len(self.rows)):
                self.rows[row][col] = self.format[col].format(self.rows[row][col])
                biggest = max(biggest, len(self.rows[row][col]))
            for row in range(0, len(self.rows)):
                self.rows[row][col] = justify[col](self.rows[row][col], biggest)

        if len(self.headers) > 0:
            for col in range(0, columns):     
                biggest = max([len(self.headers[row][col]) for row in range(0, len(self.headers))] + [len(self.rows[0][col])])
                for row in range(0, len(self.headers)):
                    self.headers[row][col] = self.headers[row][col].ljust(biggest, " ")
                for row in range(0, len(self.rows)):
                    self.rows[row][col] = self.rows[row][col].ljust(biggest, " ")
            seperator = "  "     
            output = [seperator.join(row)+"\n" for row in self.headers]
            output += "-"*(len(output[0])-1)+"\n"
        else:
            seperator = ""
            output = []

        return output + [seperator.join(row)+"\n" for row in self.rows]

def align(table):
    #pdb.set_trace()
    columns = len(table[0])
    rows = len(table)
    newtab = [[] for n in range(0, rows)]

    for col in range(0, columns):
        if type(table[0][col]) == list:
           replacement = align([table[row][col] for row in range(0, rows)])
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
            for row in range(0, rows):
                isstring = (type(table[row][col])==str)
                newtab[row].append(fstring.format(table[row][col]).ljust(biggest) if isstring else fstring.format(table[row][col]).rjust(biggest))
    return newtab


def formattable(header, table, sep=""):
    #pdb.set_trace()
    table = align(table)
    table = [sep.join(row)+"\n" for row in align(header + table)]
    table = table[0:len(header)] + [("-"*(len(table[0])-1))+"\n"] + table[len(header):]
    return table


def formattext(text, sep=""):
    return [sep.join(row)+"\n" for row in align(text)]


class TextTable(object):

    def __init__(self):
        self.rows = []
        self.headers = []

    @staticmethod
    def _aligned(table):
        columns = len(table[0])
        rows = len(table)
        newtab = [[] for n in range(0, rows)]

        for col in range(0, columns):
            if type(table[0][col]) == list:
                replacement = TextTable._aligned([table[row][col] for row in range(0, rows)])
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
                for row in range(0, rows):
                    newtab[row].append(fstring.format(table[row][col]).ljust(biggest) if (type(table[row][col])==str) else fstring.format(table[row][col]).rjust(biggest))
        return newtab


    def formated(self, sep="", sortedby=None, reverse=False, maxcolwidth=40):
        if sortedby is not None:
            def bycolumn(line): return [line[sortedby]] + line
            self.rows.sort(key=bycolumn, reverse=reverse)

        if maxcolwidth is not None: 
            for row in range(0, len(self.rows)):
                self.rows[row] = [item[0:38]+".." if (type(item)==str and len(item)>40) else item for item in self.rows[row]]

        self.rows = TextTable._aligned(self.rows)
        if len(self.headers) > 0:
            table = [sep.join(row)+"\n" for row in TextTable._aligned(self.headers + self.rows)]
            table = table[0:len(self.headers)] + [("-"*(len(table[0])-1))+"\n"] + table[len(self.headers):]
        else:
            table = [sep.join(row)+"\n" for row in self.rows]
        return table


def header(panel):
    table = TextTable()
    for descriptor, item in (("Sample: ", "Sample"), ("Run: ", "Run"), ("Panel: ", "Panel"), ("Manifest: ", "Manifest"), ("Design Studio Bedfile:", "DesignStudio"),
                             ("Known variants file: ", "Variants")):
        if item in panel["Filenames"]:
            table.rows.append([descriptor, panel["Filenames"][item]])
    table.rows += [["Minimum Depth: ", str(panel["Depth"])], ["Date of CoverMi analysis: ", time.strftime("%d/%m/%Y")], ["CoverMi version: ", "1.4"]]
    return table.formated()


def location(gr1, panel):
    if "Exons" and "Transcripts" in panel:
        loc = (panel["AllExons"] if ("AllExons" in panel) else panel["Exons"]).overlapped_by(gr1).names_as_string
        if loc =="":
            loc = "{0} intron".format((panel["Alltranscripts"] if ("AllTranscripts" in panel) else panel["Transcripts"]).overlapped_by(gr1).names_as_string)
            if loc == " intron":
                loc = "not covering a "+("" if ("AllTranscripts" in panel) else"targeted ")+"gene"
    else:
        loc = ""
    return loc


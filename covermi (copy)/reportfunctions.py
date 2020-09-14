from __future__ import print_function, absolute_import, division

import time, os, pdb

try:
    python3test = basestring
except NameError:
    basestring = str

from .include import __version__

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
                        newtab[row].append(fstring.format(table[row][col]).ljust(biggest) if isinstance(table[row][col], basestring) else fstring.format(table[row][col]).rjust(biggest))
        return newtab


    def formated(self, sep="", sortedby=None, reverse=False, page_width=99, trim_columns=()):
        if sortedby is not None:
            def bycolumn(line): return [line[sortedby]] + line
            self.rows.sort(key=bycolumn, reverse=reverse)

        rows = type(self)._aligned(self.rows, delete_empty_cols=(True if len(self.headers)==0 else False))
        if len(self.headers) > 0:
            rows = type(self)._aligned(self.headers + rows)

        for column_number, minimum_width in trim_columns:
            excess_width = len(sep.join(rows[0])) - page_width
            if excess_width > 0:
                current_column_width = len(rows[0][column_number])
                if current_column_width > minimum_width:
                    new_column_width = max(current_column_width - excess_width, minimum_width)
                    for row in rows:
                        if row[column_number][new_column_width:].strip() == "":
                            row[column_number] = row[column_number][0:new_column_width]
                        else:
                            row[column_number] = row[column_number][0:new_column_width-3]+"..."

        if len(self.headers) > 0:
            table = [sep.join(row)+"\n" for row in rows]
            table = [("-"*(len(table[0])-1))+"\n"] + table[0:len(self.headers)] + [("-"*(len(table[0])-1))+"\n"] + table[len(self.headers):]
        else:
            table = [sep.join(row)+"\n" for row in rows]
        return table


def header(panel, identifier):
    table = TextTable()
    if identifier is not None:
        table.rows.append(["Sample:", identifier.name])
        table.rows.append(["Run:", identifier.run])
    table.rows.append(["Panel:", panel.name])
    for descriptor, key in (("Manifest: ", "manifest"), ("Design Studio Bedfile:", "designstudio"), ("Known variants file: ", "variants")):
        if key in panel.files:
            table.rows.append([descriptor, os.path.basename(panel.files[key])])
    for descriptor, key in (("Panel type:", "reporttype"), ("Minimum Depth: ", "depth")):
        if key in panel.properties:
            table.rows.append([descriptor, str(panel.properties[key])])
    table.rows += [["Date of CoverMi analysis: ", time.strftime("%d/%m/%Y")], ["CoverMi version: ", __version__]]
    return table.formated()


def location(gr1, panel):
    for exons, transcripts in (("exons", "transcripts"), ("allexons", "alltranscripts")):
        loc = getattr(panel, exons).overlapped_by(gr1).names_as_string
        if not loc:
            loc = getattr(panel, transcripts).overlapped_by(gr1).names_as_string
            if loc:
                loc = "{} intron".format(loc)
        if loc:
            break
    if not loc:
        loc = "intergenic"
    return loc


#!/usr/bin/env python3

import pdb
import sys
import argparse
import os
import csv

from covermi import Panel, gff3, Gr, Entry


def genename_bed(panel, output=None):
    """ 

    Args:
        panel (str): Path of panel to be edited.
        output (str): Path to write renamed bed file to. If not
            provided then defaults to $INPUT_FILENAME_genenames.bed.
        
    Returns:
        None.
    """

    panel = Panel(panel)
    
    for needed in ("targets", "reference"):
        if needed not in panel.paths:
            sys.exit(f"Panel does not contain a {needed} file")

    reference = Gr(gff3(panel.paths["reference"], "transcripts"))
    
    if "names" in panel.paths:
        with open(panel.paths["names"]) as f_in:
            names = {l.split()[0]: l.strip() for l in f_in if l.strip()}
    else:
        names = {}
    
    if output is None:
        output = os.path.splitext(panel.paths["targets"])[0] + "_genenames.bed"
    
    if os.path.exists(output):
        sys.exit(f"Output file {output} already exists")
    
    with open(output, "wt") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        with open(panel.paths["targets"]) as f_in:
            reader = csv.reader(f_in, delimiter="\t")
            for row in reader:
                if row:
                    target = Entry(row[0], int(row[1]) + 1, int(row[2]), ".")
                    touching = reference.touched_by(target)
                    
                    genes = [names[e.name] for e in touching if e.name in names]
                    if not genes:
                        genes = [e.name for e in touching]
                        genes = genes+["******************"] if genes else "."
                    
                    writer.writerow(row[:3] + [";".join(sorted(genes))])



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("panel", help="Directory containing panel data.")
    parser.add_argument("-o", "--output", help="Output bed file path.", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        genename_bed(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()


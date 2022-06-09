#!/usr/bin/env python3

import pdb
import sys
import argparse
import os
import csv

from covermi import Panel, gff3, Gr, Entry


def genename_bed():
    """ 
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", required=False)
    parser.add_argument("-b", "--bed", help="Input bed file.", required=False)
    parser.add_argument("-a", "--annotation", help="GFF3 file containing gene names.", required=False)    
    parser.add_argument("-n", "--names", help="File containing list of required gene names.", required=False)    
    parser.add_argument("-o", "--output", help="Output bed file path.", default="genenames.bed", required=False)
    args = parser.parse_args()
    
    panel = Panel(*[path for path in (args.panel, args.bed, args.annotation, args.names) if path])
    for needed in ("amplicons", "reference"):
        if needed not in panel.paths:
            sys.exit(f"Panel does not contain a {needed} file")
    
    reference = Gr(gff3(panel.paths["reference"], "transcripts"))
    
    wanted = set()
    if "names" in panel.paths:
        with open(panel.paths["names"]) as f_in:
            for row in f_in:
                gene = row.split()
                if gene:
                    wanted.add(gene[0])
    
    if os.path.exists(args.output):
        sys.exit(f"Output file {args.output} already exists")
    
    with open(args.output, "wt") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        for amplicon in panel["amplicons"]:
            touching = [entry.name.split()[0] for entry in reference.touched_by(amplicon)]
            genes = sorted([gene for gene in touching if gene in wanted] or touching)
            writer.writerow([amplicon.chrom, amplicon.start - 1, amplicon.stop, ",".join(genes) or amplicon.name])



def main():
    try:
        genename_bed()
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()


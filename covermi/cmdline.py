import pdb
import argparse
import sys
import os
import json
from collections import Counter, defaultdict
from covermi import Panel, Cov, Plot



def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', "--bam", help="Bam path", required=True)
    parser.add_argument("-p", "--panel", help="Panel path", required=True)
    parser.add_argument("-o", "--output", help="Output path (optiomal). If not supplied then defaults to the current working directory", default=".")
    parser.add_argument("-d", "--depth", help="Minimum required depth (optional)", default=None)
    parser.add_argument("-n", "--name", help="Sample name (optional). If not supplied then defaults to the name of the bam file", default=None)

    args = parser.parse_args()

    if args.name is not None:
        name = args.name

    else:
        name = os.path.splitext(os.path.basename(args.bam))[0]

    if not os.path.isdir(args.output):
        sys.exit("Output directory {args.output} does not exist")
    output = os.path.join(args.output, name)

    depths = [30, 100, 500, 1000, 2000]
    if args.depth is not None:
        depth = int(args.depth)
        depths = sorted(set(depths + [depth]))

    else:
        depth = None

    panel = Panel(args.pane)
    if "targets" in panel:
        roi = panel.targets
    elif "exons" in panel:
        roi = panel.exons
    else:
        sys.exit("Invalid panel. Unable to work out region over which to calculate coverage")

    cov = Cov(args.bam)

    Plot(coverage=cov, panel=panel, depth=depth, title=name, output=output+".pdf")

    stats = {"coverage": {},
             "coverage_by_gene": defaultdict(dict)}
    for depth in depths:
        for i in cov.calculate(roi, depth):
            stats["coverage_by_gene"][f"{depth}x"][i.name] = int(i.percent_covered)
            stats["coverage_by_gene"]["mean_depth"][i.name] = int(i.depth)

        i = cov.calculate(roi, depth, name="Total")
        stats["coverage"][f"{depth}x"] = int(i.percent_covered)
        stats["coverage"]["mean_depth"] = int(i.depth)

    with open(output+".json", "wt") as f_out:
        json.dump(stats, f_out, sort_keys=True, indent=4)




if __name__ == "__main__":
    cmdline()

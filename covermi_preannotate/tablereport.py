import csv, os
import pdb


def create(coverage, panel, identifier, outputstem):

    if "exons" not in panel or "transcripts" not in panel or "depth" not in panel:
        return

    if "amplicons" in panel:
        targeted_range = panel.amplicons.merged
        targeted_exons = panel.exons.overlapped_by(targeted_range)
    else:
        targeted_range = panel.transcripts
        targeted_exons = panel.exons

    # Coverage by Gene        
    filepath, sample = os.path.split(outputstem)
    filepath = os.path.join(filepath, "all_samples_coverage_by_gene.tsv")
    header = ["Sample", "Gene", "Coverage of Targeted Area", "Coverage of Whole Gene", "Mean Depth"] if not os.path.exists(filepath) else []
    with open(filepath, "ab") as f:
        writer = csv.writer(f, delimiter="\t")
        if header:
            writer.writerow(header)

        for i in coverage.calculate(targeted_exons, panel.depth):
            writer.writerow([identifier.name, i.name, i.percent_covered, float(i.bases_covered*100)/panel.exons.subset(i.name).base_count, i.depth_covered])


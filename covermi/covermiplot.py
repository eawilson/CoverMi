from __future__ import print_function, absolute_import, division

try:
    from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle
    from matplotlib.backends.backend_pdf import PdfPages
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False

import pdb
from math import log10
from itertools import cycle

try: # python2
    from itertools import izip as zip
except ImportError: # python3
    pass
    
from .gr import Gr



class Plot(object):
    def __init__(self, cumulative=False):
        self.cumulative = cumulative
        self.data = []
        self.new()

    def new(self, *args):
        self.x = []
        self.y = []
        if args: self.add(*args)

    def add(self, x, y):
        self.x += [self.scalex(x)]
        self.y += [self.scaley(y)]

    def scalex(self, x):
        for beginning, end, fixed, scaling in self.data:
            if x >= beginning and x <= end:
                break
        rel_pos = x - beginning
        return x - (fixed + (rel_pos - (rel_pos//scaling)))

    def scaley(self, y):
        return (log10(y) if y > 1 else 0.0) if not self.cumulative else y / 25.0


def plot(coverage=None, panel=None, identifier=None, fn=None, transcripts=None, cumulative=False, histogram=False):
    if "transcripts" not in panel or not MATPLOTLIB:
        return

    with PdfPages(fn) as pdf:
        transcripts = panel.transcripts.subset(lambda e: e.name in transcripts) if transcripts else panel.transcripts
        for entry in transcripts.sorted_by_name:
            figure = Figure(figsize=(11.7, 4.15))
            FigureCanvasPdf(figure)
            ax = figure.gca()

            transcript = Gr().add(entry)
            amplicons = panel.amplicons.touched_by(transcript) if ("amplicons" in panel) else Gr()
            exons = panel.exons.touched_by(transcript).subset(lambda e: e.name==entry.name) if ("exons" in panel) else Gr()
            codingexons = panel.codingexons.touched_by(transcript).subset(lambda e: e.name==entry.name) if ("codingexons" in panel) else Gr()
            plot_area = transcript.combined_with(amplicons).merged
            plot_start = plot_area.data[entry.chrom][0].start
            plot_stop = plot_area.data[entry.chrom][0].stop

            max_intron_size = amplicons.base_count // amplicons.components if not amplicons.is_empty else 200
            blocks = amplicons.combined_with(exons).merged
            # data = [start, stop, fixed subtraction, scaling factor]
            extra = {"cumulative": True} if cumulative else {}
            plot = Plot(**extra)
            end_of_prev_block = plot_start-2
            fixed_total = 0
            for block in blocks:
                intron_size = block.start - end_of_prev_block - 1
                plot.data += [(end_of_prev_block+1, block.start-1, fixed_total, max(float(intron_size)//max_intron_size, 1))]
                fixed_total += max(intron_size - max_intron_size, 0)
                plot.data += [(block.start, block.stop, fixed_total, 1)]
                end_of_prev_block = block.stop
            plot.data += [(end_of_prev_block+1, plot_stop+1, fixed_total, max(float(plot_stop - end_of_prev_block)//max_intron_size, 1) if (not blocks.is_empty) else 1)]

            # coverage
            if coverage:
                plot.new(plot_start-1, 0)
                for cstart, cstop, cdepth in coverage.data[entry.chrom]:
                    if cstop >= plot_start:
                        plot.add(max(plot_start, cstart), cdepth)
                        if cstop >= plot_stop:
                            plot.add(plot_stop, cdepth)
                            break
                plot.add(plot_stop+1, 0)
                ax.plot(plot.x, plot.y, "-", drawstyle="steps-post", color="dodgerblue", linewidth=1)

            # amplicons
            for amplicon, y in zip(amplicons, cycle((-0.04, 0))):
                ax.add_patch(Rectangle((plot.scalex(amplicon.start), y), amplicon.stop-amplicon.start+1, 0.04, edgecolor="black", facecolor="black", zorder=100))

            # exons
            for regions, colour, drawlabels in ((exons, "darkgrey", True), (codingexons, "black", False)):
                xlabels = []
                xlocations = []
                for exon in regions:
                    ax.add_patch(Rectangle((plot.scalex(exon.start), -0.3), exon.stop-exon.start+1, 0.2, edgecolor="black", facecolor=colour))
                    xlabels += [exon.exon]
                    xlocations += [plot.scalex((exon.start+exon.stop)//2)]
                    if drawlabels:
                        ax.axhline(-0.2, color="black", linewidth=2, zorder=0)
                        step = (len(xlocations) // 20) + 1
                        ax.set_xticks(xlocations[::step])
                        ax.set_xticklabels(xlabels[::step], fontsize=8)
                        ax.tick_params(axis="x",length=0)

            # variants       
            if "variants_mutation" in panel:
                if histogram:
                    for exon in codingexons:
                        num_variants = panel.variants_mutation.overlapped_by(exon).weighted_components
                        ax.add_patch(Rectangle((plot.scalex(exon.start), 0), exon.stop-exon.start+1, plot.scaley(num_variants), facecolor="yellow"))
                else:
                    if cumulative:
                        yvar = 50
                    elif "depth" in panel:
                        yvar = panel.depth
                    else:
                        yvar = 0
                    x = [plot.scalex((variant.start+variant.stop)//2) for variant in plot_area.overlapped_by(panel.variants_mutation)]
                    y = [plot.scaley(yvar)] * len(x)
                    ax.plot(x, y, "x", color="black")

            # depth        
            if "depth" in panel and not cumulative:
                ax.axhline(plot.scaley(panel.depth), color="black", linewidth=0.5, linestyle=":")

            border = (plot.scalex(plot_stop) - plot_start + 3) * 6/100
            minx = plot_start - 1 - border
            maxx = plot.scalex(plot_stop) + 1 + border
            if entry.strand == "-":
                minx, maxx = (maxx, minx)

            ax.set_xlim(minx, maxx)
            ax.spines["bottom"].set_visible(False)

            ax.set_ylim(-0.3, 4.5)
            ax.set_yticks([0, 1, 2, 3, 4])
            if cumulative:
                ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=8)
                ax.set_ylabel("Proportion of time covered")       
            else:
                ax.set_yticklabels(["0", "10", "100", "1000", "10,000"], fontsize=8)
                if histogram:
                    ax.set_ylabel("Variants", fontsize=10)
                else:
                    ax.set_ylabel("Read Depth (log scale)", fontsize=10)
            ax.set_xlabel("Exon", fontsize=10)
            ax.set_title(identifier, fontsize=12)
            ax.add_patch(Rectangle((minx, 4.25), maxx-minx, 0.25, edgecolor="black", facecolor="bisque", zorder=100))
            ax.text((minx+maxx)//2, 4.355, entry.name, zorder=101, ha="center", va="center")
            
            pdf.savefig(figure)#, bbox_inches='tight')
            


def coverage_debug_plot(coverage, fn):
    if not MATPLOTLIB:
        return
    
    with PdfPages(fn) as pdf:
        for chrom in coverage.data.keys():
            cov = coverage.data[chrom]
            start = 1
            for stop in range(start, len(cov)):
                if cov[stop][1] - cov[stop][0] > 1000 and cov[stop][2] == 0:
                    figure = Figure(figsize=(11.7, 4.15))
                    FigureCanvasPdf(figure)
                    ax = figure.gca()
                    ax.plot([c[0] for c in cov[start:stop]], [c[2] for c in cov[start:stop]], "-", drawstyle="steps-post", color="dodgerblue", linewidth=1)
                    ax.set_title(chrom, fontsize=12)
                    pdf.savefig(figure)#, bbox_inches='tight')
                    start = stop + 1



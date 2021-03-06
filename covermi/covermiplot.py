from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

import pdb
from math import log10
from itertools import cycle, chain
from collections import namedtuple
    
from .gr import Gr
from .cov import bisect_left


def log(pos):
    return log10(pos) if pos > 1 else 0.0



class Plot(object):
    def __init__(self, coverage=None, panel=None, depth=None, title=None, output="plot.pdf"):
        self._pdf = PdfPages(output)
        self._title = title

        if panel is not None:
            for name in panel.names:
                amplicons = Gr(entry for entry in panel.targets if entry.name == name)
                exons = Gr(entry for entry in panel.exons if entry.name == name)
                codingexons = Gr(entry for entry in panel.codingexons if entry.name == name)
                self.figure(name, coverage, exons, codingexons, amplicons, depth)
            self.close()
    
    
    def __enter__(self):
        return self
    
    
    def __exit__(self, type_, value, traceback):
        self.close()
    
    
    def close(self):
        self._pdf.close()
    
    
    def figure(self, name, coverage, exons=Gr(), codingexons=Gr(), amplicons=Gr(), depth=None):
        figure = Figure(figsize=(11.7, 4.15))
        FigureCanvasPdf(figure)
        ax = figure.gca()
        
        blocks = amplicons.combined_with(exons).merged
        scale = Scale()
        xticklabels = []
        xticks = []
        for chrom in sorted(blocks.keys()):
            scale.initialse(blocks[chrom])
            
            # coverage
            if coverage:
                xy = [(scale.relative_start, log(0))]
                bisect = bisect_left(coverage[chrom], scale.absolute_start) # leftmost coverage.stop >= entry.start
                for cstart, cstop, cdepth in coverage[chrom][bisect:]:
                    xy.append((scale(max(scale.absolute_start, cstart)), log(cdepth)))
                    if cstop >= scale.absolute_stop:
                        xy.append((scale.relative_stop, log(cdepth)))
                        break
                xy.append((scale.relative_stop, log(0)))
                xs, ys = zip(*xy)
                ax.plot(list(xs), list(ys), "-", drawstyle="steps-post", color="dodgerblue", linewidth=1)

            # amplicons
            for amplicon, y in zip(amplicons.get(chrom, ()), cycle((-0.04, 0))):
                ax.add_patch(Rectangle((scale(amplicon.start), y), amplicon.stop-amplicon.start+1, 0.04, edgecolor="black", facecolor="black", zorder=100))

            # exons
            _xticks = []
            for exon in exons.get(chrom, ()):
                ax.add_patch(Rectangle((scale(exon.start), -0.3), exon.stop-exon.start+1, 0.2, edgecolor="black", facecolor="darkgrey"))
                _xticks += [scale((exon.start+exon.stop)//2)]
            for exon in codingexons.get(chrom, ()):
                ax.add_patch(Rectangle((scale(exon.start), -0.3), exon.stop-exon.start+1, 0.2, edgecolor="black", facecolor="black"))
            xticklabels.extend(range(len(_xticks), 0, -1) if exons and next(iter(exons)).strand == "-" else range(1, len(_xticks)+1))
            xticks.extend(_xticks)

            ## variants       
            #if "variants_mutation" in panel:
                #if "depth" in panel:
                    #yvar = panel.depth
                #else:
                    #yvar = 0
                #x = [scale((variant.start+variant.stop)//2) for variant in plot_area.overlapped_by(panel.variants_mutation)]
                #y = [log(yvar)] * len(x)
                #ax.plot(x, y, "x", color="black")

        # depth        
        if depth:
            ax.axhline(log(depth), color="black", linewidth=0.5, linestyle=":")

        border = (scale.relative_stop - scale.relative_start) * 6/100
        minx = scale.relative_start - border
        maxx = scale.relative_stop + border
        if exons and next(iter(exons)).strand == "-":
            minx, maxx = (maxx, minx)

        ax.set_xlim(minx, maxx)
        ax.spines["bottom"].set_visible(False)

        step = (len(xticks) // 20) + 1
        ax.set_xticks(xticks[::step])
        ax.set_xticklabels(xticklabels[::step], fontsize=8)
        ax.tick_params(axis="x",length=0)
        
        ax.set_ylim(-0.3, 4.5)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_yticklabels(["0", "10", "100", "1000", "10,000"], fontsize=8)
        ax.set_ylabel("Read Depth (log scale)", fontsize=10)
        
        ax.axhline(-0.2, color="black", linewidth=2, zorder=0)
        #ax.set_xlabel("Exon", fontsize=10)
        ax.set_title(self._title, fontsize=12)
        ax.add_patch(Rectangle((minx, 4.25), maxx-minx, 0.25, edgecolor="black", facecolor="bisque", zorder=100))
        ax.text((minx+maxx)//2, 4.355, name, zorder=101, ha="center", va="center")
        
        self._pdf.savefig(figure)#, bbox_inches='tight')

    
    
Correction = namedtuple("Correction", ["start", "stop", "fixed", "scaled"])



class Scale():
    def initialse(self, blocks):
        iterblocks = iter(blocks)
        block = next(iterblocks)
        fixed = block.start
        
        try:
            fixed -= self._corrections[-1].stop
            self._corrections.append(Correction(block.start - 401, block.start - 1, None, None))
        except AttributeError:
            self._corrections = []
            
        for next_block in chain(iterblocks, (None,)):
            self._corrections.append(Correction(block.start, block.stop, fixed, 1))
            if next_block is not None:
                actual_intron_size = next_block.start - block.stop - 1
                scaled_intron_size = min(actual_intron_size, 200)
                self._corrections.append(Correction(block.stop + 1, next_block.start - 1, fixed, float(scaled_intron_size) / actual_intron_size))
                if scaled_intron_size < actual_intron_size:
                    fixed += actual_intron_size - scaled_intron_size
                block = next_block
    
    @property
    def absolute_start(self):
        return self._corrections[0].start
    
    @property
    def relative_start(self):
        return self(self.absolute_start)
    
    @property
    def absolute_stop(self):
        return self._corrections[-1].stop
    
    @property
    def relative_stop(self):
        return self(self.absolute_stop)

    def __call__(self, pos):
        for correction in self._corrections:
            if correction.start <= pos <= correction.stop:
                break
        rel_pos = pos - correction.start
        return correction.start - correction.fixed + (rel_pos * correction.scaled)




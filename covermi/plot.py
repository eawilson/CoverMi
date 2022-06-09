from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

import pdb
import re
from math import log10
from itertools import cycle
from collections import namedtuple, defaultdict
    
from .gr import Gr, chromosome_order
from .cov import bisect_left



def log(pos):
    return log10(pos) if pos > 1 else 0.0



def linear(pos):
    return pos / 1000.0



class Plot(object):
    def __init__(self, panel, coverage=None, title=None, output="plot.pdf", scale="log"):
        self.pdf = PdfPages(output)
        
        genes = defaultdict(list)
        for name in panel["names"]:
            splitname = name.split()
            if splitname:
                genes[splitname[0]].extend(splitname[1:])
        
        exons = defaultdict(lambda:defaultdict(Gr))
        for entry in panel["exons"]:
            gene, transcript = entry.name.split()
            exons[gene][transcript].add(entry)
        
        codingexons = defaultdict(lambda:defaultdict(Gr))
        for entry in panel["codingexons"]:
            gene, transcript = entry.name.split()
            codingexons[gene][transcript].add(entry)
        
        for gene, transcripts in sorted(genes.items()):
            amplicons = Gr(entry for entry in panel.get("amplicons", ()) if entry.name == gene)
            variants = Gr(entry for entry in panel["variants"] if entry.name == gene)
            for transcript in sorted(transcripts or exons[gene].keys() or (None,)):
                subtitle = f"{gene} {transcript}" if transcript else gene
                self.pdf.savefig(figure(title, subtitle, coverage, exons[gene][transcript], codingexons[gene][transcript], amplicons, variants, scale))
        self.close()
    
    
    def __enter__(self):
        return self
    
    
    def __exit__(self, type_, value, traceback):
        self.close()
    
    
    def close(self):
        self.pdf.close()
    

    
def figure(title, subtitle, coverage, exons=Gr(), codingexons=Gr(), amplicons=Gr(), variants=Gr(), scale="log"):
    if scale == "log":
        yscale = log
    elif scale == "linear":
        yscale = linear
    else:
        raise ValueError(f'Invalid scale "{scale}", must be "log" or "linear"')

    figure = Figure(figsize=(11.7, 4.15))
    FigureCanvasPdf(figure)
    ax = figure.gca()
    
    allblocks = amplicons.combined_with(exons).merged
    xscale = Scale()
    xticklabels = []
    xticks = []
    for chrom in sorted(allblocks.keys(), key=chromosome_order):
        blocks = allblocks[chrom]
        xscale.initialse(blocks)
        
        # coverage
        if coverage:
            xy = [(xscale(blocks[0].start), yscale(0))]
            bisect = bisect_left(coverage[chrom], blocks[0].start) # leftmost coverage.stop >= entry.start
            for cstart, cstop, cdepth in coverage[chrom][bisect:]:
                xy.append((xscale(max(blocks[0].start, cstart)), yscale(cdepth)))
                if cstop >= blocks[-1].stop:
                    xy.append((xscale(blocks[-1].stop), yscale(cdepth)))
                    break
            xy.append((xscale(blocks[-1].stop), yscale(0)))
            x, y = zip(*xy)
            ax.plot(list(x), list(y), "-", drawstyle="steps-post", color="dodgerblue", linewidth=1)
        
        # variants
        x = []
        ymax = []
        for variant in variants.get(chrom, ()):
            try:
                pos = xscale((variant.start+variant.stop)//2)
                if x[-1] == pos:
                    ymax[-1] += 1
                    continue
            except IndexError:
                pass
            x.append(pos)
            ymax.append(1)
        if x:
            ax.plot(x, [yscale(y * 10) for y in ymax], "x", color="black")

        # amplicons
        for amplicon, y in zip(amplicons.get(chrom, ()), cycle((-0.04, 0))):
            ax.add_patch(Rectangle((xscale(amplicon.start), y), amplicon.stop-amplicon.start+1, 0.04, edgecolor="black", facecolor="black", zorder=100))

        # exons
        _xticks = []
        for exon in exons.get(chrom, ()):
            ax.add_patch(Rectangle((xscale(exon.start), -0.3), exon.stop-exon.start+1, 0.2, edgecolor="black", facecolor="darkgrey"))
            _xticks += [xscale((exon.start+exon.stop)//2)]
        for exon in codingexons.get(chrom, ()):
            ax.add_patch(Rectangle((xscale(exon.start), -0.3), exon.stop-exon.start+1, 0.2, edgecolor="black", facecolor="black"))
        xticklabels.extend(range(len(_xticks), 0, -1) if exons and next(iter(exons)).strand == "-" else range(1, len(_xticks)+1))
        xticks.extend(_xticks)

    border = xscale(blocks[-1].stop) * 6/100
    minx = -border
    maxx = xscale(blocks[-1].stop) + border
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
    if yscale == log:
        ax.set_yticklabels(["0", "10", "100", "1000", "10,000"], fontsize=8)
        ax.set_ylabel("Read Depth (log scale)", fontsize=10)
    else:
        ax.set_yticklabels(["0", "1000", "2000", "3000", "4000"], fontsize=8)
        ax.set_ylabel("Read Depth", fontsize=10)
    
    ax.axhline(-0.2, color="black", linewidth=2, zorder=0)
    #ax.set_xlabel("Exon", fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.add_patch(Rectangle((minx, 4.25), maxx-minx, 0.25, edgecolor="black", facecolor="bisque", zorder=100))
    ax.text((minx+maxx)//2, 4.355, subtitle, zorder=101, ha="center", va="center")
    
    return figure

    
    
Correction = namedtuple("Correction", ["start", "stop", "fixed", "scaled"])



class Scale():
    def initialse(self, blocks):        
        fixed = blocks[0].start 
        try:
            fixed -= self._corrections[-1].stop + 1 - self._corrections[-1].fixed + 400
        except AttributeError:
            pass
        self._corrections = []
            
        for i in range(len(blocks)):
            self._corrections.append(Correction(blocks[i].start, blocks[i].stop, fixed, 1))
            if i < len(blocks) - 1:
                actual_intron_size = blocks[i+1].start - blocks[i].stop - 1
                scaled_intron_size = min(actual_intron_size, 200)
                self._corrections.append(Correction(blocks[i].stop + 1, blocks[i+1].start - 1, fixed, float(scaled_intron_size) / actual_intron_size))
                if scaled_intron_size < actual_intron_size:
                    fixed += actual_intron_size - scaled_intron_size

    def __call__(self, pos):
        for correction in self._corrections:
            if correction.start <= pos <= correction.stop:
                break
        rel_pos = pos - correction.start
        return correction.start - correction.fixed + (rel_pos * correction.scaled)




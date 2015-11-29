from pkg_resources import resource_string
import pdb
from gr import Gr
import subprocess, os, sys, re

r_path = "unknown"

class _Encode(object):
    def __init__(self, name, strand, adj):
        self.name = name
        self.strand = strand
        self.adj = adj

    def encode(self, x, y, component):
        for beginning, end, fixed, scaling in self.adj:
            if x >= beginning and x <= end:
                break
        rel_pos = x - beginning
        adjustment = fixed + (rel_pos - (rel_pos/scaling))
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(x, y, self.name, component, self.strand, adjustment)


def plot(coverage, panel, outputstem):
    
    global r_path
    if r_path == "unknown":
        if os.name ==  "nt":
            newest_time = 0
            r_path = "unavailable"
            R_STEM = "C:\\Program Files\\R"
            for root, dirnames, filenames in os.walk(R_STEM):
                if "Rscript.exe" in filenames:
                    exe_path = "{0}\\Rscript.exe".format(root)
                    mtime = os.path.getmtime(exe_path)
                    if mtime > newest_time:
                        newest_time = mtime
                        r_path = "START /B \""+exe_path+"\" \"{0}\""
        else:
            try:
                with file(os.devnull, "wt") as DEVNULL:
                    subprocess.check_call(["which", "R"], stdout=DEVNULL, stderr=DEVNULL)
                r_path = "R -f '{0}'"
            except CalledProcessError:
                r_path = "unavailable"
        if r_path == "unavailable":
            print "WARNING R binary not found. no graphs will be plotted"

    if "Transcripts" not in panel or r_path == "unavailable":
        return

    output_file = outputstem+"_plot.pdf"
    plotname = panel["Filenames"]["Sample"] if ("Sample" in panel["Filenames"]) else panel["Filenames"]["Panel"]
    rdataframe = outputstem+"_Rdata.tsv"
    with file(rdataframe, "wt") as f:
        s = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n"
        f.write(s.format("x", "y", "name", "component", "strand", "adjustment"))

        for entry in panel["Transcripts"].all_entries:
                chrom = entry.chrom
                name = entry.name
                strand = entry.strand
                gr1 = Gr(entry)

                if "Amplicons" in panel:
                    gr1, amplicons = gr1.extended_to_include_touching_amplicons(panel["Amplicons"])
                else:
                    amplicons = Gr()
                auc = 0
                for amplicon in amplicons.all_entries:
                    auc += amplicon.stop - amplicon.start + 1                
                INTRON_SIZE = auc / amplicons.number_of_components if (amplicons.number_of_components>0) else 200

                start = gr1[chrom][0].start
                stop = gr1[chrom][0].stop
                
                exons = panel["Exons"].touched_by(gr1).subset2(name) if ("Exons" in panel) else Gr()

                blocks = amplicons.combined_with(exons).merged
                # adj = [start, stop, fixed subtraction, scaling factor]
                adj = []
                end_of_prev_block = start-2
                fixed_total = 0
                for bentry in blocks.all_entries:
                    bstart = bentry.start
                    bstop = bentry.stop
                    intron_size = bstart - end_of_prev_block - 1
                    adj.append((end_of_prev_block+1, bstart-1, fixed_total, max(float(intron_size)/INTRON_SIZE, 1)))
                    fixed_total += max(intron_size - INTRON_SIZE, 0)
                    adj.append((bstart, bstop, fixed_total, 1))
                    end_of_prev_block = bstop
                adj.append((end_of_prev_block+1, stop+1, fixed_total, max(float(stop - end_of_prev_block)/INTRON_SIZE, 1) if (not blocks.is_empty) else 1))

                line = _Encode(name, strand, adj)

                if "Amplicons" in panel:
                    for aentry in amplicons.all_entries:
                        f.write(line.encode(aentry.start, aentry.stop, "amplicon"))
                   
                f.write(line.encode(start-1, 0, "coverage"))
                for cstart, cstop, cdepth in coverage[chrom]:
                    if cstop >= start:
                        f.write(line.encode(max(start, cstart), cdepth, "coverage"))
                        if cstop >= stop:

                            f.write(line.encode(min(stop, cstop), cdepth, "coverage"))
                            break
                f.write(line.encode(stop+1, 0, "coverage"))

                if "Exons" in panel:
                    for eentry in exons.all_entries:
                        f.write(line.encode(eentry.start, eentry.stop, "exon"))
                        f.write(line.encode((eentry.start+eentry.stop)/2, str(eentry.exon), "exon_number"))
        
                if "Variants_Mutation" in panel:
                    variants = gr1.overlapped_by(panel["Variants_Mutation"])
                    vycoord = panel["Options"]["Depth"] if ("Depth" in panel["Options"]) else 0
                    for ventry in variants.all_entries:
                            f.write(line.encode((ventry.start+ventry.stop)/2, vycoord, "variants"))
        
                if "Depth" in panel["Options"]:
                    f.write(line.encode(start, panel["Options"]["Depth"], "minimum"))

                if "Other" in panel:
                    for oentry in gr1.overlapped_by(panel["Other"]).all_entries:
                        f.write(line.encode(oentry.start, oentry.stop, "other"))

    
    genericrcode = resource_string(__name__, "covermiplot.R")
    rscript = outputstem+"R"
    with file(rscript, "wt") as f:
        f.write(genericrcode.replace("INPUT", os.path.abspath(rdataframe)).replace("OUTPUT", os.path.abspath(output_file)).replace("NAME", plotname).replace("\\", "/"))
    
    with file(os.devnull, "wt") as DEVNULL:
        subprocess.call(r_path.format(rscript), shell=True, stdout=DEVNULL)
    
    os.unlink(rdataframe)
    os.unlink(rscript)



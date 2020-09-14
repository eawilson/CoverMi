from __future__ import print_function, absolute_import, division
import sys, os, tkFileDialog, Tkinter, pdb, csv, openpyxl, subprocess, socket

from .gr import vcf
from .panel import Panel
from .annotate import Filter, vep


class Igv(object):

    def __init__(self):
#    igv = subprocess.Popen(["/home/ed/IGV_2.4.5/igv.sh"], stdout=subprocess.PIPE, shell=True)
        try:
            self.sock = socket.create_connection(("127.0.0.1", 60151))
        except socket.error:    
            print("Unable to connect to igv")
            self.sock = False
        self.new()
        self.command("maxPanelHeight 1000000")

    def command(self, command):
        if self.sock:
            try:
                self.sock.sendall("{}\n".format(command))
                return self.sock.recv(4096) == "OK\n"
            except socket.error:
                self.sock = self.bamloaded = False
        return False

    def load(self, path):
        self.bamloaded = self.command("load {}".format(path))

    def new(self):
        self.command("new")
        self.bamloaded = False


def main():
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print("Please select a panel")
    panelpath = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a panel')
    if not bool(panelpath):
        sys.exit()
    panelpath = os.path.abspath(panelpath)
    print("{0} panel selected".format(os.path.basename(panelpath)))

    print("Please select the folder containing the vcf files")
    rootdir = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a folder of vcf files')
    if not bool(rootdir):
        sys.exit()
    rootdir = os.path.abspath(rootdir)
    print("{0} selected".format(os.path.basename(rootdir)))

    print("Please select location of output folder")
    outputpath = tkFileDialog.askdirectory(parent=rootwindow, initialdir=rootdir, title='Please select location of output folder')
    if not bool(outputpath):
        sys.exit()
    outputpath = os.path.abspath(os.path.join(outputpath, "{}_annotations".format(os.path.basename(rootdir))))
    try:
        os.mkdir(outputpath)
    except OSError:
        raise RuntimeError("Error {} already exists".format(outputpath) if os.path.exists(outputpath) else "Error creating {}".format(outputpath))
    print("{0} selected".format(os.path.basename(outputpath)))

    panel = Panel(panelpath, verbose=True)

    igv = Igv()

    for fn in os.listdir(rootdir):
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "variants"
        identifier, extension = os.path.splitext(fn)
        vcfpath = os.path.join(rootdir, fn)
        if extension == ".vcf" and os.path.isfile(vcfpath):
            print(identifier)
            bampath = os.path.join(rootdir, "{}.bam".format(identifier))
            if os.path.exists(bampath) and os.path.exists("{}.bai".format(bampath)):
                igv.load(bampath)

            ws.append(["Gene", "Chrom", "Coordinate", "Change", "VAF", "Qual", "Filters", "Depth", "Filtered", "Transcript", "HGVSc", "HGVSp", "Impact", 
                       "Consequence", "Sift", "Polyphen", "MAF", "dbSNP", "PUBMED", "COSMIC", "IGV"])
            variants = list(vcf(vcfpath))
            annotations = vep(variants, panel=panel)
            Filter(panel.properties.get("filters", "")).annotate(annotations)
            imagelist = set()
            for index, v in enumerate(annotations, start=2):

                if v["passedfilters"]:
                    if igv.bamloaded:
                        centre = (v["start"]+v["stop"])//2
                        location = "{}:{}".format(v["chrom"], centre)
                        jpgpath = os.path.join(outputpath, "{}_{}.jpg".format(identifier, location.replace(":", "_")))
                        if jpgpath not in imagelist:
                            igv.command("goto {}".format(location))
                            igv.command("snapshot {}".format(jpgpath))
                            imagelist.add(jpgpath)
                        img = openpyxl.drawing.image.Image(jpgpath)
                        imgsheet = wb.create_sheet(str(index))
                        img.anchor(imgsheet["A1"])
                        imgsheet.add_image(img)
                        hyperlink = '=HYPERLINK("#{}!A1", "IGV")'.format(index)
                    else:
                        hyperlink = "=HYPERLINK(\"http://localhost:60151/load?file={}&locus={}:{}\")".format(bampath, v["chrom"], v["pos"])
                else:
                    hyperlink = ""

                ws.append([v.get("gene_symbol", ""),
                           str(v["chrom"]),
                           v["pos"],
                           "{}/{}".format(v["ref"], v["alt"]),
                           int(v["vaf"] * 1000) / 10.0,
                           v["qual"],
                           v["filters"],
                           v["depth"],
                           "PASS" if v["passedfilters"] else "FAIL",
                           v.get("transcript_id", ""),
                           v.get("hgvsc", ""),
                           v.get("hgvsp", ""),
                           v.get("impact"),
                           ", ".join(v.get("consequence_terms", ())),
                           v.get("sift_prediction", ""),
                           v.get("polyphen_prediction", ""),
                           v.get("maf", ""),
                           ", ".join(v.get("dbsnp", ())),
                           ", ".join(v.get("pubmed", ())),
                           ", ".join(v.get("cosmic", ())),
                           hyperlink,
                          ])

            wb.save(os.path.join(outputpath, "{}.annotation.xlsx".format(identifier)))
            for jpgpath in imagelist:
                os.unlink(jpgpath)
            igv.new()








if __name__ == "__main__":
    main()

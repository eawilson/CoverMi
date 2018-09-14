import sys, os, tkFileDialog, Tkinter, pdb, csv
from collections import defaultdict
from itertools import izip

def main():
    rootwindow = Tkinter.Tk()
    rootwindow.withdraw()

    print("Please select directory containing knownCanonical, knownGene, knownToRefSeq and RefFlat files downloaded from UCSC")
    path = tkFileDialog.askdirectory(parent=rootwindow, title='Please select a folder')
    if not bool(path):
        sys.exit()

    known_2_exons = {}
    with open(os.path.join(path, "knownGene.txt"), "rb") as f:
        for row in csv.reader(f, delimiter="\t"):
            known_2_exons[row[0]] = row[8]+row[9] # exon start, stop positions

    known_2_refseq = {}
    with open(os.path.join(path, "knownToRefSeq.txt"), "rb") as f:
        for row in csv.reader(f, delimiter="\t"):
            known_2_refseq[row[0]] = row[1]

    canonicals = set()
    canonical_exons = set()
    with open(os.path.join(path, "knownCanonical.txt"), "rb") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[4] in known_2_refseq:
                canonicals.add(known_2_refseq[row[4]])
                canonical_exons.add(known_2_exons[row[4]])

    gene_2_transcripts = defaultdict(set)
    ref_2_exons = {}
    transcript_length = {}
    with open(os.path.join(path, "refFlat.txt"), "rb") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[1] in canonicals:
                gene_2_transcripts[row[0]].add(row[1])
                ref_2_exons[row[1]] = row[9]+row[10]
                cstart = int(row[6]) + 1
                cstop = int(row[7])
                length = 0
                for estart, estop in izip(row[9].strip(",").split(","), row[10].strip(",").split(",")):
                    exonlength = min(int(estop), cstop) - max(int(estart) + 1, cstart) + 1
                    if exonlength > 0:
                        length += exonlength
                transcript_length[row[1]] = length

    opath = "canonical.txt"
    if os.path.exists(opath):
        raise("File {} already exists".format(opath))

    with open("canonical.txt", "wb") as f:
        writer = csv.writer(f, delimiter="\t")
        for gene, alltranscripts in gene_2_transcripts.items():
            transcripts = [transcript for transcript in alltranscripts if transcript in canonicals]
            if len(transcripts) > 1:
                transcripts = [transcript for transcript in transcripts if ref_2_exons[transcript] in canonical_exons]
            if len(transcripts) > 1:
                transcripts = sorted(transcripts, key=lambda x: transcript_length[x], reverse=True)
                if transcript_length[transcripts[0]] > transcript_length[transcripts[1]]:
                    transcripts = transcripts[:1]
            if len(transcripts) == 1:
                writer.writerow((gene, transcripts[0]))

    print "Written file {}".format(opath)
    print "This file needs to be placed in all panel directories"


if __name__ == "__main__":
    main()

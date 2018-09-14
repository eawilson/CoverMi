from covermi2 import annotate
from cStringIO import StringIO
from collections import namedtuple
import pdb
from collections import defaultdict
import csv

features = {"ACT_SITE": "Active site",
            "METAL": "Metal binding",
            "BINDING": "Binding site",
            "SITE": "Site",
            "CA_BIND": "Calcium binding",
            "ZN_FING": "Zinc finger",
            "DNA_BIND": "DNA binding",
            "NP_BIND": "Nucleotide binding",
            "DOMAIN": "Domain",
            "REGION": "Region",
            "MOTIF": "Motif",
           }

#COILED
#MUTAGEN
#HELIX
#PROPEP
#VAR_SEQ
#LIPID
#CHAIN
#TURN
#CARBOHYD
#CONFLICT
#REPEAT
#TRANSMEM
#DISULFID
#CROSSLNK
#TOPO_DOM
#STRAND
#MOD_RES
#SIGNAL
#PEPTIDE

class NullList(object):
    def __iadd__(self, other):
        return self


genes = []
for genefile in ("/home/ed/Desktop/icbd/ipd/IPD v1/ipd_gene_list.txt", "/home/ed/Desktop/icbd/icbd/ICBD v1/guess_gene_list.txt"):
    with open(genefile, "rt") as f:
        for gene in f:
            gene = gene.split()[0].strip()
            if gene: genes += [gene]

nulllist = NullList()
domains = defaultdict(list)
for gene, text in annotate.uniprottext(genes, cachepath="cache"):
    print gene
    for row in StringIO(text):
        if row[:2] == "FT":
            info = [row[33:].strip()]
            try:
                feature, start, stop = row[2:33].split()
            except ValueError:
                lastinfo += info
                continue
            if feature in features:
                domains[gene] += [[features[feature], int(start), int(stop), info]]
                lastinfo = info
            else:
                lastinfo = nulllist


with open("icbd_ipd_domains.tsv", "wb") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Gene", "Feature", "AA Start", "AA Stop", "Description"])
    for gene, domains in domains.items():
        for domain in domains:
            domain[3] = " ".join(domain[3]).split("{")[0].strip(". ")
            writer.writerow([gene] + domain)














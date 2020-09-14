from covermi.gr import reference, Gr
import unittest


refflat = """HBB	NM_000518	chr11	-	5246695	5248301	5246827	5248251	3	5246695,5247806,5248159,	5246956,5248029,5248301,""".split("\n")

ensembl = """#!genome-build GRCh37.p13
#!genome-version GRCh37
11	ensembl_havana	gene	5246694	5250625	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
11	ensembl_havana	transcript	5246694	5248301	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; tag "basic";
11	ensembl_havana	exon	5248160	5248301	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "1"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; exon_id "ENSE00001829867"; exon_version "2"; tag "basic";
11	ensembl_havana	CDS	5248160	5248251	.	-	0	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "1"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; protein_id "ENSP00000333994"; protein_version "3"; tag "basic";
11	ensembl_havana	start_codon	5248249	5248251	.	-	0	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "1"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; tag "basic";
11	ensembl_havana	exon	5247807	5248029	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "2"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; exon_id "ENSE00001057381"; exon_version "1"; tag "basic";
11	ensembl_havana	CDS	5247807	5248029	.	-	1	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "2"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; protein_id "ENSP00000333994"; protein_version "3"; tag "basic";
11	ensembl_havana	exon	5246694	5246956	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "3"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; exon_id "ENSE00001600613"; exon_version "2"; tag "basic";
11	ensembl_havana	CDS	5246831	5246956	.	-	0	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "3"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; protein_id "ENSP00000333994"; protein_version "3"; tag "basic";
11	ensembl_havana	stop_codon	5246828	5246830	.	-	0	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; exon_number "3"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; tag "basic";
11	ensembl_havana	five_prime_utr	5248252	5248301	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; tag "basic";
11	ensembl_havana	three_prime_utr	5246694	5246827	.	-	.	gene_id "ENSG00000244734"; gene_version "2"; transcript_id "ENST00000335295"; transcript_version "4"; gene_name "HBB"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "HBB-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS7753"; havana_transcript "OTTHUMT00000142977"; havana_transcript_version "2"; tag "basic";""".split("\n")


class TestReference(unittest.TestCase):

    def test_everything(self):
        gr = Gr(reference(refflat, "exons", targets=["HBB"]))
        for entry in gr:
            print repr(entry)

        gr = Gr(reference(ensembl, "exons", targets=["HBB"]))
        for entry in gr:
            print repr(entry)













if __name__ == '__main__':
    unittest.main()

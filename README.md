## Overview

CoverMi provides coverage checking for next generation sequencing data. It has been designed for use with NGS panels run on the Illumina MiSeq platform although it has also been successfully used with whole genome data produced by the Illumina Hiseq. It will produce a report detailing coverage by gene, exon and known variant (if a list of known variants is provided) and will graph the results. 

## Requirements

*Python >= 2.7.10 or Python >= 3.5.2

*Matplotlib package (optional, used for graphing)

*Biopython package (optional, used for random access to indexed BAM files)

*Requests package (used to access NCBI AND UniProt databases)

*VEP http://www.ensembl.org/vep (PERL scripts used for variant annotation)

## Installation

pip install covermi or run python setup.py install from within the covermi root directory

## Documentation

In order to perform the necessary analysis CoverMi needs access to a collection of files that detail the NGS panel being analysed. All the files that make up each panel are placed inside a directory that is named after the panel in question. The names of the files are unimportant as their identity is determined from the structure of the file contents. The following list details the different files that make up a panel. Other than for the reference genome all files are optional, however the more detailed the panel information provided the more detailed the report will be. Both GRCh37 and GRCh38 human genomes are supported, as are Ensembl and RefSeq transcripts.

*refFlat.txt file detailing the RefSeq transcripts for the human genome (GRCh37 or GRCh38), downloadable from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html) or
gtf file detailing the Ensembl transcripts for the human genome (GRCh37 or GRCh38), downloadable from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html)

*Illumina Manifest file, Illumina DesignStudio bedfile or other bedfile detailing targeted regions of genome.

*List of genes/transcripts over which coverage is to be measured. This is a text file with each line containing either gene name or gene name and transcript name separated by a space. If the transcript is not specified then a transcript will be selected based on rank (N prefixed for RefSeq or merged Ensembl/Havana for Ensembl, principal isoform as specified by APPRIS (if file present in panel directory), longest coding region, longest transcript, lowest accession id).

*List of known variants of interest in the form of a [COSMIC](https://cancer.sanger.ac.uk/cosmic) or [HGMD](http://www.hgmd.cf.ac.uk) file.

*Disease names file, this is automatically created when covermi is first run on a panel. It is of the format disease name as specified in variants file, equals, name to be used in covermi reports. It should be hand edited if needed, if the line ends in an equals with no specified name afterwards then all variants associated with that disease will be ignored.

*Options file, this is automatically created when covermi is first run on a panel. It can be hand edited but this should be rarely needed unless filters are needed for variant annotation.

Covermi is started from the command line with the command covermi and asks the user to select the panel to be used and  either a single bam file or a directory of bam files to be analysed. The name of the bam file (minus anything following a trailing underscore) is taken as the sample name and the containing directory is taken as the run name.

## License

MIT, see LICENSE.txt for further details.

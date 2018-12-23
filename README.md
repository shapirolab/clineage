# Short Tandem Repeat stutter model inferred from direct measurement of in vitro stutter noise

#### Author: Ofir Raz <ofir.raz@weizmann.ac.il> <br> License: GNU v2

[Introduction](#introduction)
[Installation](#installation)
[Reference Genome Configuration](#reference)
[Quick Start](#quick-start)
[Help](#help)
[Citation](#citation)

## Introduction
This is a standalone fork of the Cell Lineage project of Ehud Shapiro's lab meant to enable independent usage of the STR genotyping tools developed as part of this project.

Very short (mono and di-repeats) STRs are particularly difficult to genotype. This tool is meant to aid their genotyping in conditions of severe PCR stutter distortions.


## Installation

    sudo apt-get install git
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda create --name cl python=3.5
    source activate cl

    git clone -b standalone https://github.com/ofirr/clineage.git
    cd clineage
    conda install -c bioconda --yes --file conda_requirements.txt
    pip install -r requirements.txt
    pip install -r requirements-unpackaged.txt

    ./manage.py migrate


## Reference Genome Configuration
Configure CHROMOSOMES_PATH in settings or local_settings to your reference genome folder.

Example for human hg19 reference:
In the settings file, CHROMOSOMES_PATH is set by default to the "reference_genome" subdirectory.
We create this subdirectory and within it, folders for the Taxa and for the reference version:

    mkdir -p reference_genomes/Human/hg19
    cd reference_genomes/Human/hg19

We download a given chromosome from ucsc and extract it:

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz
    tar xvzf chrY.fa.gz

We strip the extracted chromosome file of newlines and its header:

	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < chrY.fa > chrY.tmp
	tail -n +3 chrY.tmp > chrY.txt

Cleanup and step back to the main clineage directory:

	rm chrY.tmp
	cd ../../..


## Quick Start
```
python mutation_calling.py
         --fastqr1       input_R1.fastq
         --fastqr2       input_R2.fastq
         --regions       str_regions.bed
         --outputfile    output.tab
```

## Help
Please contact us at ofir.raz@weizmann.ac.il



## Citation
If you found this tool useful, we would appreciate it if you could cite our manuscript:
**[Short Tandem Repeat stutter model inferred from direct measurement of in vitro stutter noise](https://doi.org/10.1101/065110)**
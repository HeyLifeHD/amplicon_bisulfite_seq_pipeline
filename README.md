# amplicon_bisulfite_seq_pipeline
Pipeline for processing amplicon bisulfite sequencing results.

Getting Started
===============

## Depedencies
- [miniconda](https://conda.io/miniconda.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [cutadapt](https://cutadapt.readthedocs.io/)
- [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [multiqc](https://multiqc.info/)

## Installation
Clone the pipeline by issuing:
```bash
git clone --recursive https://github.com/HeyLifeHD/amplicon_bisulfite_seq_pipeline/
```

## Input
Change the directories of your working directory, adapter sequences and bismark reference genome in the config file. Put your fastq.gz files from an amplicon sequencing experiment into a a folder called "01_fastq_raw/" inside your specified working directory. In addition specify the location of the config file in the Snakefile. Now you are ready to go!

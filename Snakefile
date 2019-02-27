#Snakmake Pipeline


import os
from os import path

SNAKEDIR = path.dirname(workflow.snakefile)
configfile: path.join(SNAKEDIR, "config.yml")
  
workdir: config["workdir_top"]

import sys

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")


import glob
import re

SAMPLES, = glob_wildcards("01_fastq_raw/{sam}_R1.fastq.gz")

NB_SAMPLES = len(SAMPLES)

message(str(NB_SAMPLES))


for sample in SAMPLES:
  message("Sample " + str(sample) + " will be processed")

rule all:
    input:
        #expand("04_methyl/CpG_context_{sample}_R1_val_1_bismark_bt2_pe.txt", sample=SAMPLES)
        #"pipeline_bisSeq.html"
        expand("04_methyl/{sample}_R1_val_1_bismark_bt2_pe.bismark.cov", sample=SAMPLES)

rule multiQC:
    message:
        "Perform MultiQC"
    shell:"""
        multiqc --title "Raw FastQC reports" --filename Raw_FastQC_report --outdir MultiQC/ ./
        """ 

rule trim:
    message:
        "Preprocessing of fastq files"
    input:
        fwd = "01_fastq_raw/{sample}_R1.fastq.gz",
        rev = "01_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    shell:"""
        trim_galore --paired --nextera  {input.fwd} {input.rev} -o 02_fastq_trimmed/
        touch {output.fwd}
        touch {output.rev}
        """

rule align:
    message:
        "Align trimmed fastq files"
    input:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    output:
        bam = "03_bismark/{sample}_R1_val_1_bismark_bt2_pe.bam"
    params: 
        ref= config["Bismark_reference"]
    shell:"""
        bismark -o 03_bismark/ --bam --ambiguous --non_directional \
        {params.ref} -1 {input.fwd} -2 {input.rev}
        touch {output.bam}
        """

rule methyl_extract:
    message:
        "Extract methylation"
    input:
        bam = "03_bismark/{sample}_R1_val_1_bismark_bt2_pe.bam"
    output:
        sam = "04_methyl/{sample}_R1_val_1_bismark_bt2_pe.bismark.cov"
        
    shell:"""
        bismark_methylation_extractor -p --comprehensive \
        --merge_non_CpG --no_overlap -o 04_methyl/ {input.bam} --bedGraph
        gunzip 04_methyl/{sample}_R1_val_1_bismark_bt2_pe.bismark.cov.gz
        touch {output.sam}
        """

rule r_qc:
    message:
        "Perform QC"
    input:
        cov = "04_methyl/",
        bed = config["bed"],
        meta = config["meta"]
    output:
        html = "pipeline_bisSeq.html"
    shell:"""
        Rscript -e "rmarkdown::render('pipeline_bisSeq.Rmd')" {input.cov} {input.bed} 10 {input.meta} "~/"
        touch {output.html}
        """

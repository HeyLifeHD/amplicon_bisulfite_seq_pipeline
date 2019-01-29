#Snakmake Pipeline


import os
from os import path

configfile: "/home/epicwl/c010-datasets/Internal/amplicon_seq_pipeline/Pipeline/config.yml"
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
        expand("04_methyl/CpG_context_{sample}_bismark_bt2_pe.txt", sample=SAMPLES)

rule multiQC:
    message:
        "Perform MultiQC"
    shell:"""
        multiqc --title "Raw FastQC reports" --filename Raw_FastQC_report --outdir MultiQC/ ./
        """ 

rule cutadapt:
    message:
        "Run cutadapt"
    input:
        adapters = config["adapters"],
        fwd = "01_fastq_raw/{sample}_R1.fastq.gz",
        rev = "01_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd = "02_fastq_trimmed/{sample}_R1_trimmed.fastq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_trimmed.fastq.gz"
    shell:"""
        cutadapt -g file:{input.adapters} -q 20,20 -o {output.fwd} {input.fwd}
        cutadapt -g file:{input.adapters} -q 20,20 -o {output.rev} {input.rev}
        """

rule align:
#maybe path to bowtie still has to be added
    message:
        "Align trimmed fastq files"
    input:
        fwd = "02_fastq_trimmed/{sample}_R1_trimmed.fastq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        bam =  "03_bismark/{sample}_R1_trimmed_bismark_bt2_pe.bam"
    params: 
        ref= config["Bismark_reference"]
    shell:"""
        bismark -N 0 -q -o 03_bismark/ --bam --ambiguous --non_directional \
        {params.ref} -1 {input.fwd} -2 {input.rev}
        touch {output.bam}
        """

rule methyl_extract:
    message:
        "Extract methylation"
    input:
        bam = "03_bismark/{sample}_R1_trimmed_bismark_bt2_pe.bam"
    output:
        sam = "04_methyl/CpG_context_{sample}_bismark_bt2_pe.txt"
        
    shell:"""
        bismark_methylation_extractor -p --comprehensive \
        --merge_non_CpG --no_overlap -o 04_methyl/ {input.bam} --bedGraph
        touch {output.sam}
        """

#rule make_bedgraph:
#    message:
#        "Make bedgraphs"
#    input:
#        sam = "04_methyl/CpG_context_{sample}_bismark_bt2_pe.txt"
#    output:
#        bedgraph = "05_bedgraph/{sample}.bedgraph"
#    shell:"""
#        bismark2bedGraph --dir 05_bedgraph/ -o {output.bedgraph} {input.sam}
#        """



#rule all: ## run the whole pipeline
#    input:
#        cutadapt = rules.cutadapt.output.fwd,
#        align = rules.align.output.bam,
#        methyl_extract = rules.methyl_extract.output.sam,
#        make_bedgraph = rules.make_bedgraph.output.bedgraph,


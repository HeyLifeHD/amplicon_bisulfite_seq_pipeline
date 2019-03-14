# amplicon_bisulfite_seq_pipeline
Pipeline for processing of amplicon bisulfite sequencing results.


## Authors
* Joschka Hey (@HeyLifeHD)
* Maximilian Schönung (@MaxSchoenung)

## Usage

### Step 1: Install workflow

To use this pipeline simply download and extract the [latest release](https://github.com/HeyLifeHD/amplicon_bisulfite_seq_pipeline).
```bash
git clone --recursive https://github.com/HeyLifeHD/amplicon_bisulfite_seq_pipeline/
```
For further modifications of this pipeline, fork this repository and consider providing applicable modifications via pull request. 
In case you use this workflow in a publication, dont forget to give credits to the authors by citing the URL of this repository.

### Step 2: Configure workflow

Configure the pipeline by editing the `config.yaml` file. Specify the `workdir_top`, `Bismark_reference` and number of `threads`to use. To run the report you should also specify a .bed file containing the regions of interest as well as a .txt file containing the sample names. 
Last create an folder called `01_fastq_raw` in the working directory in which you move the .fastq files to be processed. 

### Step 3: Execute workflow

All you need to execute this pipeline is to install Snakemake via the [Conda package manager](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda). Software needed by this workflow is automatically deployed into isolated environments by Snakemake.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda all

Execute the report of the analysis via

    snakemake --use-conda r_qc


## Depedencies
Dependencies are automatically deployed into isolated environments by Snakmake in case you execute the pipeline with --use-conda.

- [miniconda](https://conda.io/miniconda.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [cutadapt](https://cutadapt.readthedocs.io/)
- [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [multiqc](https://multiqc.info/)

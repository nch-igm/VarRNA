# Variant calling and classification with RNA-Seq data. 



Setup
=======

## Installation

A conda environment is used to manage *almost* all packages. Please refer to the following documentation to install micromamba (or conda), if not already installed: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html.

Create the micromamba environment and activate it.
```bash
micromamba env create -f dependencies/mamba_environment.yml
micromamba activate varrna
```

## Dependencies

### ANNOVAR

ANNOVAR is used to add annotations to our VCFs using several databases. You will need to install this yourself since it requires a user agreement license. Please download annovar from their [download page](https://annovar.openbioinformatics.org/en/latest/user-guide/download/), and use the link provided from registering to download to the following directory:

```bash
wget <link/to/annovar/tar/file> -P dependencies/
tar -xvzf dependencies/annovar.latest.tar.gz -C dependencies/
```

All the databases that annovar uses will be downloaded in the ```get_resourches.sh``` script, below.

### Salmon

[Salmon](https://github.com/COMBINE-lab/salmon) read quantification is used to get TPM values. We'll need to install this ourselves, create the index, and then point to it's installation location in our config file.

```bash
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz -P dependencies/
tar -xvzf dependencies/salmon-1.9.0_linux_x86_64.tar.gz -C dependencies/
```

Note, config/config.yaml references this installation to run salmon. If you download a different version or copy it to a different location, you will need to update that line in the ```config/config.yaml``` file (under dependencies).



## Resources

A number of resources are needed for variant calling and annotation. Download these by running:
```bash
bash get_resources.sh
```

# Usage

## Input data

- Input RNA BAM files

Input BAM files are expected to be processed (e.g., run through fastqc, reads are trimmed and aligned).

- Input RNA FASTQ files

If you are starting with FASTQs you can use the optional Alignment scripts first before running the pipeline to get your BAM files. See Alignment folder. 

Steps:
``` bash
# Download STAR
bash get_star.sh 
# Build the genome. The sjdbOverhang is dependent on sequencing library and should be your <max(read length) - 1>.
bash star_genome_build.sh
# Alignment using the STAR 2-pass method. I used a python script to execute the following command. But this will need to be edited to add your own sample paths.
bash star_alignment.sh <sample>
```

- Config file

You will need to modify the ```config/samples.csv``` file to add the sample names and paths to your input RNA BAM files.


## Outputs
Most files generated are temporary. Final results include filtered and annotated RNA VCFs, and features TSV matrix with variant labels predicted as one of: Germline, Somatic, or Artifact:
```
results/<sample>/VCFs/annotated/<sample>.annotated.vcf.gz
results/<sample>/Predictions/<sample>.annotated_predictions.csv
```


## How to run

Making the `qsub_logfiles` directory will allow snakemake to create a separate logfile for each rule. If this step is ommited it will write all logs to one file.
```
cd workflow
mkdir qsub_logfiles
qsub scripts/scheduler.sh
```




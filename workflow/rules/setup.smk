import sys, os
import pandas as pd
from snakemake.utils import min_version

min_version("6.0")

###### Config file and sample sheet #####
configfile: "../config/config.yaml"
rna_samples = pd.read_csv(config["samples"], comment="#").set_index("sample", drop=False)

##### Wildcard constraints #####
wildcard_constraints:
    sample="|".join(rna_samples.index),
    chrom = "|".join(["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"])

##### Helper functions #####
def get_rna_sample_bams(wildcards):
    """Get input RNA BAM file of a given sample."""
    return(rna_samples.loc[wildcards.sample, "file_path"])

def get_chroms():
    """Get list of chromosomes to paralellize variant calling."""
    return(["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"])


# def get_haplotype_caller_params(wildcards):
#     """Get all parameters for GATK HaplotypeCaller."""
#     return(
#         config["params"]["gatk"]["HaplotypeCaller"]["misc"] + 
#         " -L {chrom} ".format(chrom=wildcards.chrom) + 
#         config["params"]["gatk"]["HaplotypeCaller"]["annotations"]
#     )

def get_bcftools_query_values():
    """Get all values to query and format them for bcftools query."""
    return(config["querying"]["default"].replace(" ", "\\t"))

def get_all_inputs():
    """Get input files for rule all."""
    return(
        expand("../results/{sample}/Predictions/{sample}.annotated_predictions.csv", sample=rna_samples.index),
    )








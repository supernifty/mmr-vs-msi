
## Look for MMR mutations
This pipeline takes VCF files and finds candidates for MMR deficiency. 
This is achieved by:
* filtering input VCFs on MMR genes 
* annotating using VEP
* filtering on GnomAD, Impact, and SIFT

The result is a list of samples that are candidates for MMR deficient.

## Look for MSI
Also want to estimate the level of MSI in each sample.
Method:
* Generate a bed file from a list of MSI sites
* Are there any indels in these regions?

Answer the question:
* How likely is it that this genome is affected by MSI?

## Sources
* Refseq exons from UCSC

## Installation
```
conda create --name msi python=3.6
source activate msi
conda install --name msi snakemake
conda install --name msi cyvcf2
```

## Setup
In cfg:
* cluster.json: slurm configuration
* config.yaml: points to reference files, and all samples to be processed. 
* genes.txt: list of MMR genes to consider
* msi.bed: list of MSI regions to consider
* sample-metadata.csv: metadata used for patient ID and sample type
* samples: list of samples used by make_inputs

In reference:
* genome.fa
* refseq

## Running
First time:
```
src/make_inputs.py # to create symlinks in the in directory and makes config.yaml
```

Subsequently:
```
run.sh # kicks off the snakemake script
```

## Implementation
* Snakefile contains the main pipelinee

## TODO/Ideas
* gridss breakends are doubled if they are both within the gene (e.g. deletion)
* haplotype to find double deactivation
* calculate hypermutation -> total number of mutations per megabase
* cadd score
* filter on dbsnp, 1000GP, GnomAd

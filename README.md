
## Look for MMR mutations

## Sources
* Exon list if from Ensembl biomart

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
* config.yaml: points to reference files, and all samples to be processed
* genes.txt: list of MMR genes to consider
* sample-metadata.csv: metadata used for patient ID and sample type
* samples: list of samples used by make_inputs

In reference:
* genome.fa
* refseq

## Running
First time:
* run src/make_inputs.sh to create symlinks in the in directory

Subsequently:
* run.sh kicks off the snakemake script

## Implementation
* Snakefile contains the main pipelinee

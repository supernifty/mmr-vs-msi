
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
cd software/ && wget https://github.com/supernifty/mutational_signature/archive/0.1.tar.gz && tar xvfz 0.1.tar.gz && pip install -r mutational_signature-0.1/requirements.txt && cd -
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

### Generating reference/msi.candidates.bed
* Use https://github.com/supernifty/tandem_repeats_nim
```
./main --repeat:1 --min:6 < genome.fa >tb1.out
./main --repeat:2 --min:6 < genome.fa >tb2.out
./main --repeat:3 --min:6 < genome.fa >tb3.out
./main --repeat:4 --min:6 < genome.fa >tb4.out
./main --repeat:5 --min:6 < genome.fa >tb5.out
```

* Use src/annotate.py
```
python src/annotate.py --name exon tb1.out tb2.out tb3.out tb4.out tb5.out < hg19.exons > ./msi.exons
python src/annotate.py --name bethesda msi.exons < bethesda.bed | sort -k1,1 -k2,2n > reference/msi.candidates.bed
```

## Dependencies
* pindel

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

## Outputs
* msi.cluster.tumours: for each exonic MS loci, which tumour samples have indel mutations?
* msi.sample.cluster: for each exonic MS loci, which samples have indel mutations?
* msi.summary: MMR status, MSI context info, and MSI mutation summary.
* msi.cluster.?.tsv: loci mutations for different ranges.

* msi.samples.?.png: clustering based on the mutated MS loci in msi.sample.cluster
* msi.cluster.?.png: clustering bsaed on mutated loci but only considering indels of particular length
* msi.tumours.?.png: clustering based on msi.cluster.tumours

## TODO/Ideas
* gridss breakends are doubled if they are both within the gene (e.g. deletion)
* haplotype to find double deactivation
* calculate hypermutation -> total number of mutations per megabase
* cadd score
* filter on dbsnp, 1000GP, GnomAd
* mutational signatures

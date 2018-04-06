#!/bin/bash

# parameters:
# $1: sample name
# $2: reference
# $3: regions

SAMPLE="$1"
REFERENCE="$2"
REGIONS="$3"
INPUT_DIR="in"
OUTPUT_DIR="out"

set -e

echo "starting $SAMPLE at $(date)"

VEPPATH=/vlsci/UOM0040/shared/km/programs/ensembl-vep
CACHE=/vlsci/UOM0040/shared/km/programs/ensembl-vep/data/
THREADS=4

export PERL5LIB=$PERL5LIB:/vlsci/VR0002/kmahmood/Programs/vep/vep_87/ensembl-tools-release-87/scripts/variant_effect_predictor
  
module load bedtools-intel/2.27.1
module load perl/5.18.0
module load samtools-intel/1.5

#################
# mmr 
#################

# strategy
# - filter vcf on mmr exons
# - filter on quality (?)
# - run vep
# - filter on gnomad
# - filter on deleteriousness (?)
# - classify as MMRd

#################
# - filter vcf on mmr exons
bedtools intersect -a $INPUT_DIR/${SAMPLE}.varscan.vcf -b $REGIONS -header -wa -u > ${OUTPUT_DIR}/${SAMPLE}.mmr.vcf

#################
# - run vep

echo "starting vep for $SAMPLE at $(date)"

INPUT=${OUTPUT_DIR}/${SAMPLE}.mmr.vcf
OUTPUT_VEP=${OUTPUT_DIR}/${SAMPLE}.mmr.varscan.vep.vcf

$VEPPATH/vep --cache --refseq --offline --dir_cache  $CACHE --fasta $REFERENCE \
    -i $INPUT --sift b --polyphen b --symbol --numbers --biotype --total_length --hgvs \
    --format vcf -o $OUTPUT_VEP --force_overwrite --vcf --exclude_predicted \
    --af_gnomad \
    --fields Consequence,IMPACT,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,HGVSc,HGVSp,cDNA_position,CDS_position,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,PICK \
    --fork $THREADS --flag_pick

echo "finished $SAMPLE at $(date)"

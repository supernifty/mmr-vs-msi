#!/bin/bash

# parameters:
# $1: reference
# $2: regions
# $3: input
# $4: output

REFERENCE="$1"
REGIONS="$2"
INPUT="$3"
OUTPUT="$4"

set -e

module load bedtools-intel/2.27.1

# - filter vcf on msi regions
bedtools intersect -a "$INPUT" -b $REGIONS -header -wa -u > "$OUTPUT"

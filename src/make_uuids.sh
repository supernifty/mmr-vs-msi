#!/usr/bin/env bash

module load samtools-intel/1.5

for f in /scratch/VR0211/pan-prostate/out/*.validation; do
  # take the last line and the 9th column
  uuid=$(tail -1 $f | cut -f9)
  filename=$(basename $f)
  sample=${filename/.validation/}
  uuid_bam=$(samtools view -H /scratch/VR0211/pan-prostate/out/$sample.mapped.bam | grep "^@RG" | tail -1 | cut -f3 | sed 's/SM://')
  if [ "$uuid" != "$uuid_bam" ]; then
    (>&2 echo "WARN: different UUIDs for $sample: $uuid $uuid_bam")
  fi
  echo "$sample,$uuid_bam"
done

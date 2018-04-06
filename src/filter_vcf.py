#!/usr/bin/env python

import collections
import sys

import cyvcf2

AF_THRESHOLD=0.0001
SIFT_THRESHOLD=0.05

VEP="Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|PICK".split('|')

gnomad_index = VEP.index('gnomAD_AF')
impact_index = VEP.index('IMPACT')
sift_index = VEP.index('SIFT')

stats = collections.defaultdict(int)

vcf_in = cyvcf2.VCF('-')

sys.stdout.write(vcf_in.raw_header)

for variant in vcf_in:
  csq_list = variant.INFO.get('CSQ')
  include = False
  for csq in csq_list.split(','):
    csq_fields = csq.split('|')
    af = csq_fields[gnomad_index]
    impact = csq_fields[impact_index]
    sift = csq_fields[sift_index]

    # check impact
    if impact == 'LOW':
      stats['low_impact'] += 1
      continue

    if impact == 'MODIFIER':
      stats['modifier_impact'] += 1
      continue

    if sift != '' and sift.startswith('tolerated'):
      stats['tolerated_sift'] += 1
      continue

    # check allele frequency
    if af == '':
      include = True
      stats['empty_af'] += 1
      break

    if float(af) <= AF_THRESHOLD:
      include = True
      stats['low_af'] += 1
      break

    else:
      stats['high_af'] += 1

  if include:
      sys.stdout.write(str(variant))

sys.stderr.write(str(stats))

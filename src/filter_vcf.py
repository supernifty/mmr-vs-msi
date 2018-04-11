#!/usr/bin/env python

import collections
import sys

import cyvcf2

VERBOSE=True

AF_THRESHOLD=0.0001
SIFT_THRESHOLD=0.05

VEP="Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|PICK".split('|')

def log(msg):
    if VERBOSE:
        sys.stderr.write('{}\n'.format(msg))

def main(caller):
  gnomad_index = VEP.index('gnomAD_AF')
  impact_index = VEP.index('IMPACT')
  sift_index = VEP.index('SIFT')
  
  stats = collections.defaultdict(int)
  
  vcf_in = cyvcf2.VCF('-')
  
  sys.stdout.write(vcf_in.raw_header)
  
  for variant in vcf_in:
    # gridss is unfiltered
    if caller in ['gridss']:
      log('writing {} {}'.format(variant.CHROM, variant.POS))
      sys.stdout.write(str(variant))
      stats['written'] += 1
      continue # next variant

    csq_list = variant.INFO.get('CSQ')
    include = False
    log('assessing variant {} {}...'.format(variant.CHROM, variant.POS))
    for j, csq in enumerate(csq_list.split(',')): # each transcript
      csq_fields = csq.split('|')
      af = csq_fields[gnomad_index]
      impact = csq_fields[impact_index]
      sift = csq_fields[sift_index]
  
      # check impact
      if impact == 'LOW':
        stats['low_impact'] += 1
        log('skipped transcript {} due to LOW impact'.format(j))
        continue # next transcript
  
      if impact == 'MODIFIER':
        stats['modifier_impact'] += 1
        log('skipped transcript {} due to MODIFIER impact'.format(j))
        continue # next transcript
  
      if sift != '' and sift.startswith('tolerated'):
        stats['tolerated_sift'] += 1
        log('skipped transcript {} due to tolerated sift: {}'.format(j, sift))
        continue # next transcript
  
      # check allele frequency
      if af == '':
        stats['empty_af'] += 1
        include = True
        break # print
  
      if float(af) <= AF_THRESHOLD:
        stats['low_af'] += 1
        include = True
        break # print
  
      else:
        stats['high_af'] += 1
        log('skipped transcript {} due to high af: {}'.format(j, af))
        continue # next transcript
  
    if include:
        log('writing {} {}'.format(variant.CHROM, variant.POS))
        sys.stdout.write(str(variant))
        stats['written'] += 1
  
  sys.stderr.write(str(stats))

if __name__ == '__main__':
  main(sys.argv[1])

#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import os.path
import sys

import cyvcf2

def main(vcfs):
  logging.info('starting...')

  line = ['sample', 'chrom', 'pos', 'ref', 'alt', 'qual', 'af', 'dp', 'gene', 'oncogene', 'exon', 'annotation']
  sys.stdout.write('{}\n'.format('\t'.join(line)))
  for vcf in vcfs:
    logging.info('processing %s...', vcf)
    sample = os.path.basename(vcf).split('.')[0]

    vcf_in = cyvcf2.VCF(vcf)
    for variant in vcf_in:
      try:
        annotation = variant.INFO['VD'] # vagrent default annotation e.g. CCNL2|CCDS30558.1|r.984delA|c.?|p.?|protein_coding:exon:three_prime_UTR:deletion:3_prime_UTR_variant:transcript_variant|SO:0000010:SO:0000147:SO:0000205:SO:0000159:SO:0001624:SO:0001576
        detail = annotation.split('|')[5]
      except KeyError:
        detail = ''
      try:
        oncogene = variant.INFO['msi_all_oncogene']
      except KeyError:
        oncogene = '' 
      gene = variant.INFO['msi_exon'].split('_')[0]
      #line = [sample, variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), '{:.2f}'.format(variant.INFO['AF']), str(variant.INFO['DP']), gene, oncogene, variant.INFO['msi_exon'], detail]
      #TODO af and dp
      line = [sample, variant.CHROM, str(variant.POS), variant.REF, ','.join(variant.ALT), str(variant.QUAL), '0', '0', gene, oncogene, variant.INFO['msi_exon'], detail]
      sys.stdout.write('{}\n'.format('\t'.join(line)))
    logging.info('processing %s: done', vcf)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcfs', required=True, nargs='+', help='input vcfs')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.vcfs)

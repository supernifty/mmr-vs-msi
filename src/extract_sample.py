#!/usr/bin/env python
'''
  extracts just the calls for the provided sample name from matched calls
  note that this does NOT take into account the genotype
'''

import argparse
import logging
import sys

import cyvcf2

UUIDS="cfg/samples.uuids"

def main(sample, qual, pindel, filter_homref, strelka):
  logging.debug('reading vcf from stdin. qual filter %i, pindel filter %.2f', qual, pindel)
  vcf_in = cyvcf2.VCF('-')
  if sample not in vcf_in.samples:
    logging.info('{} not in samples: {}'.format(sample, vcf_in.samples))
    found = False
    sample_name = sample.split('.')[0]
    for uuid in open(UUIDS, 'r'):
      fields = uuid.strip('\n').split(',')
      if fields[0] == sample_name and fields[1] in vcf_in.samples:
        sample = fields[1]
        logging.info('using {} instead'.format(sample))
        found = True
        break
    if not found:
      logging.error('unable to determine correct sample name')
      return
  vcf_in.set_samples(sample)
  sys.stdout.write(vcf_in.raw_header)

  allowed = 0
  denied = 0
  for variant in vcf_in:
    ok = True
    # gridss
    if 'QUAL' in variant.FORMAT:
      variant_qual = variant.format('QUAL')
      if variant_qual < qual:
        ok = False

    # pindel
    if all([x in variant.FORMAT for x in ('PP', 'NP', 'PR', 'NR')]):
      if variant.format('PP') / variant.format('PR') < pindel and variant.format('NP') / variant.format('NR') < pindel:
        ok = False

    # strelka_indel
    if all([x in variant.FORMAT for x in ('TAR', 'TIR')]):
      tir = sum(variant.format('TIR')[0])
      tar = sum(variant.format('TAR')[0])
      if tar == 0 or tir / (tir + tar) < strelka:
        ok = False

    # strelka
    if variant.INFO.get('AF') is not None:
      if variant.INFO.get('AF') < strelka:
        ok = False

    # platypus
    # check gt 0,1,2,3==HOM_REF (0/0), HET (0/1), UNKNOWN (./.), HOM_ALT (1/1)
    if filter_homref:
      gt = variant.gt_types[0]
      if gt == 0:
        ok = False

    if ok:
      sys.stdout.write(str(variant))
      allowed += 1
    else:
      denied += 1

    if (allowed + denied) % 100000 == 0:
      logging.debug('%i processed. %i allowed.', allowed + denied, allowed)
    
  logging.info('done. wrote {}. skipped {}. total {}'.format(allowed, denied, allowed + denied))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract samples from VCF')
  parser.add_argument('--sample', required=True, help='sample to extract')
  parser.add_argument('--minqual', required=False, type=int, default=-1, help='minimum sample quality')
  parser.add_argument('--minpindel', required=False, type=float, default=0.01, help='minimum proportion pindel call to mapped reads')
  parser.add_argument('--minstrelka', required=False, type=float, default=0.01, help='minimum proportion strelka call to mapped reads')
  parser.add_argument('--filter_homref', action='store_true', help='filter homref')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.sample, args.minqual, args.minpindel, args.filter_homref, args.minstrelka)

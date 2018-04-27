#!/usr/bin/env python
'''
  extracts just the calls for the provided sample name from matched calls
'''

import argparse
import logging
import sys

import cyvcf2

QUAL=500
UUIDS="cfg/samples.uuids"

def main(sample):
  logging.debug('reading vcf from stdin')
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
    qual = variant.format('QUAL')
    if qual >= QUAL:
      sys.stdout.write(str(variant))
      allowed += 1
    else:
      denied += 1
    
  logging.info('done. wrote {}. skipped {}.'.format(allowed, denied))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract samples from VCF')
  parser.add_argument('--sample', required=True, help='sample to extract')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.sample)

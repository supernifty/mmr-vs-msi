#!/usr/bin/env python
'''
  subtract calls from a vcf
'''

import argparse
import logging
import sys

import cyvcf2

def main(germline):
  logging.info('reading vcf from %s...', germline)

  baseline = set()
  for variant in cyvcf2.VCF(germline):
    key = '{chrom}\t{pos}\t{ref}\t{alt}'.format(chrom=variant.CHROM, pos=variant.POS, ref=variant.REF, alt=','.join(variant.ALT))
    baseline.add(key)
  logging.info('reading vcf from %s: %i variants added', germline, len(baseline))

  logging.info('reading source vcf from stdin...')
  filtered = 0
  total = 0
  source = cyvcf2.VCF('-')
  sys.stdout.write(source.raw_header)
  for total, variant in enumerate(source):
    key = '{chrom}\t{pos}\t{ref}\t{alt}'.format(chrom=variant.CHROM, pos=variant.POS, ref=variant.REF, alt=','.join(variant.ALT))
    if key in baseline:
      filtered += 1
    else:
      sys.stdout.write(str(variant))

  logging.info('done. filtered out %i variants from a total of %i', filtered, total + 1)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Subtract calls from a source')
  parser.add_argument('--subtract', required='true', help='vcf to subtract')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.subtract)

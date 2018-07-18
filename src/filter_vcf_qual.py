#!/usr/bin/env python
'''
  extracts just the calls for the provided sample name from matched calls
'''

import argparse
import logging
import sys

import cyvcf2

UUIDS="cfg/samples.uuids"

def main(qual):
  logging.info('reading vcf from stdin. qual filter %i', qual)
  vcf_in = cyvcf2.VCF('-')
  sys.stdout.write(vcf_in.raw_header)

  allowed = 0
  denied = 0
  for variant in vcf_in:
    ok = variant.QUAL is None or variant.QUAL >= qual
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
  parser.add_argument('--qual', required=True, type=int, default=-1, help='minimum sample quality')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.qual)

#!/usr/bin/env python
'''
  tool to detect artifacts by counting common mutations
'''

import argparse
import collections
import logging
import os
import os.path
import sys

import cyvcf2
import intervaltree

def main(threshold, common_in, out, position_only):
  # TODO position_only = true only supported
  common = set()
  first = True
  logging.info('reading common variants from %s', common_in)
  for line in open(common_in, 'r'):
    if first:
      first = False
      continue
    fields = line.strip('\n').split('\t')
    prop = float(fields[3])
    if prop >= threshold:
      common.add('{}\t{}'.format(fields[0], fields[1]))

  logging.info('reading vcf from stdin')
  vcf = cyvcf2.VCF('-')
  filtered = total = 0
  out.write(vcf.raw_header)
  for total, variant in enumerate(vcf):
    if '{}\t{}'.format(variant.CHROM, variant.POS) in common:
      filtered += 1
    else:
      out.write(str(variant))
    
  logging.info('filtered %i of %i', filtered, total + 1)

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='Filter common variations')
  parser.add_argument('--threshold', type=float, default=0.2, help='min proportion to include')
  parser.add_argument('--common', required=True, help='filename of common variations')
  parser.add_argument('--position_only', action='store_true', help='list of vcfs')
  args = parser.parse_args()
  main(args.threshold, args.common, sys.stdout, args.position_only)

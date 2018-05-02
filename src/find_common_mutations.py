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

def main(vcfs, out, threshold, position_only):
  logging.info('reading %i vcf files. threshold is %.2f (seen at least %.1f times)...', len(vcfs), threshold, int(threshold * len(vcfs)))
  counts = collections.defaultdict(int)
  for idx, vcf in enumerate(vcfs):
    logging.info('processing %s: %i of %i...', vcf, idx + 1, len(vcfs))
    for variant in cyvcf2.VCF(vcf):
      if position_only:
        identity = '{chrom}\t{pos}'.format(chrom=variant.CHROM, pos=variant.POS)
      else:
        identity = '{chrom}\t{pos}\t{ref}\t{alt}'.format(chrom=variant.CHROM, pos=variant.POS, ref=variant.REF, alt=variant.ALT[0])
      counts[identity] += 1
    logging.info('processing %s: %i of %i. total %i distinct variants', vcf, idx + 1, len(vcfs), len(counts))

  logging.info('generating result...')
  if position_only:
    out.write('{}\t{}\t{}\t{}\n'.format('Chrom', 'Pos', 'Count', 'Proportion'))
  else:
    out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('Chrom', 'Pos', 'Ref', 'Alt', 'Count', 'Proportion'))
  included = excluded = 0
  for key in sorted(counts.keys()):
    fields = key.split('\t')
    if counts[key] >= threshold: 
      proportion = counts[key] / len(vcfs)
      if position_only:
        out.write('{}\t{}\t{}\t{:.2f}\n'.format(fields[0], fields[1], counts[key], proportion))
      else:
        out.write('{}\t{}\t{}\t{}\t{}\t{:.2f}\n'.format(fields[0], fields[1], fields[2], fields[3], counts[key], proportion))
      included += 1
    else:
      excluded += 1

  if included + excluded == 0:
    logging.info('no variants found')
  else:
    logging.info('wrote %i common variants and excluded %i rare variants out of %i total. That\'s %.1f%%', included, excluded, included + excluded, 100 * included / (included + excluded))

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='Count variations')
  parser.add_argument('--vcfs', required=True, nargs='+', help='list of vcfs')
  parser.add_argument('--position_only', action='store_true', help='list of vcfs')
  parser.add_argument('--threshold', type=int, default=2, help='min count to include')
  args = parser.parse_args()
  main(args.vcfs, sys.stdout, args.threshold, args.position_only)

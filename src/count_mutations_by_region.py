#!/usr/bin/env python

import argparse
import collections
import logging
import os
import os.path
import sys

import cyvcf2
import intervaltree

def main(bed, vcfs, out):
  logging.info('reading bed file...')
  intervals = {}
  for idx, line in enumerate(open(bed, 'r')):
    chr, start, finish, annotation = line.strip('\n').split('\t', 3)
    if chr not in intervals:
      intervals[chr] = intervaltree.IntervalTree()
    intervals[chr][int(start):int(finish)] = (annotation, set())
    if (idx + 1) % 100000 == 0:
      logging.debug('%i lines...', idx + 1)

  logging.info('reading %i vcf files...', len(vcfs))
  samples = set()
  totals = collections.defaultdict(int)
  for idx, vcf in enumerate(vcfs):
    logging.info('processing %s: %i of %i...', vcf, idx + 1, len(vcfs))
    sample = os.path.basename(vcf).split('.')[0]
    samples.add(sample)
    for variant in cyvcf2.VCF(vcf):
      if variant.CHROM in intervals:
        for interval in intervals[variant.CHROM][variant.POS]:
          interval.data[1].add(sample)
          totals[sample] += 1
    
  logging.info('generating result...')
  out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Chrom', 'Begin', 'End', 'Annotation', 'Len', 'Samples', '\t'.join(sorted(samples))))
  for chrom in sorted(intervals):
    for interval in sorted(intervals[chrom]):
      annotation, matches = interval.data
      out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, interval.begin, interval.end, annotation, interval.end - interval.begin, len(matches), '\t'.join(['1' if sample in matches else '0' for sample in sorted(samples)])))

  out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('TOTAL', '0', '0', 'ALL', '0', len(samples), '\t'.join([str(totals[sample]) for sample in sorted(samples)])))

  logging.info('done')

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  parser = argparse.ArgumentParser(description='Given a bed file and list of vcfs, count variations by region')
  parser.add_argument('--bed', required=True, help='list of regions')
  parser.add_argument('--vcfs', nargs='+', help='list of vcfs')
  args = parser.parse_args()
  main(args.bed, args.vcfs, sys.stdout)

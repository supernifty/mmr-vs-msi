#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
  do this by enumerating all indels in germline and tumour samples in MS regions.
  if there is a difference, mark it as evidence for MSI 
'''

import argparse
import collections
import logging
import os.path
import sys

import cyvcf2
import intervaltree


def count_as_str(count):
  return ' '.join(['{}:{}'.format(x, count[x]) for x in sorted(count.keys())])

def net(count):
  return sum([x * count[x] for x in count])

def main(regions, tumour_vcf, normal_vcf, screen, summary):
  logging.info('starting...')

  logging.debug('reading %s...', regions)
  count = 0
  intervals = {}
  for count, line in enumerate(open(regions, 'r')):
    chr, start, finish, annotation = line.strip('\n').split('\t', 3)
    if chr not in intervals:
      intervals[chr] = intervaltree.IntervalTree()
    intervals[chr][int(start):int(finish)] = {'annot': annotation, 'normal': collections.defaultdict(int), 'tumour': collections.defaultdict(int)}
    if (count + 1) % 100000 == 0:
      logging.debug('%i lines...', count + 1)
  logging.debug('reading %s: %i lines processed', regions, count + 1)

  logging.debug('reading %s...', normal_vcf)
  count = 0
  added = 0
  normal = {}
  for count, variant in enumerate(cyvcf2.VCF(normal_vcf)):
    indel_len = len(variant.ALT[0]) - len(variant.REF)
    if indel_len != 0:
      if variant.CHROM in intervals:
        for interval in intervals[variant.CHROM][variant.POS]:
          added += 1
          interval.data['normal'][indel_len] += 1
  logging.debug('%s: done reading %i variants, %i added', normal_vcf, count + 1, added)

  logging.debug('reading %s...', tumour_vcf)
  count = 0
  added = 0
  for count, variant in enumerate(cyvcf2.VCF(tumour_vcf)):
    indel_len = len(variant.ALT[0]) - len(variant.REF)
    if indel_len != 0:
      if variant.CHROM in intervals:
        for interval in intervals[variant.CHROM][variant.POS]:
          interval.data['tumour'][indel_len] += 1
          added += 1
  logging.debug('%s: done reading %i variants, %i added', tumour_vcf, count + 1, added)
  
  
  logging.debug('finding differences...')
  same = 0
  different = 0
  for chrom in sorted(intervals.keys()):
    for interval in sorted(intervals[chrom]):
      if count_as_str(interval.data['normal']) == count_as_str(interval.data['tumour']):
        same += 1
      else:
        different += 1
        net_diff = net(interval.data['tumour']) - net(interval.data['normal'])
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, interval.begin, interval.end, interval.data['annot'], count_as_str(interval.data['normal']), count_as_str(interval.data['tumour']), net_diff))

  if summary is not None:
    with open(summary, 'w') as fh:
      sample = os.path.basename(tumour_vcf).split('.')[0]
      ratio = different / (same + different)
      fh.write('{}\t{}\t{}\t{:.3f}\n'.format(sample, same, different, ratio))
    
  logging.info('done. %i same. %i different', same, different)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--tumour', required=True, help='tumour vcf')
  parser.add_argument('--normal', required=True, help='normal vcf')
  parser.add_argument('--screen', required=False, help='screen to filter common mutations')
  parser.add_argument('--regions', required=True, help='msi regions')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--summary', required=False, help='write summary to this file')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.regions, args.tumour, args.normal, args.screen, args.summary)

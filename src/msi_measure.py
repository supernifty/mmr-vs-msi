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

def find_normal(metadata, sample_prefix):
  prefix, sample = sample_prefix.rsplit('/', 1)
  
  patient = None
  normals = {}
  # CMHS1,CMHP1,299bT,BT,N,,,
  for line in metadata:
    fields = line.strip('\n').split(',')
    if fields[0] == sample:
      patient = fields[1]
    if fields[4] == 'Y':
      normals[fields[1]] = fields[0]

  if patient in normals:
    return '/'.join([prefix, normals[patient]])
  else:
    logging.warn('Sample {} with patient {} not found in metadata'.format(sample, patient))
    return None

def main(regions, tumour_vcfs, screen, summary, metadata):
  logging.info('starting...')

  logging.debug('reading regions from %s...', regions)
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

  # find normal from metadata
  tumour_sample, rest = tumour_vcfs[0].split('.', 1)
  normal_sample = find_normal(open(metadata, 'r'), tumour_sample)

  # make a list of indels found in each msi region, for the normal vcf
  for tumour_vcf in tumour_vcfs:
    tumour_sample, rest = tumour_vcf.split('.', 1)
    normal_vcf = '.'.join([normal_sample, rest])

    logging.debug('reading %s...', normal_vcf)
    
    count = 0
    added = 0
    for count, variant in enumerate(cyvcf2.VCF(normal_vcf)):
      indel_len = len(variant.ALT[0]) - len(variant.REF)
      if indel_len != 0:
        if variant.CHROM in intervals:
          for interval in intervals[variant.CHROM][variant.POS]:
            added += 1
            interval.data['normal'][indel_len] += 1
  
    logging.debug('%s: done reading %i variants, %i added', normal_vcf, count + 1, added)

  # make a list of indels found in each msi region, for the tumour vcf
  for tumour_vcf in tumour_vcfs:
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
  
  # find msi regions with different indel counts in the tumour relative to the normal
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
  parser.add_argument('--metadata', required=True, help='metadata file')
  parser.add_argument('--tumours', required=True, nargs='+', help='tumour vcf')
  parser.add_argument('--screen', required=False, help='screen to filter common mutations')
  parser.add_argument('--regions', required=True, help='msi regions')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--summary', required=False, help='write summary to this file')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.regions, args.tumours, args.screen, args.summary, args.metadata)

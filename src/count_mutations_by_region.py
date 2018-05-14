#!/usr/bin/env python

import argparse
import collections
import csv
import logging
import os
import os.path
import sys

import cyvcf2
import intervaltree

def main(bed, vcfs, out, use_lengths, minlength, maxlength, minsamples, maxsamplespct, categories, category):
  logging.info('reading bed file...')
  intervals = {}
  for idx, line in enumerate(open(bed, 'r')):
    chr, start, finish, annotation = line.strip('\n').split('\t', 3)
    if chr not in intervals:
      intervals[chr] = intervaltree.IntervalTree()
    intervals[chr][int(start):int(finish)] = (annotation, dict())
    if (idx + 1) % 100000 == 0:
      logging.debug('%i lines...', idx + 1)

  allowed = set()
  if categories is not None:
    with open(categories, 'rt') as csvfile:
      for row in csv.reader(csvfile):
        if row[3] == category:
          allowed.add(row[0])

  logging.info('reading %i vcf files...', len(vcfs))
  samples = set()
  totals = collections.defaultdict(int)
  length_ok = length_notok = 0
  for idx, vcf in enumerate(vcfs):
    logging.info('processing %s: %i of %i...', vcf, idx + 1, len(vcfs))
    sample = os.path.basename(vcf).split('.')[0]
    if len(allowed) > 0 and sample not in allowed:
      continue
    samples.add(sample)
    for variant in cyvcf2.VCF(vcf):
      variant_length = len(variant.ALT[0]) - len(variant.REF)
      if minlength <= variant_length <= maxlength:
        length_ok += 1
        if variant.CHROM in intervals:
          for interval in intervals[variant.CHROM][variant.POS]: # all overlapping intervals
            if args.use_lengths:
              interval.data[1][sample] = variant_length # assuming only one variant per interval
            else:
              interval.data[1][sample] = 1
            totals[sample] += 1 # variants per sample
      else:
        length_notok += 1

  logging.info('generating result from %i passed variants (%i outside length restriction)...', length_ok, length_notok)
  out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Chrom', 'Begin', 'End', 'Annotation', 'Len', 'Samples', '\t'.join(sorted(samples))))
  regions_ok = regions_skipped = 0
  for chrom in sorted(intervals):
    for interval in sorted(intervals[chrom]):
      annotation, matches = interval.data
      if len(matches) >= minsamples and len(matches) <= len(vcfs) * maxsamplespct / 100:
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, interval.begin, interval.end, annotation, interval.end - interval.begin, len(matches), '\t'.join([str(matches[sample]) if sample in matches else '0' for sample in sorted(samples)])))
        regions_ok += 1
      else:
        regions_skipped += 1

  logging.info('wrote %i and skipped %i regions with samples < %i or > %i%%', regions_ok, regions_skipped, minsamples, maxsamplespct)

  out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('TOTAL', '0', '0', 'ALL', '0', len(samples), '\t'.join([str(totals[sample]) for sample in sorted(samples)])))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Given a bed file and list of vcfs, count variations by region')
  parser.add_argument('--bed', required=True, help='list of regions')
  parser.add_argument('--use_lengths', action='store_true',  help='write lengths of indels rather than 1')
  parser.add_argument('--minsamples', required=False, default=0, type=int, help='min number of samples required to print region')
  parser.add_argument('--maxsamplespct', required=False, default=100, type=int, help='max % of samples allowed to print region')
  parser.add_argument('--minlength', required=False, default=-100000, type=int, help='min length of an indel')
  parser.add_argument('--maxlength', required=False, default=100000, type=int, help='max length of an indel')
  parser.add_argument('--vcfs', nargs='+', help='list of vcfs')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--categories', required=False, help='CSV of additional info about the sample')
  parser.add_argument('--category', required=False, help='filter on category')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.bed, args.vcfs, sys.stdout, args.use_lengths, args.minlength, args.maxlength, args.minsamples, args.maxsamplespct, args.categories, args.category)

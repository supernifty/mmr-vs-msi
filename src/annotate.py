#!/usr/bin/env python

# Usage
# python annotate.py --name exon tb1.chr22.bed < hg19.chr22.exons > tb1.chr22.annotated.bed

# chr1    11868   12227   LOC102725121_exon_0_0_chr1_11869_f      0       +

import argparse
import collections
import logging
import sys

import intervaltree


def main(name, annotations, beds, target, annotate_none):
  logging.info('reading annotations...')
  trees = {}
  for idx, line in enumerate(annotations):
    fields = line.strip('\n').split('\t')
    if len(fields) < 4:
      logging.warn('Too few fields on line %i: %s', idx, line.strip('\n'))
      continue
    chrom = fields[0]
    start = int(fields[1])
    finish = int(fields[2])
    data = fields[3]
    if chrom not in trees:
      trees[chrom] = intervaltree.IntervalTree()
      logging.info('added chromosome %s', chrom)
    trees[chrom][start:finish] = data
  logging.info('reading annotations: done')

  # now find overlaps
  found = 0
  skipped = 0
  missing = set()

  for bed in beds:
    logging.info('reading %s', bed)
    with open(bed, 'r') as src:
      for idx, line in enumerate(src):
        if idx % 1000000 == 0:
          logging.debug('read %i lines from %s', idx, bed)
        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        start = int(fields[1])
        finish = int(fields[2])
        data = fields[3]
        if chrom not in trees:
          overlaps = []
          if chrom not in missing:
            logging.warn('%s not found in annotation', chrom)
            missing.add(chrom)
        else:
          overlaps = trees[chrom][start:finish]
  
        if len(overlaps) > 0:
          found += 1
          # put all overlaps in the exon field
          combined_overlaps = ','.join(['{}[{}:{}]'.format(overlap.data, overlap.begin, overlap.end) for overlap in overlaps])
          sys.stdout.write('{}\t{}\t{}\t{};{}={}\n'.format(chrom, start, finish, data, name, combined_overlaps))
        else:
          if annotate_none:
            sys.stdout.write('{}\t{}\t{}\t{};{}=none\n'.format(chrom, start, finish, data, name))
          else:
            sys.stdout.write('{}\t{}\t{}\t{}\n'.format(chrom, start, finish, data))
          skipped += 1

  logging.info("done: found %i overlaps skipped %i", found, skipped)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='annotate bed file')
  parser.add_argument('--name', required=True, help='name of annotation field')
  parser.add_argument('--annotate_none', action='store_true', help='add none annotation if not found')
  parser.add_argument('beds', nargs='+', help='input bed files to annotate')
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  
  args = parser.parse_args()
  main(args.name, sys.stdin, args.beds, sys.stdout, args.annotate_none)


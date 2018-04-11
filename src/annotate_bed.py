#!/usr/bin/env python

'''
  adds annotation to first bed from second bed

  Usage:
  $0 name ann.bed < main.bed > main_plus_ann.bed
'''


import logging
import sys

import intervaltree

def main(name, bed, src, target):
  logging.info('parsing {}...'.format(bed))
  tree = {}
  for idx, line in enumerate(open(bed, 'r')):
    chrom, start, finish, annotation = line.strip('\n').split('\t')[:4]
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
    tree[chrom][int(start):int(finish)] = annotation
  logging.info('parsing {}: done.'.format(bed))

  logging.info('reading from stdin...')
  yes = no = 0
  for idx, line in enumerate(src):
    chrom, start, finish, annotation = line.strip('\n').split('\t')[:4]
    if chrom not in tree:
      sys.stdout.write(line)
      no += 1
      continue
    overlaps = tree[chrom][int(start):int(finish)]
    if len(overlaps) == 0:
      sys.stdout.write(line)
      no += 1
    else:
      values = set([overlap.data for overlap in overlaps])
      value = ','.join(values)
      if annotation.endswith(';'):
        sys.stdout.write('{}\t{}\t{}\t{}{}={}\n'.format(chrom, start, finish, annotation, name, value))
      else:
        sys.stdout.write('{}\t{}\t{}\t{};{}={}\n'.format(chrom, start, finish, annotation, name, value))
      yes += 1

  logging.info('%i overlaps. %i with no overlap.', yes, no)

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1], sys.argv[2], sys.stdin, sys.stdout)

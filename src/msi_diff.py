#!/usr/bin/env python
'''
  find mutations that segregate provided clusters
'''

import argparse
import logging
import sys

def main(fh, sample1, sample2, out):
  logging.info('reading cluster from stdin...')
  header = fh.readline().strip('\n').split('\t')
  empty = same = only1 = only2 = both = 0
  s1idx = []
  for sample in sample1:
    s1idx.append(header.index(sample))
  s2idx = []
  for sample in sample2:
    if sample in sample1:
      logging.debug('skipping %s in sample2', sample)
      continue
    s2idx.append(header.index(sample))

  sys.stdout.write('Cat\t{}\t{}\t{}\n'.format(sample1[0], sample2[0], '\t'.join(header[0:6])))
  for line in fh:
    fields = line.strip('\n').split('\t')
    if fields[0] == 'TOTAL':
      continue
    if all([fields[x] == '0' for x in s1idx]) and all([fields[x] == '0' for x in s2idx]):
      empty += 1
      cat = None
    elif fields[s1idx[0]] == fields[s2idx[0]] and all([fields[s1idx[0]] == fields[x] for x in s1idx]) and all(fields[s2idx[0]] == fields[x] for x in s2idx):
      same += 1
      cat = 'same'
    elif any([fields[x] != '0' for x in s1idx]) and all([fields[x] == '0' for x in s2idx]):
      only1 += 1
      cat = sample1[0]
    elif all([fields[x] == '0' for x in s1idx]) and any([fields[x] != '0' for x in s2idx]):
      only2 += 1
      cat = sample2[0]
    else:
      both += 1
      cat = 'both'

    if cat is not None:
      sys.stdout.write('{}\t{}\t{}\t{}\n'.format(cat, ','.join([fields[x] for x in s1idx]), ','.join([fields[x] for x in s2idx]), '\t'.join(fields[0:6])))

  sys.stdout.write('#summary\n#empty\t{}\n#same\t{}\n#only {}\t{}\n#only {}\t{}\n#both\t{}\n'.format(empty, same, sample1[0], only1, sample2[0], only2, both))

  logging.info('Summary: Empty: {}. Same: {}. Only {}: {}. Only {}: {}. Both: {}\n'.format(empty, same, sample1[0], only1, sample2[0], only2, both))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--sample1', required='true', nargs='+', help='samples in cluster')
  parser.add_argument('--sample2', required='true', nargs='+', help='samples in cluster')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.sample1, args.sample2, sys.stdout)


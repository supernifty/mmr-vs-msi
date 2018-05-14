#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys

def main(fh, match, out):
  logging.info('starting...')
  header = fh.readline().strip('\n').split('\t')
  sys.stdout.write('{}\tSamples\n'.format('\t'.join(header[0:6])))
  for line in fh:
    if match is None or match in line:
      fields = line.strip('\n').split('\t')
      samples = [header[i + 6] for i, x in enumerate(fields[6:]) if x != '0']
      sys.stdout.write('{}\t{}\n'.format('\t'.join(fields[0:6]), ','.join(samples)))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--match', required=False, help='filter on match')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.match, sys.stdout)

#!/usr/bin/env python
'''
  limit rows
'''

import argparse
import logging
import sys

def main(include):
  logging.info('starting...')
  line = sys.stdin.readline()
  sys.stdout.write(line)

  for line in sys.stdin:
    sample = line.strip('\n').split('\t')[0]
    if sample in include:
      sys.stdout.write(line)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--include', nargs='+', help='samples to include')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.include)

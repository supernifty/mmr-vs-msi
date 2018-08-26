#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys

def main(genes_fh):
  logging.info('starting...')
  genes = set()
  for line in open(genes_fh, 'r'):
    genes.add(line.strip('\n'))

  logging.debug('%i genes to include', len(genes))

  for line in sys.stdin:
    if len(line.strip('\n')) == 0:
      break

  logging.debug('reading candidates...')
  candidates = sys.stdin.readline() # candidates...

  header = sys.stdin.readline().strip('\n').split('\t')

  # we include 0..4
  genes_found = set()
  genes_available = set()
  columns = []
  for i in range(5, len(header)):
    gene = header[i].split('|')[0]
    genes_available.add(gene)
    if gene in genes:
      columns.append(i)
      genes_found.add(gene)

  logging.debug('found %i genes of %i requested from %i available', len(genes_found), len(genes), len(genes_available))

  row = header[0:5] + [header[x] for x in columns]
  sys.stdout.write('{}\n'.format('\t'.join(row)))

  logging.debug('reading data...')
  lines = 0
  for lines, line in enumerate(sys.stdin):
    fields = line.strip('\n').split('\t')
    row = fields[0:5] + [fields[x] for x in columns]
    sys.stdout.write('{}\n'.format('\t'.join(row)))

  logging.info('done reading %i rows', lines + 1)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--genes', required=True, help='list of genes to keep')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genes)

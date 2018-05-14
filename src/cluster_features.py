#!/usr/bin/env python
'''
  find mutations that segregate provided clusters
'''

import argparse
import logging
import sys

import pandas as pd

def main(fh, cluster1, cluster2=None):
  logging.info('reading cluster from stdin...')
  data = pd.read_csv(fh, sep='\t', dtype={'Chrom': object, 'Begin': int, 'End': int, 'Annotation': object, 'Len': int, 'Samples': int})

  # for each loci, want to count how many times it's affected for samples in cluster, and for those not in cluster
  Xt = data[:-1].T
  X = Xt[6:] # skip the last line; each sample is a row
  y = data.columns.values[6:] # skip Chrom..Samples

  # columns to include in each cluster
  for index, row in Xt.iterrows(): # each loci
    logging.debug(index)
    #if row[0] != 0:# and sum(v == 0 for v in row[1:]):
    logging.info('%i %i %s', row[0], sum(row[1:]), Xt.loc[index][:6])

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--cluster1', required='true', nargs='+', help='samples in cluster')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.cluster1)


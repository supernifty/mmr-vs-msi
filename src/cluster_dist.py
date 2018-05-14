#!/usr/bin/env python
'''
  cluster on ms location
'''

import argparse
import logging
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

import scipy.cluster.hierarchy

def main(fh, target):
  logging.info('starting...')

  #Chrom   Begin   End     Annotation      Len     Samples CMHS1   CMHS101 CMHS104 CMHS105 CMHS106 CMHS107 CMHS109 CMHS110...
  data = pd.read_csv(fh, sep='\t', dtype={'Chrom': object, 'Begin': int, 'End': int, 'Annotation': object, 'Len': int, 'Samples': int})
  logging.debug('input data has %i rows and %i columns', data.shape[0], data.shape[1])

  X = data[:-1].T[6:] # skip the last line; each sample is a row
  y = data.columns.values[6:] # skip Chrom..Samples

  logging.debug('X has %i rows and %i columns: %s', X.shape[0], X.shape[1], X.head())
  logging.debug('y is %s', y)

  logging.info('distance to sample {}...'.format(y[0]))
  
  sys.stdout.write('{}\t{}\t{}\t{}\n'.format('Sample', 'Distance', 'Overlap', 'Total'))
  for idx, row in X.iterrows(): # each sample
    distance = 0
    overlap = 0
    count = 0
    for pos in range(len(row)): # each loci
      candidate = row[pos]
      source = X.iloc[0][pos]
      if candidate == 0 and source != 0 or source == 0 and candidate != 0:
        distance += 1
      if candidate != 0:
        count += 1
      if candidate != 0 and source != 0:
        overlap += 1
    sys.stdout.write('{}\t{}\t{}\t{}\n'.format(idx, distance, overlap, count))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Cluster MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--image', default='dend.png', help='dend file')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.image)

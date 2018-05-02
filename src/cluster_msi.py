#!/usr/bin/env python
'''
  cluster on ms location
'''

import argparse
import collections
import csv
import logging
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

#import sklearn
#import sklearn.cluster
#import sklearn.preprocessing

import scipy.cluster.hierarchy

def main(fh, target, categories):
  logging.info('reading cluster results from stdin...')

  vcf_categories = collections.defaultdict(set)
  sample_categories = {}
  # "CMHS376","CMHP153","PP102_TURPBX","Primary","Library","ContEst","","Inferred barcode"
  if categories is not None:
    with open(categories, 'rt') as csvfile:
      for row in csv.reader(csvfile):
        vcf_categories[row[3]].add(row[0])
        sample_categories[row[0]] = row[3]

  #Chrom   Begin   End     Annotation      Len     Samples CMHS1   CMHS101 CMHS104 CMHS105 CMHS106 CMHS107 CMHS109 CMHS110...
  data = pd.read_csv(fh, sep='\t', dtype={'Chrom': object, 'Begin': int, 'End': int, 'Annotation': object, 'Len': int, 'Samples': int})
  logging.debug('input data has %i rows and %i columns', data.shape[0], data.shape[1])

  X = data[:-1].T[6:] # skip the last line; each sample is a row
  y = data.columns.values[6:] # skip Chrom..Samples

  if len(sample_categories) > 0:
    y_all = ['{} {}'.format(s, sample_categories[s]) for s in y]
  else:
    y_all = y

  logging.debug('X has %i rows and %i columns: %s', X.shape[0], X.shape[1], X.head())
  logging.debug('y is %s', y)

  logging.info('clustering...')
  #km = sklearn.cluster.KMeans(n_clusters=4, init='k-means++', max_iter=100, n_init=1, verbose=True)
  #km.fit(X)

  Z = scipy.cluster.hierarchy.linkage(X, method='ward', metric='euclidean')

  logging.info('plotting...')
  plt.figure(figsize=(30, 12))
  plt.title('Hierarchical Clustering Dendrogram - All samples')
  plt.xlabel('sample index')
  plt.ylabel('distance')
  scipy.cluster.hierarchy.dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
    labels=y_all
  )
  plt.savefig(target)

  if len(vcf_categories) > 0:
    for category in vcf_categories:
      y_cat = [s for s in y if s in vcf_categories[category]]
      if len(y_cat) > 1:
        logging.info('clustering %i samples in %s', len(y_cat), category)
        logging.debug(list(y_cat))
        X_cat = X[X.index.isin(y_cat)]
        Z = scipy.cluster.hierarchy.linkage(X_cat, method='ward', metric='euclidean')
        plt.figure(figsize=(30, 12))
        plt.title('Hierarchical Clustering Dendrogram - {}'.format(category))
        plt.xlabel('sample index')
        plt.ylabel('distance')
        scipy.cluster.hierarchy.dendrogram(
          Z,
          leaf_rotation=90.,  # rotates the x axis labels
          leaf_font_size=8.,  # font size for the x axis labels
          labels=y_cat
        )
        plt.savefig(target.replace('.png', '.{}.png'.format(category)))
      else:
        logging.info('skipping %s with %i samples', category, len(y_cat))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Cluster MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--image', default='dend.png', help='dend file')
  parser.add_argument('--categories', required=False, help='CSV of additional info about the sample')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.image, args.categories)

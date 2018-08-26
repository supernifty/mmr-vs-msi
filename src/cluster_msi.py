#!/usr/bin/env python
'''
  cluster on ms location
'''

import argparse
import collections
import csv
import logging
import math
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

import scipy.cluster.hierarchy

import sklearn.decomposition
import sklearn.manifold
import sklearn.metrics

FIGSIZE=(36,36)
FIGSIZE_TREE=(30,12)

matplotlib.rcParams.update({'font.size': 36})

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

  # distance test
  #d = 0
  #e = 0
  #for i in range(X.shape[1]):
  #  d += abs(X.iloc[0][i] - X.iloc[1][i])
  #  e += math.pow(X.iloc[0][i] - X.iloc[1][i], 2)
  #logging.debug('{} vs {}: {}, {}'.format(y_all[0], y_all[1], d, math.sqrt(e)))

  # pca
  logging.info('pca...')
  pca = sklearn.decomposition.PCA(n_components=2)
  projection = pca.fit_transform(X)
  plt.figure(figsize=FIGSIZE)
  plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
  for i, sample in enumerate(y_all):
    plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  plt.savefig(target.replace('.png', '.pca.png'))

  # also write pca as tsv
  pca_tsv = target.replace('.png', '.pca.tsv')
  with open(pca_tsv, 'w') as fh:
    fh.write('sample\tprojection_0\tprojection_1\n') 
    for i, sample in enumerate(y_all):
      fh.write('{sample}\t{projection_0:.2f}\t{projection_1:.2f}\n'.format(sample=sample, projection_0=projection[i, 0], projection_1=projection[i, 1])) 

  # mds
  logging.info('mds...')
  distances = sklearn.metrics.pairwise_distances(X, metric='cosine')
  mds = sklearn.manifold.LocallyLinearEmbedding(n_neighbors=8, n_components=2, method='modified', eigen_solver='dense')
  projection = mds.fit_transform(distances)
  plt.figure(figsize=FIGSIZE)
  plt.scatter(projection[:, 0], projection[:, 1])
  for i, sample in enumerate(y_all):
    plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  plt.savefig(target.replace('.png', '.mds.png'))

  # tsne
  logging.info('tsne...')
  tsne = sklearn.manifold.TSNE(n_components=2, init='pca', random_state=0)
  projection = tsne.fit_transform(X)
  plt.figure(figsize=FIGSIZE)
  plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
  for i, sample in enumerate(y_all):
    plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  plt.savefig(target.replace('.png', '.tsne.png'))


  #km = sklearn.cluster.KMeans(n_clusters=4, init='k-means++', max_iter=100, n_init=1, verbose=True)
  #km.fit(X)

  logging.info('hierarchy...')

  # hierarchical cluster
  Z = scipy.cluster.hierarchy.linkage(X, method='ward', metric='euclidean')

  logging.info('plotting...')
  plt.figure(figsize=FIGSIZE_TREE)
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
        plt.figure(figsize=FIGSIZE_TREE)
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

        # pca
        logging.info('pca...')
        pca = sklearn.decomposition.PCA(n_components=2)
        projection = pca.fit_transform(X_cat)
        plt.figure(figsize=FIGSIZE)
        plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
        for i, sample in enumerate(y_cat):
          plt.annotate(sample, (projection[i, 0], projection[i, 1]))
        plt.savefig(target.replace('.png', '.{}.pca.png'.format(category)))


        # also write pca as tsv
        pca_tsv = target.replace('.png', '.{}.pca.tsv'.format(category))
        with open(pca_tsv, 'w') as fh:
          fh.write('sample\tprojection_0\tprojection_1\n') 
          for i, sample in enumerate(y_cat):
            fh.write('{sample}\t{projection_0:.2f}\t{projection_1:.2f}\n'.format(sample=sample, projection_0=projection[i, 0], projection_1=projection[i, 1])) 

      else:
        logging.info('skipping %s with %i samples', category, len(y_cat))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Cluster MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--image', default='cluster.png', help='dend file')
  parser.add_argument('--categories', required=False, help='CSV of additional info about the sample')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.image, args.categories)

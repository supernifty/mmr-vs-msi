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

import scipy.cluster.hierarchy

import sklearn.decomposition
import sklearn.manifold
import sklearn.metrics

def main(fh, target, genes_fn):
  logging.info('reading cluster results from stdin...')

  genes = set()
  for line in open(genes_fn, 'r'):
    genes.add(line.strip())
  logging.info('%i genes', len(genes))

  X = []
  ys = []
  cs = []
  cats = set()
  first = True
  included_genes = set()
  for line in fh:
    fields = line.strip('\n').split('\t')
    if len(fields) < 6:
      continue
    if first:
      first = False
      for idx in range(len(fields)):
        if '|' in fields[idx]:
          gene, caller = fields[idx].split('|')
          if gene in genes:
            included_genes.add(idx)
      continue
    # Sample  Patient SampleType      MSH2Stain       Burden  A1CF|gridss ...
    ys.append(fields[0])
    cs.append(fields[2])
    cats.add(fields[2])
    X.append([fields[idx] for idx in included_genes])

  y_all = ['{} {}'.format(ys[i], cs[i]) for i in range(len(ys))]

  # pca
  logging.info('pca...')
  pca = sklearn.decomposition.PCA(n_components=2)
  projection = pca.fit_transform(X)
  plt.figure(figsize=(18, 18))
  plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
  for i, sample in enumerate(y_all):
    plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  plt.savefig(target.replace('.png', '.pca.png'))

  # mds
  #logging.info('mds...')
  #distances = sklearn.metrics.pairwise_distances(X, metric='cosine')
  #mds = sklearn.manifold.LocallyLinearEmbedding(n_neighbors=8, n_components=2, method='modified', eigen_solver='dense')
  #projection = mds.fit_transform(distances)
  #plt.figure(figsize=(18, 18))
  #plt.scatter(projection[:, 0], projection[:, 1])
  #for i, sample in enumerate(y_all):
  #  plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  #plt.savefig(target.replace('.png', '.mds.png'))

  # tsne
  #logging.info('tsne...')
  #tsne = sklearn.manifold.TSNE(n_components=2, init='pca', random_state=0)
  #projection = tsne.fit_transform(X)
  #plt.figure(figsize=(18, 18))
  #plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
  #for i, sample in enumerate(y_all):
  #  plt.annotate(sample, (projection[i, 0], projection[i, 1]))
  #plt.savefig(target.replace('.png', '.tsne.png'))

  logging.info('hierarchy...')

  # hierarchical cluster
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

  if len(cats) > 0:
    for category in cats:
      y_cat = [ys[i] for i in range(len(ys)) if cs[i] == category]
      if len(y_cat) > 1:
        logging.info('clustering %i samples in %s', len(y_cat), category)
        logging.debug(list(y_cat))
        X_cat = [X[i] for i in range(len(X)) if cs[i] == category]
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

        # pca
        logging.info('pca...')
        pca = sklearn.decomposition.PCA(n_components=2)
        projection = pca.fit_transform(X_cat)
        plt.figure(figsize=(18, 18))
        plt.scatter(projection[:, 0], projection[:, 1], alpha=0.5)
        for i, sample in enumerate(y_cat):
          plt.annotate(sample, (projection[i, 0], projection[i, 1]))
        plt.savefig(target.replace('.png', '.{}.pca.png'.format(category)))
      else:
        logging.info('skipping %s with %i samples', category, len(y_cat))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Cluster MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--image', default='dend.png', help='dend file')
  parser.add_argument('--genes', help='genes filter')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.image, args.genes)

#!/usr/bin/env python
'''
  graph proportions and genes
'''

import argparse
import collections
import csv
import logging
import sys

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
plt.rc('savefig', dpi=300)

import numpy as np

FIGSIZE=(24,24)

def rand_jitter(arr):
    stdev = .03 * (max(arr)-min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def main(output, labels, cosmic_fn):
  cosmic = set()
  if cosmic_fn is not None:
    logging.info('reading cosmic...')
    for line in open(cosmic_fn, 'r'):
      cosmic.add(line.split('\t')[0])

  logging.info('reading from stdin...')
  plt.figure(figsize=FIGSIZE)
  # Gene    Group1_Affected Group1_Unaffected       Group1_Prop     Group2_Affected Group2_Unaffected       Group2_Prop     Prop_Diff
  header = None
  points = collections.defaultdict(list)
  xs = []
  ys = []
  ls = []
  for idx, row in enumerate(csv.reader(sys.stdin, delimiter='\t')):
    if header is None:
      header = row
      continue
    point = (float(row[3]), float(row[6]))
    points[point].append(row[0])
    xs.append(point[0])
    ys.append(point[1])
    ls.append(row[0])

  xs = rand_jitter(xs)
  ys = rand_jitter(ys)

  plt.scatter(xs, ys, alpha=0.5)
  plt.title('Comparison of proportion of samples containing a mutation by gene')
  plt.xlabel('Proportion of samples affected in group 1 (MMRd)')
  plt.ylabel('Proportion of samples affected in group 2 (MMRp)')
  plt.xlim(-0.1, 1.1)
  plt.ylim(-0.1, 1.1)

  plt.annotate("", xy=(0, 0), xycoords='data', xytext=(1, 1), textcoords='data',
              arrowprops=dict(arrowstyle="-",
                              edgecolor = "red",
                              linewidth=2,
                              alpha=0.6,
                              connectionstyle="arc3,rad=0."), 
              )

  if labels:
    for i, gene in enumerate(ls):
      if gene in cosmic:
        plt.annotate(gene, (xs[i], ys[i]), color='#000030')
      else:
        plt.annotate(gene, (xs[i], ys[i]), color='#606060', alpha=0.1)

  logging.info('saving image...')
  plt.savefig(output)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot affected proportions between groups')
  parser.add_argument('--output', default='affected.png', help='output file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--labels', action='store_true', help='label genes')
  parser.add_argument('--cosmic', required=False, help='cosmic tsv')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.output, args.labels, args.cosmic)

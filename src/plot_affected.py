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

FIGSIZE=(12,12)

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
  # Gene    Group1_Affected Group1_Unaffected       Group1_Prop     Group2_Affected Group2_Unaffected       Group2_Prop     Prop_Diff   p_value
  header = None
  xs = []
  ys = []
  ls = []
  pvalues = []
  for idx, row in enumerate(csv.reader(sys.stdin, delimiter='\t')):
    if header is None:
      header = row
      continue
    point = (float(row[3]), float(row[6]))
    xs.append(point[0])
    ys.append(point[1])
    #pvalues.append(float(row[8])) # raw p value
    pvalues.append(float(row[9])) # adjusted p value
    ls.append(row[0])

  xs = rand_jitter(xs)
  ys = rand_jitter(ys)

  plt.scatter(xs, ys, c=pvalues, alpha=0.5, cmap='jet_r')
  plt.title('Proportion of tumours with exonic microsatellite INDEL by gene (MMR proficient versus MMR deficient)')
  plt.xlabel('Proportion of MMR deficient tumours with exonic microsatellite INDEL')
  plt.ylabel('Proportion of MMR proficient tumours with exonic microsatellite INDEL')
  plt.xlim(-0.1, 1.1)
  plt.ylim(-0.1, 1.1)

  cb = plt.colorbar()
  cb.ax.get_yaxis().labelpad = 15
  cb.ax.set_ylabel('FDR corrected p-value of significant difference', rotation=270)

  plt.annotate("", xy=(0, 0), xycoords='data', xytext=(1, 1), textcoords='data',
              arrowprops=dict(arrowstyle="-",
                              edgecolor = "red",
                              linewidth=2,
                              alpha=0.6,
                              connectionstyle="arc3,rad=0."), 
              )

  if labels:
    for i, gene in enumerate(ls):
      if gene.split('_')[0] in cosmic or pvalues[i] < 0.05:
        plt.annotate(gene, (xs[i], ys[i]), color='#000030')
      else:
        plt.annotate(gene, (xs[i], ys[i]), color='#606060', alpha=0.2)

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

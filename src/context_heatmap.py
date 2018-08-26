#!/usr/bin/env python
'''
  take context info and generate a heatmap
  Sample  Patient Type    All     Exon    Onco    OncoAll OncoExon        Mono    Bi      Tri     Tetra   Short   Medium  Long    Bethesda        total_gridss    total_pindel    total_varscan   rt_A    rt_AAAC ...
'''

import argparse
import logging
import math
import sys

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
plt.rc('savefig', dpi=300)

import numpy as np
import pandas as pd
import seaborn as sns

FIGSIZE_FACTOR=2

rc={'font.size': 4 * FIGSIZE_FACTOR, 'axes.labelsize': 4 * FIGSIZE_FACTOR, 'legend.fontsize': 4.0 * FIGSIZE_FACTOR, 'axes.titlesize': 4 * FIGSIZE_FACTOR, 'xtick.labelsize': 4 * FIGSIZE_FACTOR, 'ytick.labelsize': 4 * FIGSIZE_FACTOR}
sns.set(rc=rc)

import scipy.cluster.hierarchy

def rotate_repeat(r):
  '''
    normalize repeats to start from earliest base e.g. GCAA -> AAGC
  '''
  best = 0
  for idx, el in enumerate(r):
    if r[idx] < r[best]:
      best = idx

  return r[best:] + r[:best]

def main(target_image, in_fh, category, prefix='rt_', log=False, normalize_max=False, min_len=0, max_len=1000, highlight=None, normalize_custom=None, threshold_proportion=0, rotate_context=False, normalize_sample=False):
  logging.info('reading from stdin...')
  header = None
  df = None
  features = [] # indexes to read from
  input_feature_names = [] # initial feature names
  final_feature_names = [] # final feature names
  for line in in_fh:
    fields = line.strip('\n').split('\t')

    if header is None:
      header = fields
      for idx in range(len(fields)):
        if fields[idx].startswith(prefix):
          if (len(prefix) + min_len) <= len(fields[idx]) <= (len(prefix) + max_len):
            features.append(idx)
            feature_name = fields[idx][len(prefix):]
            input_feature_names.append(feature_name)
            if rotate_context:
              feature_name = rotate_repeat(feature_name)
              if feature_name not in final_feature_names:
                final_feature_names.append(feature_name)
            else:
              final_feature_names.append(feature_name)
          else:
            logging.debug('skipping feature %s due to length', fields[idx])
      #features = features[:10]
      df = pd.DataFrame(columns=['sample'] + final_feature_names)
      #df.columns.name = 'context'
      continue

    if category is None or category == fields[header.index('Type')]:
      if category is None:
        sample = '{} {} {}'.format(fields[header.index('Sample')], fields[header.index('Patient')], fields[header.index('Type')])
      else:
        sample = fields[header.index('Sample')]

      if highlight is not None and fields[header.index('Sample')] in highlight:
        sample = '***{}***'.format(sample)

      if log:
        row_features = [math.log(1 + float(fields[idx])) for idx in features]
      else:
        row_features = [float(fields[idx]) for idx in features]

      # compress into final feature names
      if rotate_context:
        final_features = [0] * len(final_feature_names)
        for idx, feature in enumerate(row_features):
          target = final_feature_names.index(rotate_repeat(input_feature_names[idx]))
          final_features[target] += feature
      else:
        final_features = row_features

      df = df.append(pd.Series(data=[sample] + final_features, index=df.columns), ignore_index=True)
    else:
      pass
      #logging.debug('skipped %s', line)
  
  df = df.convert_objects(convert_numeric=True)
  df.set_index('sample', inplace=True)

  # thresholds
  if threshold_proportion > 0:
    threshold = int(threshold_proportion * df.shape[0])
    logging.info('applying threshold of %i (%.2f of %i samples)', threshold, threshold_proportion, df.shape[0])
    for feature in final_feature_names:
      if rotate_context:
        feature = rotate_repeat(feature)
      if df[feature].sum() < threshold:
        logging.info('dropping %s with low count', feature)
        df.drop(feature, axis=1, inplace=True)

  logging.info(df.head())

  if normalize_max:
    logging.info('normalizing by variant type count...')
    df = df / df.max()
    #df = ( df - df.min() ) / ( df.max() - df.min() )
    # and remove na columns
    df.dropna(1, inplace=True)

  elif normalize_custom is not None:
    logging.info('normalizing by %s...', normalize_custom)
    updated = 0
    for line in open(normalize_custom, 'r'):
      fields = line.strip('\n').split('\t')
      if min_len <= len(fields[0]) <= max_len:
        if fields[0] in df.columns:
          df[fields[0]] /= int(fields[1])
          updated += 1
        else:
          logging.info('skipped normalizing %s', fields[0])
    logging.info('normalizing by %s: updated %i', normalize_custom, updated)

  elif normalize_sample:
    logging.info('normalizing by sample variant count...')
    df = df.div(df.sum(axis=1), axis=0)

  df = df.transpose()

  logging.info('(features, samples) = %s', df.shape)

  if df.shape[1] < 2:
    logging.warn('not enough samples. exiting')
    sys.exit(0)
  if df.shape[0] < 1:
    logging.warn('not enough features. exiting')
    sys.exit(0)

  logging.info('plotting...')
  plot = sns.clustermap(df, figsize=(max(FIGSIZE_FACTOR * 8, FIGSIZE_FACTOR * df.shape[1] / 10), max(FIGSIZE_FACTOR * 8, FIGSIZE_FACTOR * df.shape[0] / 10)))
  plot.savefig(target_image)
    
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Heatmap')
  parser.add_argument('--image', required=True, help='png output')
  parser.add_argument('--category', required=False, help='filter on type')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--normalize_max', action='store_true', help='normalize using max count of each feature')
  parser.add_argument('--normalize_sample', action='store_true', help='normalize using max count of each sample')
  parser.add_argument('--normalize_custom', required=False, help='normalize using a provided list')
  parser.add_argument('--log', action='store_true', help='apply log to counts')
  parser.add_argument('--prefix', required=False, default='rt_', help='feature prefix')
  parser.add_argument('--min_len', required=False, default=0, type=int, help='feature prefix')
  parser.add_argument('--max_len', required=False, default=1000, type=int, help='feature prefix')
  parser.add_argument('--highlight', required=False, nargs='+', help='samples to highlight')
  parser.add_argument('--threshold', required=False, default=0, type=float, help='minimum mutation count as a proportion of sample count')
  parser.add_argument('--rotate_context', action='store_true', help='rotate contexts')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.image, sys.stdin, args.category, args.prefix, args.log, args.normalize_max, args.min_len, args.max_len, args.highlight, args.normalize_custom, args.threshold, args.rotate_context, args.normalize_sample)

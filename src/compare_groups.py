#!/usr/bin/env python
'''
  Find common and unique genes in groups
'''

import argparse
import collections
import csv
import logging
import sys

import scipy.stats
import statsmodels.stats.multitest

def main(group1_list, group2_list):
  logging.info('reading from stdin...')
  # input is of the form:
  # sample gene...
  group1 = set(group1_list)
  group2 = set(group2_list)

  header = None
  count = {}
  for idx, row in enumerate(csv.reader(sys.stdin, delimiter='\t')):
    if header is None:
      header = row
      for gene in header[1:]:
        count[gene] = collections.defaultdict(int)
      continue

    if row[0] in group1:
      group = '1'
    elif row[0] in group2:
      group = '2'
    else:
      logging.info('skipping sample %s', row[0])
      continue

    for gene_col in range(1, len(header)):
      gene = header[gene_col]
      if row[gene_col] == '0': # unaffected
        count[gene]['{}_unaffected'.format(group)] += 1
      else:
        count[gene]['{}_affected'.format(group)] += 1

  # now write results 
  sys.stdout.write('Gene\tGroup1_Affected\tGroup1_Unaffected\tGroup1_Prop\tGroup2_Affected\tGroup2_Unaffected\tGroup2_Prop\tProp_Diff\tp_value\tp_adjusted\n')

  result = []
  pvalues = []
  for gene in sorted(header[1:]):
    if count[gene]['1_affected'] + count[gene]['1_unaffected'] > 0:
      g1_prop = count[gene]['1_affected'] / (count[gene]['1_affected'] + count[gene]['1_unaffected'])
    else:
      logging.warn('Gene {} has no samples in group 1', gene)
      continue

    if count[gene]['2_affected'] + count[gene]['2_unaffected'] > 0:
      g2_prop = count[gene]['2_affected'] / (count[gene]['2_affected'] + count[gene]['2_unaffected'])
    else:
      logging.warn('Gene {} has no samples in group 2', gene)
      continue

    contingency = [[count[gene]['1_affected'], count[gene]['2_affected']], [count[gene]['1_unaffected'], count[gene]['2_unaffected']]]
    if contingency[0] == [0, 0] or contingency[1] == [0, 0]: # both groups have no affected, or both groups have no unaffected
      p_value = 1
    else:
      p_value = scipy.stats.chi2_contingency(contingency)[1]

    pvalues.append(p_value)
    result.append('{}\t{}\t{}\t{:.3f}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.7f}'.format(gene, count[gene]['1_affected'], count[gene]['1_unaffected'], g1_prop, count[gene]['2_affected'], count[gene]['2_unaffected'], g2_prop, g1_prop - g2_prop, p_value))

  logging.info('calculating adjusted p-values...')
  adjusteds = statsmodels.stats.multitest.multipletests(pvalues, method='fdr_bh')[1]
  
  logging.info('writing results...')
  for output, adjusted in zip(result, adjusteds):
    sys.stdout.write('{}\t{:.7f}\n'.format(output, adjusted))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--group1', required=True, nargs='+', help='samples in group 1')
  parser.add_argument('--group2', required=True, nargs='+', help='samples in group 2')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.group1, args.group2)

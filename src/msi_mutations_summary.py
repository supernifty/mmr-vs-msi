#!/usr/bin/env python
'''
  takes msi details and counts prevalence by indel length and ms type
'''

import argparse
import collections
import logging
import sys

import scipy
import scipy.stats

MIN_COUNT = 10
LONG_INDEL = 5
Z_SIG = 2

def main(fh_in, out, proportions, statistics):
  logging.info('starting...')

  # Chrom   Begin   End     Annotation      Len     Samples CMHS1
  first = True
  keys = collections.defaultdict(int)
  header = None
  samples = None
  summary = {}
  for line in sys.stdin:
    fields = line.strip('\n').split('\t')
    if fields[0] == 'TOTAL':
      continue
    if first:
      header = fields
      samples = fields[6:]
      for sample in samples:
        summary[sample] = collections.defaultdict(int)
      first = False
      continue
    
    # extract annotations of interest
    annotations = {}
    for kv in fields[3].split(';'):
      k, v = kv.split('=')
      annotations[k] = v

    # each sample
    for idx, sample in enumerate(samples):
      indel_length = fields[6 + idx]
      if indel_length != '0':
        if int(indel_length) >= LONG_INDEL:
          indel_length = '+'
        elif int(indel_length) <= -LONG_INDEL:
          indel_length = '-'
        key = '{}{}'.format(annotations['repeat'], indel_length)
        keys[key] += 1
        summary[sample][key] += 1

  filtered_keys = [key for key in keys if keys[key] >= MIN_COUNT]

  logging.info('%i keys %i filtered keys', len(keys), len(filtered_keys))

  # write info for each sample
  out.write('Sample\t{}\n'.format('\t'.join(sorted(filtered_keys))))
  for sample in samples:
    if proportions:
      total = sum([summary[sample][key] for key in keys])
      out.write('{}\t{}\n'.format(sample, '\t'.join(['{:.2f}'.format(summary[sample][key] / total) for key in sorted(filtered_keys)])))
    else: 
      out.write('{}\t{}\n'.format(sample, '\t'.join([str(summary[sample][key]) for key in sorted(filtered_keys)])))
  
  # calculate stats for each key
  if statistics is not None:
    logging.info('statistics...')
    with open(statistics, 'w') as stats_fh:
      for key in sorted(filtered_keys):
        values = []
        for sample in samples:
          if proportions:
            total = sum([summary[sample][key] for key in keys])
            values.append(summary[sample][key] / total)
          else:
            values.append(summary[sample][key])
  
        zs = scipy.stats.zscore(values) # only valid if n is large enough
        zp = p_values = (1 - scipy.special.ndtr(zs)) * len(keys) * len(samples) # bonferroni
        #stats_fh.write("{}\t{}\n".format(key, ','.join(['{}/{:.2f}'.format(samples[x], zs[x]) for x in range(len(samples)) if zs[x] > Z_SIG or zs[x] < -Z_SIG])))
        stats_fh.write("{}\t{}\n".format(key, ','.join(['{}/{:.4f}'.format(samples[x], zp[x]) for x in range(len(samples)) if zp[x] < 0.05])))
    
    logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Summarise MSI characteristics')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--proportions', action='store_true', help='use proportions')
  parser.add_argument('--statistics', required=False, help='write stats')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout, args.proportions, args.statistics)

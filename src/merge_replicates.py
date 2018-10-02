#!/usr/bin/env python
'''
  merge replicates from affected_genes using a majority rules
  produces:
    sample gene(s)
    name   count...
  usage:
    src/merge_replicates.py --replicates A1,A2,A3 B1,B2,B3 < msi.affected_genes.tsv >msi.affected_genes.2.tsv 2>log
'''

import argparse
import collections
import csv
import logging
import sys

def main(replicates):
  replicate_target = {}
  if replicates is not None:
    for replicate in replicates:
      samples = replicate.split(',')
      target = samples[0]
      for sample in samples:
        replicate_target[sample] = target
    logging.info('%i replicate targets', len(replicate_target))
      
  logging.info('reading from stdin...')
  header = None
  results = collections.defaultdict(list) # sample -> list of results
  for idx, row in enumerate(csv.reader(sys.stdin, delimiter='\t')): # each msi location
    if header is None:
      header = row
      continue
    sample = row[0]
    if sample in replicate_target:
      sample = replicate_target[sample]
    results[sample].append(row[1:])

  # write results
  logging.info('writing results...')
  sys.stdout.write('{}\n'.format('\t'.join(header)))
  for sample in sorted(results):
    stats = {'all_affected': 0, 'all_unaffected': 0, 'discordant': 0}
    result = results[sample]
    if len(result) == 1:
      sys.stdout.write('{}\t{}\n'.format(sample, '\t'.join(result[0])))
    else: # merge results
      merged = []
      for items in zip(*result):
        unaffected = sum([item == '0' for item in items])
        if unaffected > len(items) / 2:
          merged.append('0')
        else:
          merged.append('1')
        # stats
        if unaffected == 0:
          stats['all_affected'] += 1
        elif unaffected == len(items):
          stats['all_unaffected'] += 1
        else:
          stats['discordant'] += 1
      total = sum([stats[x] for x in stats])
      logging.info('replicate stats for %s: %s, %s', sample, stats, ' '.join(['{}: {:.1f}%'.format(stat, stats[stat] * 100 / total) for stat in stats]))
      sys.stdout.write('{}\t{}\n'.format(sample, '\t'.join(merged)))

  logging.info('done. wrote %i', len(results))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--replicates', nargs='*', help='list of comma separated replicates')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.replicates)


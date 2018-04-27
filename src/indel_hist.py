#!/usr/bin/env python
'''
  given an input vcf, generate a histogram of indel lengths
'''

import argparse
import collections
import logging
import sys

import cyvcf2

def main(out):
  logging.info('reading from stdin...')
  stats_all = collections.defaultdict(int)
  stats_exon = collections.defaultdict(int)
  stats_onco = collections.defaultdict(int)
  stats_exon_onco = collections.defaultdict(int)
  for variant in cyvcf2.VCF('-'):
    net = len(variant.ALT[0]) - len(variant.REF)
    stats_all[net] += 1
    if variant.INFO.get('msi_exon') is not None:
      stats_exon[net] += 1
      if variant.INFO.get('msi_oncogene') is not None:
        stats_exon_onco[net] += 1
    if variant.INFO.get('msi_oncogene') is not None:
      stats_onco[net] += 1

  out.write('Change\tAll\tExon\tOnco\tExonOnco\n')
  for stat in sorted(stats_all.keys()):
    out.write('{change}\t{total}\t{exon}\t{onco}\t{exon_onco}\n'.format(change=stat, total=stats_all[stat], exon=stats_exon[stat], onco=stats_onco[stat], exon_onco=stats_exon_onco[stat]))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Indel histogram')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(sys.stdout)

#!/usr/bin/env python
'''
  visualize uniqueness
'''

import argparse
import collections
import logging
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
  commonness = []
  sample_singletons = collections.defaultdict(int)
  for line in fh:
    fields = line.strip('\n').split('\t')
    if fields[0] in ('Chrom', 'TOTAL'):
      if fields[0] == 'Chrom':
        header = fields
      continue
    commonness.append(int(fields[5]))

    if int(fields[5]) == 1: # a singleton
      for idx in range(6, len(fields)):
        if fields[idx] != '0':
          sample = header[idx]
          sample_singletons[sample] += 1

  # commonness: histogram of how common each loci is
  plt.figure(figsize=(24, 18))
  plt.title('How common is each variant? How many samples is each variant found in?')
  plt.xlabel('Number of Samples')
  plt.ylabel('Site Count')
  plt.hist(commonness, bins=range(1, max(commonness) + 1))
  plt.savefig(target.replace('.png', '.loci.png'))

  # singletons: how many unique sites each sample has
  plt.figure(figsize=(24, 18))
  plt.title('How many singletons does each sample have?')
  plt.xlabel('Singletons')
  plt.ylabel('Number of Samples')
  shared = [sample_singletons[x] for x in sample_singletons]
  plt.hist(shared, bins=range(1, max(shared) + 1))
  plt.savefig(target.replace('.png', '.samples.png'))

  sys.stdout.write('Sample\tSingletons\n')
  for sample in sorted(sample_singletons.keys()):
    sys.stdout.write('{}\t{}\n'.format(sample, sample_singletons[sample]))

  if len(vcf_categories) > 0:
    for category in vcf_categories:
      pass

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Cluster MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--image', default='unique.png', help='hist file')
  parser.add_argument('--categories', required=False, help='CSV of additional info about the sample')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, args.image, args.categories)

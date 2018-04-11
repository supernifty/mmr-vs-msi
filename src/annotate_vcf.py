#!/usr/bin/env python

'''
  Intersects a vcf file with a bed and annotates it with bed values

  Usage:
    $0 bed < vcf > annotated_filtered_vcf
'''

import logging
import sys

import cyvcf2
import intervaltree

def main(bed):
  logging.info('parsing {}...'.format(bed))
  tree = {}
  annotations = set()
  for idx, line in enumerate(open(bed, 'r')):
    chrom, start, finish, annotation = line.strip('\n').split('\t')
    if chrom not in tree:
      tree[chrom] = intervaltree.IntervalTree()
    tree[chrom][int(start):int(finish)] = annotation
    for namevalue in annotation.split(';'):
      name, _ = namevalue.split('=')
      annotations.add(name)
    if idx % 100000 == 0:
      logging.debug('parsing {}: {} lines parsed'.format(bed, idx))
  logging.info('parsing {}: done. annotations: {}'.format(bed, ', '.join(list(annotations))))

  logging.info('parsing vcf from stdin...'.format(bed))
  vcf_in = cyvcf2.VCF('-')
  for annotation in annotations:
    vcf_in.add_info_to_header({'ID': 'msi_{}'.format(annotation), 'Description': 'details of msi {}'.format(annotation), 'Type':'Character', 'Number': '1'})

  accept = reject = 0
  sys.stdout.write(vcf_in.raw_header)
  for variant in vcf_in:
    # does it overlap the bed?
    if variant.CHROM in tree:
      overlap = tree[variant.CHROM].search(variant.POS)
      if len(overlap) == 0:
        reject += 1
      else:
        accept += 1
        for namevalue in list(overlap)[0].data.split(';'):
          name, value = namevalue.split('=')
          variant.INFO["msi_{}".format(name)] = value
        sys.stdout.write(str(variant))
    else:
      reject += 1
  logging.info('included %i variants. rejected %i variants.', accept, reject)

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1])

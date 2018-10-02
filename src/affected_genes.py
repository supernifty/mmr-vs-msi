#!/usr/bin/env python
'''
  simplify msi.cluster.tumour to counts by gene
  produces:
    sample gene(s)
    name   count...
  usage:
    src/affected_genes.py < msi.cluster.tumurs >msi.affected_genes.tsv 2>log
'''

import argparse
import collections
import csv
import logging
import sys

def find_genes(annotations, annotate):
  result = set()
  # determine the gene from annotation repeat=TTCT;length=8;exon=DDX11L1_exon_0_0_chr1_11874_f[11873:12227],LOC102725121_exon_0_0_chr1_11869_f[11868:12227];exon=DDX11L1_exon_0_3_11874_f[11873:12227];ncRNA=1[11873:12227]
  additional = set()
  for ann in annotations.split(';'):
    if '=' in ann:
      name, values = ann.split('=', 1)
      if name == 'exon':
        for value in values.split(','):
          gene, _ = value.split('_', 1)
          result.add(gene)
      elif annotate is not None and name in annotate:
        additional.add(name)
    else:
      logging.debug('skipped %s', ann)
  if len(additional) == 0:
    return result

  annotated_result = set()
  for gene in result:
    annotated_result.add('{}_{}'.format(gene, '_'.join(list(additional))))
  return annotated_result

def main(annotate):
  logging.info('reading from stdin...')

  counts = {}
  genes_seen = set()
  header = None
  # Chrom   Begin   End     Annotation      Len     Samples 0151010101_T    0151100011_T...
  for idx, row in enumerate(csv.reader(sys.stdin, delimiter='\t')): # each msi location
    if header is None:
      header = row
      for sample in header[6:]:
        counts[sample] = collections.defaultdict(int)
      continue

    # genes associated with this msi location
    genes = find_genes(row[3], annotate)
    
    # now update the counts of all samples
    for sample_col in range(6, len(header)):
      sample = header[sample_col]
      for gene in genes: # all genes associated with this loci
        if row[sample_col] != '0':
          counts[sample][gene] += 1
          genes_seen.add(gene)

    if idx % 10000 == 0:
      logging.debug('processed %i lines. %i genes seen...', idx, len(genes_seen))

  # now write results
  sorted_genes = sorted(genes_seen)
  sys.stdout.write('Sample\t{}\n'.format('\t'.join(sorted_genes)))
  for sample in sorted(counts.keys()):
    line = '\t'.join([str(counts[sample][gene]) for gene in sorted_genes])
    sys.stdout.write('{}\t{}\n'.format(sample, line))

  logging.info('done. %i genes seen.', len(genes_seen))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--annotate', nargs='+', help='additional annotations to include from input file')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.annotate)

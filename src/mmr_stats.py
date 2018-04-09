#!/usr/bin/env python

import collections
import os
import sys

import cyvcf2

def to_type(is_normal):
  if is_normal == 'Y':
    return 'Normal'
  else:
    return 'Tumour'

SAMPLE_METADATA='cfg/sample-metadata.csv'

VEP="Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|PICK".split('|')
gene_index = VEP.index('SYMBOL')

def main(files):
  stats = collections.defaultdict(int)
  of_interest = collections.defaultdict(list)
  gene_count = collections.defaultdict(int)

  for f in files:
    count = 0
    for variant in cyvcf2.VCF(f):
      sample = os.path.basename(f).split('.')[0]
      caller = os.path.basename(f).split('.')[2]
      gene = variant.INFO.get('CSQ').split('|')[gene_index]
      of_interest[sample].append(gene)
      gene_count['{}|{}'.format(sample, gene)] += 1
      count += 1
    stats[count] += 1 # sample

  metadata = {}
  for line in open(SAMPLE_METADATA, 'r'):
    fields = line.strip('\n').split(',')
    metadata[fields[0]] = (fields[1], fields[4])


  sys.stdout.write('VariantCount\tNumberOfSamples\n{}'.format('\n'.join([ '{}\t{}'.format(count, stats[count]) for count in sorted(stats)])))

  # TODO woah
  result = [(sample, metadata[sample][0], to_type(metadata[sample][1]), ', '.join(['{}:{}'.format(gene, gene_count['{}|{}'.format(sample, gene)]) for gene in of_interest[sample]])) for sample in of_interest]

  sys.stdout.write('\n\nCandidates ({}):\nSample\tPatient\tSampleType\tGenes\n{}\n'.format(
    len(of_interest), 
    '\n'.join(['{}\t{}\t{}\t{}'.format(x[0], x[1], x[2], x[3]) for x in sorted(result, key=lambda x: x[1])]))
  )

if __name__ == '__main__':
  main(sys.argv[1:])

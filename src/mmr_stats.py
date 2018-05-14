#!/usr/bin/env python
'''
  collects gene counts from input vcfs and bed files
'''

import argparse
import collections
import logging
import os
import sys

import cyvcf2

SAMPLE_METADATA='cfg/sample-metadata.csv'
HTSDB_METADATA='cfg/htsdb.csv'
IHC_METADATA='cfg/samples.ihc'

VEP="Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|PICK".split('|')
gene_index = VEP.index('SYMBOL')

def to_type(metadata, htsdb):
  if metadata[4] == 'Y':
    return 'Normal'
  else:
    return htsdb[3]

def main(show_all, files):
  stats = collections.defaultdict(int)
  all_samples = set()
  of_interest = collections.defaultdict(set)
  gene_count = collections.defaultdict(int)
  sample_count = collections.defaultdict(int)
  gene_callers = set()

  for f in files:
    logging.debug('evaluating %s', f)
    sample = os.path.basename(f).split('.')[0]
    all_samples.add(sample)
    caller = os.path.basename(f).split('.')[2]
    if f.endswith('.vcf'):
      for variant in cyvcf2.VCF(f):
        gene = variant.INFO.get('CSQ').split('|')[gene_index]
        of_interest[sample].add(gene)
        gene_count['{}|{}|{}'.format(sample, gene, caller)] += 1
        gene_callers.add('{}|{}'.format(gene, caller))
        sample_count[sample] += 1
    elif f.endswith('.bed'):
      for line in open(f, 'r'):
        fields = line.strip('\n').split('\t')
        #2       47629205        47711367        MSH2    2       42995001        49739000        2       -0.699884539104366
        gene = fields[3]
        of_interest[sample].add(gene)
        gene_count['{}|{}|{}'.format(sample, gene, caller)] += 1
        gene_callers.add('{}|{}'.format(gene, caller))
        sample_count[sample] += 1
    else:
      logging.error('Unrecognized extension for %s. Should be vcf or bed.', f)

  for sample in sample_count:
    stats[sample_count[sample]] += 1

  metadata = {}
  for line in open(SAMPLE_METADATA, 'r'):
    fields = line.strip('\n').split(',')
    metadata[fields[0]] = fields

  htsdb = {}
  for line in open(HTSDB_METADATA, 'r'):
    fields = [field.strip('"') for field in line.strip('\n').split(',')]
    htsdb[fields[0]] = fields

  stain = collections.defaultdict(str)
  for line in open(IHC_METADATA, 'r'):
    fields = line.strip('\n').split('\t')
    stain[fields[1]] = fields[2] # sample = status

  # overall summary
  sys.stdout.write('VariantCount\tNumberOfSamples\n{}'.format('\n'.join([ '{}\t{}'.format(count, stats[count]) for count in sorted(stats)])))

  # each sample
  sys.stdout.write('\n\nCandidates ({}):\nSample\tPatient\tSampleType\tMSH2Stain\tBurden\t{}\n'.format(len(of_interest), '\t'.join([gc for gc in sorted(gene_callers)])))
  if show_all:
    samples = sorted(all_samples)
  else:
    samples = sorted(of_interest)

  for sample in samples:
    genes = [gene_count['{}|{}'.format(sample, gc)] for gc in sorted(gene_callers)]
    burden = sum(genes)
    sys.stdout.write('{sample}\t{patient}\t{sample_type}\t{stain}\t{burden}\t{genes}\n'.format(sample=sample, patient=metadata[sample][1], sample_type=to_type(metadata[sample], htsdb[sample]), stain=stain[sample], burden=burden, genes='\t'.join([str(g) for g in genes])))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--all', action='store_true', help='show all sample details')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('inputs', nargs='+', help='vcf files and bed files to extract reads')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.all, args.inputs)

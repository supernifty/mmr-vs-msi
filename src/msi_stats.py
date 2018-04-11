#!/usr/bin/env python

import collections
import operator
import os
import sys

import cyvcf2

def main(mmr_result, vcfs):
  # pull in mmr result
  #Sample  Patient SampleType      Genes
  #CMHS1   CMHP1   Tumour  PMS1:1, POLE:1
  mmr = {}
  for line in open(mmr_result, 'r'):
    fields = line.strip('\n').split('\t')
    if len(fields) >= 4:
      mmr[fields[0]] = fields

  # evaluatae msi result
  samples = collections.defaultdict(int)
  for filename in vcfs:
    # for now just count indels
    sample = os.path.basename(filename).split('.')[0] # CMHS1
    count = 0
    for variant in cyvcf2.VCF(filename):
      count += 1
    samples[sample] += count # total number of msi events

  sys.stdout.write('Sample\tMSI\tPatient\tSampleType\t{}\n'.format('\t'.join(mmr['Sample'][3:])))
  for sample, count in sorted(samples.items(), key=operator.itemgetter(1)): # each sample and its overall count
    if sample in mmr:
      patient = mmr[sample][1]
      sample_type = mmr[sample][2]
      genes = mmr[sample][3:]
    else:
      sample_type = genes = patient = ''
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(sample, samples[sample], patient, sample_type, '\t'.join(genes)))
  
if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2:])

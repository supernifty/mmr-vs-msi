#!/usr/bin/env python

import logging
import os
import sys

import cyvcf2

def main(mmr_result, vcfs):
  logging.info('{} vcfs to process...'.format(len(vcfs)))
  # pull in mmr result
  #Sample  Patient SampleType      Genes
  #CMHS1   CMHP1   Tumour  PMS1:1, POLE:1
  mmr = {}
  for line in open(mmr_result, 'r'):
    fields = line.strip('\n').split('\t')
    if len(fields) >= 4:
      mmr[fields[0]] = fields

  # evaluate msi result
  callers = set()
  counts = {}
  for idx, filename in enumerate(vcfs):
    logging.info('processing file {} of {}: {}...'.format(idx + 1, len(vcfs), filename))
    # for now just count variants
    sample = os.path.basename(filename).split('.')[0] # e.g. CMHS1
    caller = os.path.basename(filename).split('.')[2] # e.g. gridss
    callers.add(caller)
    if sample not in counts:
      counts[sample] = {'all': 0, 'exon': 0, 'oncogene': 0, 'repeat_1': 0, 'repeat_2': 0, 'repeat_3': 0, 'repeat_4': 0}
    if caller not in counts[sample]:
      counts[sample][caller] = 0

    variant_count = 0
    for variant_count, variant in enumerate(cyvcf2.VCF(filename)):
      counts[sample]['all'] += 1
      counts[sample][caller] += 1
      # check other variant categories from INFO msi_repeat=T;msi_length=15
      if variant.INFO.get('msi_exon') is not None:
        counts[sample]['exon'] += 1
      if variant.INFO.get('msi_oncogene') is not None:
        counts[sample]['oncogene'] += 1
      repeat_length = len(variant.INFO['msi_repeat'])
      counts[sample]['repeat_{}'.format(repeat_length)] += 1
    logging.info('{} variants processed.'.format(variant_count + 1))

  sys.stdout.write('Sample\tPatient\tType\tAll\tExon\tOnco\tMono\tBi\tTri\tTetra\t{}\t{}\n'.format('\t'.join(sorted(list(callers))), '\t'.join(mmr['Sample'][3:])))
  for sample in sorted(counts, key=lambda k: counts[k]['all']): # each sample and its overall count
    if sample in mmr:
      patient = mmr[sample][1]
      sample_type = mmr[sample][2]
      genes = mmr[sample][3:]
    else:
      sample_type = genes = patient = ''
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, patient, sample_type, counts[sample]['all'], counts[sample]['exon'], counts[sample]['oncogene'], counts[sample]['repeat_1'], counts[sample]['repeat_2'], counts[sample]['repeat_3'], counts[sample]['repeat_4'], '\t'.join(str(counts[sample][caller]) for caller in sorted(list(callers))), '\t'.join(genes)))

  logging.info('done')
  
if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1], sys.argv[2:])

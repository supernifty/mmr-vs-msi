#!/usr/bin/env python
'''
  counts variants in provided vcf files
'''

import logging
import os
import sys

import cyvcf2

def main(mmr_result, count_result, vcfs):
  logging.info('{} vcfs to process...'.format(len(vcfs)))
  # pull in mmr result
  #Sample  Patient SampleType      Genes
  #CMHS1   CMHP1   Tumour  PMS1:1, POLE:1
  mmr = {}
  for line in open(mmr_result, 'r'):
    fields = line.strip('\n').split('\t')
    if len(fields) >= 4:
      mmr[fields[0]] = fields

  # pull in mutation result
  # this is unfiltered count of variants across the whole genome
  #Sample  Patient SampleType      Genes
  #CMHS1   CMHP1   Tumour  PMS1:1, POLE:1
  totals = {}
  total_callers = []
  first = True
  for line in open(count_result, 'r'):
    if first:
      total_callers = line.strip('\n').split('\t')[1:]
      first = False
    else:
      fields = line.strip('\n').split('\t')
      if len(fields) >= 4:
        totals[fields[0]] = fields[1:]

  # evaluate msi result
  callers = set()
  counts = {}
  for idx, filename in enumerate(vcfs):
    logging.info('processing file {} of {}: {}...'.format(idx + 1, len(vcfs), filename))
    # for now just count variants
    sample = os.path.basename(filename).split('.')[0] # e.g. CMHS1
    caller = os.path.basename(filename).split('.')[2] # e.g. gridss
    callers.add(caller)
    # all the things we want to count
    if sample not in counts:
      counts[sample] = {'all': 0, 'exon': 0, 'all_oncogene': 0, 'oncogene': 0, 'oncogene_exon': 0, 'repeat_1': 0, 'repeat_2': 0, 'repeat_3': 0, 'repeat_4': 0, 'short': 0, 'medium': 0, 'long': 0, 'bethesda': 0}
    if caller not in counts[sample]:
      counts[sample][caller] = 0

    variant_count = 0
    for variant_count, variant in enumerate(cyvcf2.VCF(filename)):
      counts[sample]['all'] += 1 # all variants for this sample across all callers
      counts[sample][caller] += 1 # all variants for this sample and caller -> total_{caller}
      # check other variant categories from INFO msi_repeat=T;msi_length=15
      if variant.INFO.get('msi_exon') is not None:
        counts[sample]['exon'] += 1
      if variant.INFO.get('msi_bethesda') is not None:
        counts[sample]['bethesda'] += 1
      if variant.INFO.get('msi_oncogene') is not None:
        counts[sample]['oncogene'] += 1
        if variant.INFO.get('msi_exon') is not None:
          counts[sample]['oncogene_exon'] += 1
      if variant.INFO.get('msi_all_oncogene') is not None:
        counts[sample]['all_oncogene'] += 1
      repeat_length = len(variant.INFO['msi_repeat'])
      counts[sample]['repeat_{}'.format(repeat_length)] += 1
      ms_length = int(variant.INFO['msi_length'])
      if ms_length < 20:
        counts[sample]['short'] += 1
      elif ms_length < 100:
        counts[sample]['medium'] += 1
      else:
        counts[sample]['long'] += 1
    logging.info('{} variants processed.'.format(variant_count + 1))

  sys.stdout.write('Sample\tPatient\tType\tAll\tExon\tOnco\tOncoAll\tOncoExon\tMono\tBi\tTri\tTetra\tShort\tMedium\tLong\tBethesda\t{}\t{}\t{}\n'.format('\t'.join([
    'total_{}'.format(x) for x in total_callers]), # 
    '\t'.join(sorted(list(callers))), 
    '\t'.join(mmr['Sample'][3:])))
  for sample in sorted(counts, key=lambda k: counts[k]['all']): # each sample and its overall count
    if sample in mmr:
      patient = mmr[sample][1]
      sample_type = mmr[sample][2]
      genes = mmr[sample][3:]
    else:
      sample_type = genes = patient = ''
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, patient, sample_type, counts[sample]['all'], counts[sample]['exon'], counts[sample]['all_oncogene'], counts[sample]['oncogene'], counts[sample]['oncogene_exon'], counts[sample]['repeat_1'], counts[sample]['repeat_2'], counts[sample]['repeat_3'], counts[sample]['repeat_4'], counts[sample]['short'], counts[sample]['medium'], counts[sample]['long'], counts[sample]['bethesda'], 
      '\t'.join([str(totals[sample][x]) for x in range(len(total_callers))]), 
      '\t'.join([str(counts[sample][caller]) for caller in sorted(list(callers))]), 
      '\t'.join(genes)))

  logging.info('done')
  
if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1], sys.argv[2], sys.argv[3:])

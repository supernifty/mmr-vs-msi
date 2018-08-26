#!/usr/bin/env python
'''
  counts variants in provided vcf files and merges with mmr results
'''

import argparse
import logging
import os
import sys

import cyvcf2

def rotate_repeat(r):
  '''
    normalize repeats to start from earliest base e.g. GCAA -> AAGC
  '''
  best = 0
  for idx, el in enumerate(r):
    if r[idx] < r[best]:
      best = idx

  return r[best:] + r[:best]

def normalize_complement(r):
  rc = reverse_complement(r)
  if r <= rc:
    return r
  else:
    return rc

RC = {'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n', 'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
def reverse_complement(repeat):
  return ''.join([RC[x] for x in repeat][::-1])
  
def transcribe_strand(repeat, exon, stats):
  # ZKSCAN3_exon_5_0_chr6_28333203_f[28333202:28336954],ZKSCAN3...
  fr = exon.split(',')[0].split('_')[-1][0]
  if fr == 'r':
    stats['strand_r'] += 1
    return reverse_complement(repeat)
  elif fr == 'f':
    stats['strand_f'] += 1
    return repeat
  else:
    logging.warn('unable to determine strand for exon: %s', exon)
    return repeat

def main(mmr_result, count_result, vcfs, rotate_context, transcribed_strand, complement):
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
  repeat_types = set()
  indel_lens = set()
  repeat_indels = set()
  for idx, filename in enumerate(vcfs): # each vcf
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

    stats = {'variant_count': 0, 'skipped': 0, 'strand_f': 0, 'strand_r': 0}

    locations_seen = set()
    for stats['variant_count'], variant in enumerate(cyvcf2.VCF(filename)): # each variant
      indel_len = len(variant.ALT[0]) - len(variant.REF)
      if indel_len == 0: # not an indel
        continue

      location = '{}/{}'.format(variant.CHROM, variant.POS)
      if location in locations_seen:
        stats['skipped'] += 1
        continue # don't double count variants at same location
      locations_seen.add(location)
      counts[sample]['all'] += 1 # all variants for this sample across all callers
      counts[sample][caller] += 1 # all variants for this sample and caller -> total_{caller}
      indel_lens.add(indel_len)
      il = 'il_{}'.format(indel_len)
      if il not in counts[sample]:
        counts[sample][il] = 0
      counts[sample][il] += 1

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

      # complement
      if complement:
        repeat_type = normalize_complement(variant.INFO['msi_repeat']) 
      else:
        repeat_type = variant.INFO['msi_repeat']

      # repeat type
      if rotate_context:
        repeat_type = rotate_repeat(repeat_type)

      # transcribed strand
      if transcribed_strand and variant.INFO.get('msi_exon') is not None:
        repeat_type = transcribe_strand(repeat_type, variant.INFO.get('msi_exon'), stats)

      repeat_types.add(repeat_type)
      rt = 'repeat_type_{}'.format(repeat_type)
      if rt not in counts[sample]:
        counts[sample][rt] = 0
      counts[sample][rt] += 1

      # repeat type and indel len
      repeat_indel = '{}{}'.format(repeat_type, indel_len)
      repeat_indels.add(repeat_indel)
      rtil = 'rtil_{}'.format(repeat_indel)
      if rtil not in counts[sample]:
        counts[sample][rtil] = 0
      counts[sample][rtil] += 1

    logging.info(', '.join(['{}: {}'.format(k, stats[k]) for k in stats]))

  sys.stdout.write('Sample\tPatient\tType\tAll\tExon\tOnco\tOncoAll\tOncoExon\tMono\tBi\tTri\tTetra\tShort\tMedium\tLong\tBethesda\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
    '\t'.join(['total_{}'.format(x) for x in total_callers]), # 
    '\t'.join(['rt_{}'.format(x) for x in sorted(list(repeat_types))]),
    '\t'.join(['il_{}'.format(x) for x in sorted(list(indel_lens))]),
    '\t'.join(['rtil_{}'.format(x) for x in sorted(list(repeat_indels))]),
    '\t'.join(sorted(list(callers))), 
    '\t'.join(mmr['Sample'][3:])))
  for sample in sorted(counts, key=lambda k: counts[k]['all']): # each sample and its overall count
    if sample in mmr:
      patient = mmr[sample][1]
      sample_type = mmr[sample][2]
      genes = mmr[sample][3:]
    else:
      sample_type = genes = patient = ''
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sample, patient, sample_type, counts[sample]['all'], counts[sample]['exon'], counts[sample]['all_oncogene'], counts[sample]['oncogene'], counts[sample]['oncogene_exon'], counts[sample]['repeat_1'], counts[sample]['repeat_2'], counts[sample]['repeat_3'], counts[sample]['repeat_4'], counts[sample]['short'], counts[sample]['medium'], counts[sample]['long'], counts[sample]['bethesda'], 
      '\t'.join([str(totals[sample][x]) for x in range(len(total_callers))]), 
      '\t'.join([str(counts[sample]['repeat_type_{}'.format(x)]) if 'repeat_type_{}'.format(x) in counts[sample] else '0' for x in sorted(list(repeat_types))]),
      '\t'.join([str(counts[sample]['il_{}'.format(x)]) if 'il_{}'.format(x) in counts[sample] else '0' for x in sorted(list(indel_lens))]),
      '\t'.join([str(counts[sample]['rtil_{}'.format(x)]) if 'rtil_{}'.format(x) in counts[sample] else '0' for x in sorted(list(repeat_indels))]),
      '\t'.join([str(counts[sample][caller]) for caller in sorted(list(callers))]), 
      '\t'.join(genes)))

  logging.info('done')
  
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--mmr', required='true', help='mmr results')
  parser.add_argument('--counts', required='true', help='overall counts')
  parser.add_argument('--vcfs', required='true', nargs='+', help='vcf files')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--transcribed_strand', action='store_true', help='use mutation on transcribed strand')
  parser.add_argument('--complement', action='store_true', help='normalize reverse complements')
  parser.add_argument('--rotate_context', action='store_true', help='rotate contexts')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.mmr, args.counts, args.vcfs, args.rotate_context, args.transcribed_strand, args.complement)

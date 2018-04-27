#!/usr/bin/env python
'''
  copies across all input files
'''

import logging
import os
import sys

LOGGING_LEVEL=logging.INFO
METADATA='cfg/sample-metadata.csv'
SAMPLES='cfg/samples'
TUMOURS='cfg/tumours'
SOURCE='/scratch/VR0211/pan-prostate/out'

def run(cmd):
  logging.debug('executing "{}"...'.format(cmd))
  result = os.system(cmd)
  if result == 0:
    logging.debug('execution successful')
  else:
    logging.error('non-zero return from "{}"'.format(cmd))

def is_normal(sample, db):
  for line in open(METADATA, 'r'):
    fields = line.strip('\n').split(',')
    if fields[0] == sample:
      return fields[4] == 'Y'

  logging.error('{} not found in metadata'.format(sample))

def make_db():
  samples = set()
  for line in open(SAMPLES, 'r'):
    samples.add(line.strip('\n'))
  normal_to_patient = {}
  patient_to_tumour = {}
  for line in open(METADATA, 'r'):
    fields = line.strip('\n').split(',')
    if fields[4] == 'Y': # normal
      normal_to_patient[fields[0]] = fields[1]
    else: # tumour
      if fields[0] in samples: # only point to sample that is in our sample set
        patient_to_tumour[fields[1]] = fields[0]
 
  return (normal_to_patient, patient_to_tumour)

def find_tumour(sample, db):
  if sample not in db[0]:
    logging.warn('Failed to find patient for sample {}.'.format(sample))
    return None

  patient = db[0][sample]
  if patient not in db[1]:
    logging.warn('Failed to find tumour for sample {}.'.format(sample))
    return None

  return db[1][patient]

def main():
  logging.info('starting...')
  run("src/make_uuids.sh > cfg/samples.uuids")
  db = make_db()
  samples = set()
  tumours = set()
  for line in open('cfg/samples', 'r'):
    sample = line.strip('\n')
    # varscan
    if os.path.isfile('in/{sample}.varscan.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.varscan.vcf'.format(sample=sample))
    else:
      run('ln -s {source}/{sample}.varscan.high.vcf in/{sample}.varscan.vcf'.format(source=SOURCE, sample=sample)) # snvs
  
    # varscan indels
    if os.path.isfile('in/{sample}.varscan_indel.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.varscan_indel.vcf'.format(sample=sample))
    else:
      run('ln -s {source}/{sample}.varscan_indel.high.vcf in/{sample}.varscan_indel.vcf'.format(source=SOURCE, sample=sample)) # indels
  
    # gridss
    if os.path.isfile('in/{sample}.gridss.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.gridss.vcf'.format(sample=sample))
    else:
      if is_normal(sample, db):
        # gridss - normals 
        tumour = find_tumour(sample, db)
        if tumour is None:
          logging.warn('Skipping {}: unable to find tumour'.format(sample))
          continue
        run('src/make_gridss.py --sample {sample}.mapped.bam < {source}/{tumour}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=SOURCE, sample=sample, tumour=tumour))
      else:
        # gridss - tumours
        run('src/make_gridss.py --sample {sample}.mapped.bam < {source}/{sample}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=SOURCE, sample=sample))
  
    # hmmcopy
    if os.path.isfile('in/{sample}.hmmcopy'.format(sample=sample)):
      logging.debug('skipping {sample}.hmmcopy'.format(sample=sample))
    else:
      if os.path.isfile('{source}/{sample}.hmmcopy/normal_segments.txt'.format(source=SOURCE, sample=sample)):
        run('ln -s {source}/{sample}.hmmcopy/normal_segments.txt in/{sample}.hmmcopy'.format(source=SOURCE, sample=sample))
      elif os.path.isfile('{source}/{sample}.hmmcopy/somatic_segments.txt'.format(source=SOURCE, sample=sample)):
        run('ln -s {source}/{sample}.hmmcopy/somatic_segments.txt in/{sample}.hmmcopy'.format(source=SOURCE, sample=sample))

    if not is_normal(sample, db):
      tumours.add(sample)
    
    samples.add(sample)

  logging.debug('making config.yaml file...')
  samples_list = ',\n  '.join(["'{}'".format(sample) for sample in sorted(samples)])
  tumours_list = ',\n  '.join(["'{}'".format(sample) for sample in sorted(tumours)])
  with open('cfg/config.yaml', 'w') as fh:
    for line in open('cfg/config.yaml.template', 'r'):
      line = line.replace('SAMPLES', samples_list)
      line = line.replace('TUMOURS', tumours_list)
      fh.write(line)
  logging.info('done')

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=LOGGING_LEVEL)
  main()

#!/usr/bin/env python
'''
  copies across all input files
'''

import logging
import os
import sys

LOGGING_LEVEL=logging.DEBUG

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
  tumour_to_patient = {}
  patient_to_tumour = {}
  patient_to_normal = {}
  for line in open(METADATA, 'r'):
    fields = line.strip('\n').split(',')
    if fields[4] == 'Y': # normal
      normal_to_patient[fields[0]] = fields[1]
      patient_to_normal[fields[1]] = fields[0]
    else: # tumour
      tumour_to_patient[fields[0]] = fields[1]
      if fields[0] in samples: # only point to sample that is in our sample set
        patient_to_tumour[fields[1]] = fields[0]
 
  return (normal_to_patient, patient_to_tumour, tumour_to_patient, patient_to_normal)

def find_normal(sample, db):
  if sample not in db[2]:
    logging.warn('Failed to find patient for sample {}.'.format(sample))
    return None

  patient = db[2][sample]
  if patient not in db[3]:
    logging.warn('Failed to find normal for sample {}.'.format(sample))
    return None

  return db[3][patient]

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
  
    # varscan indel
    if os.path.isfile('in/{sample}.varscan_indel.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.varscan_indel.vcf'.format(sample=sample))
    else:
      #run('ln -s {source}/{sample}.varscan_indel.high.vcf in/{sample}.varscan_indel.vcf'.format(source=SOURCE, sample=sample)) # indels
      run('ln -s {source}/{sample}.varscan_indel.vcf in/{sample}.varscan_indel.vcf'.format(source=SOURCE, sample=sample)) # indels

    # pindel
    if os.path.isfile('in/{sample}.pindel.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.pindel.vcf'.format(sample=sample))
    else:
      if is_normal(sample, db):
        tumour = find_tumour(sample, db)
        if tumour is None:
          logging.warn('Skipping {}: unable to find tumour'.format(sample))
          continue
        run('gunzip < {source}/{tumour}.pindel-1.1.2.vcf.gz | src/extract_sample.py --minpindel 0.01 --sample NORMAL --verbose > in/{sample}.pindel.vcf'.format(source=SOURCE, tumour=tumour, sample=sample))
      else:
        run('gunzip < {source}/{sample}.pindel-1.1.2.vcf.gz | src/extract_sample.py --minpindel 0.01 --sample TUMOUR --verbose > in/{sample}.pindel.vcf'.format(source=SOURCE, sample=sample))
  
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
        run('src/extract_sample.py --sample {sample}.mapped.bam --minqual 500 --verbose < {source}/{tumour}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=SOURCE, sample=sample, tumour=tumour))
      else:
        # gridss - tumours
        run('src/extract_sample.py --sample {sample}.mapped.bam --minqual 500 --verbose < {source}/{sample}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=SOURCE, sample=sample))
  
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

  # exclude any tumour that doesn't have a normal
  for tumour in list(tumours): # list lets us modify it
    normal = find_normal(tumour, db)
    if normal is None:
      tumours.remove(tumour)
      continue
    if normal not in samples:
      logging.info('Excluding tumour %s because normal %s is unavailable', tumour, normal)
      tumours.remove(tumour)
      continue

  logging.info('making config.yaml file. %i tumours, %i samples.', len(tumours), len(samples))
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

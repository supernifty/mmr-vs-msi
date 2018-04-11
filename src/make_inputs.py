#!/usr/bin/env python

import logging
import os
import sys

LOGGING_LEVEL=logging.INFO
METADATA='cfg/sample-metadata.csv'
SAMPLES='cfg/samples'

def run(cmd):
  logging.info('executing "{}"...'.format(cmd))
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
  patient = db[0][sample]
  return db[1][patient]

def main():
  logging.info('starting...')
  run("src/make_uuids.sh > cfg/samples.uuids")
  db = make_db()
  samples = set()
  for line in open('cfg/samples', 'r'):
    sample = line.strip('\n')
    # varscan
    if os.path.isfile('in/{sample}.varscan.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.varscan.vcf'.format(sample=sample))
    else:
      run('ln -s /scratch/VR0211/pan-prostate/out/{sample}.varscan.vcf in/'.format(sample=sample)) # snvs
  
    # varscan indels
    if os.path.isfile('in/{sample}.varscan_indel.vcf'.format(sample=sample)):
      logging.debug('skipping {sample}.varscan_indel.vcf'.format(sample=sample))
    else:
      run('ln -s /scratch/VR0211/pan-prostate/out/{sample}.varscan_indel.vcf in/'.format(sample=sample)) # indels
  
    if is_normal(sample, db):
      # gridss - normals 
      tumour = find_tumour(sample, db)
      run('src/make_gridss.py --sample {sample}.mapped.bam < /scratch/VR0211/pan-prostate/out/{tumour}.gridss.high.vcf > in/{sample}.gridss.vcf'.format(sample=sample, tumour=tumour))
    else:
      # gridss - tumours
      run('src/make_gridss.py --sample {sample}.mapped.bam < /scratch/VR0211/pan-prostate/out/{sample}.gridss.high.vcf > in/{sample}.gridss.vcf'.format(sample=sample))
  
    samples.add(sample)
  
  logging.debug('making config.yaml file...')
  samples_list = ',\n  '.join(["'{}'".format(sample) for sample in sorted(samples)])
  with open('cfg/config.yaml', 'w') as fh:
    for line in open('cfg/config.yaml.template', 'r'):
      line = line.replace('SAMPLES', samples_list)
      fh.write(line)
  logging.info('done')

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=LOGGING_LEVEL)
  main()

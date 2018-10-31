#!/usr/bin/env python
'''
  copies across all input files
'''

import argparse
import logging
import os
import sys

LOGGING_LEVEL=logging.DEBUG

METADATA='cfg/sample-metadata.csv'
SAMPLES='cfg/samples'
TUMOURS='cfg/tumours'

#INCLUDE=('make_uuids', 'varscan', 'varscan_indel', 'gridss', 'hmmcopy', 'pindel', 'platypus', 'platypus_crc', 'strelka', 'strelka_indel', 'caveman', 'config')
#INCLUDE=('strelka', 'strelka_indel', 'config')
INCLUDE=('config',)

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

def main(source):
  logging.info('starting...')
  if 'make_uuids' in INCLUDE:
    run("src/make_uuids.sh > cfg/samples.uuids")
  db = make_db()
  samples = set()
  tumours = set()
  for line in open('cfg/samples', 'r'):
    sample = line.strip('\n')
    # varscan
    if 'varscan' in INCLUDE:
      if os.path.isfile('in/{sample}.varscan.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.varscan.vcf'.format(sample=sample))
      else:
        run('ln -s {source}/{sample}.varscan.high.vcf in/{sample}.varscan.vcf'.format(source=source, sample=sample)) # snvs
  
    # varscan indel
    if 'varscan_indel' in INCLUDE:
      if os.path.isfile('in/{sample}.varscan_indel.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.varscan_indel.vcf'.format(sample=sample))
      else:
        #run('ln -s {source}/{sample}.varscan_indel.high.vcf in/{sample}.varscan_indel.vcf'.format(source=source, sample=sample)) # indels
        run('ln -s {source}/{sample}.varscan_indel.vcf in/{sample}.varscan_indel.vcf'.format(source=source, sample=sample)) # indels

    # pindel
    if 'pindel' in INCLUDE:
      if os.path.isfile('in/{sample}.pindel.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.pindel.vcf'.format(sample=sample))
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          run('gunzip < {source}/{tumour}.pindel-1.1.2.vcf.gz | src/extract_sample.py --minpindel 0.01 --sample NORMAL --verbose > in/{sample}.pindel.vcf'.format(source=source, tumour=tumour, sample=sample))
        else:
          run('gunzip < {source}/{sample}.pindel-1.1.2.vcf.gz | src/extract_sample.py --minpindel 0.01 --sample TUMOUR --verbose > in/{sample}.pindel.vcf'.format(source=source, sample=sample))
  
    # gridss
    if 'gridss' in INCLUDE:
      if os.path.isfile('in/{sample}.gridss.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.gridss.vcf'.format(sample=sample))
      else:
        if is_normal(sample, db):
          # gridss - normals 
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          run('src/extract_sample.py --sample {sample}.mapped.bam --minqual 500 --verbose < {source}/{tumour}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=source, sample=sample, tumour=tumour))
        else:
          # gridss - tumours
          run('src/extract_sample.py --sample {sample}.mapped.bam --minqual 500 --verbose < {source}/{sample}.gridss.vcf > in/{sample}.gridss.vcf'.format(source=source, sample=sample))
  
    # hmmcopy
    if 'hmmcopy' in INCLUDE:
      if os.path.isfile('in/{sample}.hmmcopy'.format(sample=sample)):
        logging.debug('skipping {sample}.hmmcopy'.format(sample=sample))
      else:
        if os.path.isfile('{source}/{sample}.hmmcopy/normal_segments.txt'.format(source=source, sample=sample)):
          run('ln -s {source}/{sample}.hmmcopy/normal_segments.txt in/{sample}.hmmcopy'.format(source=source, sample=sample))
        elif os.path.isfile('{source}/{sample}.hmmcopy/somatic_segments.txt'.format(source=source, sample=sample)):
          run('ln -s {source}/{sample}.hmmcopy/somatic_segments.txt in/{sample}.hmmcopy'.format(source=source, sample=sample))

    if not is_normal(sample, db):
      tumours.add(sample)
    
    samples.add(sample)

    # platypus /scratch/VR0211/pan-prostate/out/CMHS94.platypus.somatic.vcf
    if 'platypus' in INCLUDE:
      if False and os.path.isfile('in/{sample}.platypus.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.platypus.vcf'.format(sample=sample)) # TODO this is wrong
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          run('src/extract_sample.py --sample {sample}.mapped.bam --filter_homref --verbose < {source}/{tumour}.platypus.vcf > in/{sample}.platypus.vcf'.format(source=source, sample=sample, tumour=tumour))
        else: # tumour
          run('src/extract_sample.py --sample {sample}.mapped.bam --filter_homref --verbose < {source}/{sample}.platypus.vcf > in/{sample}.platypus.vcf'.format(source=source, sample=sample))

    if 'platypus_crc' in INCLUDE:
      if False and os.path.isfile('in/{sample}.platypus.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.platypus.vcf'.format(sample=sample)) # TODO this is wrong
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          run('src/extract_sample.py --sample {sample}.sorted.dups --filter_homref --verbose < {source}/{tumour}.platypus.joint.vcf.gz > in/{sample}.platypus.vcf'.format(source=source, sample=sample, tumour=tumour))
        else: # tumour
          run('src/extract_sample.py --sample {sample}.sorted.dups --filter_homref --verbose < {source}/{sample}.platypus.somatic.vcf.gz > in/{sample}.platypus.vcf'.format(source=source, sample=sample))

    if 'strelka' in INCLUDE:
      if os.path.isfile('in/{sample}.strelka.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.strelka.vcf'.format(sample=sample)) 
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          #run('src/extract_sample.py --sample NORMAL --verbose < {source}/{tumour}.strelka.somatic.snvs.af.norm.vcf.gz > in/{sample}.strelka.vcf'.format(source=source, sample=sample, tumour=tumour))
          # use the actual germline calls
          run('src/extract_sample.py --sample {sample} --filter_homref --verbose < {source}/{sample}.strelka.germline.filter_gt.vep.vcf.gz > in/{sample}.strelka.vcf'.format(source=source, sample=sample, tumour=tumour))
        else: # copy the tumour
          run('src/extract_sample.py --sample TUMOR --verbose < {source}/{sample}.strelka.somatic.snvs.af.norm.vcf.gz > in/{sample}.strelka.vcf'.format(source=source, sample=sample))

    if 'strelka_indel' in INCLUDE:
      if os.path.isfile('in/{sample}.strelka_indel.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.strelka_indel.vcf'.format(sample=sample)) 
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          #run('src/extract_sample.py --sample NORMAL --verbose < {source}/{tumour}.strelka.somatic.indels.vep.vcf.gz > in/{sample}.strelka_indel.vcf'.format(source=source, sample=sample, tumour=tumour))
          # use the actual germline calls
          run('src/extract_sample.py --sample {sample} --filter_homref --verbose < {source}/{sample}.strelka.germline.filter_gt.vep.vcf.gz > in/{sample}.strelka_indel.vcf'.format(source=source, sample=sample, tumour=tumour))
        else: # tumour
          run('src/extract_sample.py --sample TUMOR --verbose < {source}/{sample}.strelka.somatic.indels.vep.vcf.gz > in/{sample}.strelka_indel.vcf'.format(source=source, sample=sample))

    # caveman CMHS1.caveman-1.1.2.flagged.muts.vcf
    # caveman marks all germline 0|0 so it will contain no variants
    if 'caveman' in INCLUDE:
      if os.path.isfile('in/{sample}.caveman.vcf'.format(sample=sample)):
        logging.debug('skipping {sample}.caveman.vcf'.format(sample=sample)) 
      else:
        if is_normal(sample, db):
          tumour = find_tumour(sample, db)
          if tumour is None:
            logging.warn('Skipping {}: unable to find tumour'.format(sample))
            continue
          run('src/extract_sample.py --sample NORMAL --filter_homref --verbose < {source}/{tumour}.caveman-1.1.2.flagged.muts.vcf > in/{sample}.caveman.vcf'.format(source=source, sample=sample, tumour=tumour))
        else: # tumour
          run('src/extract_sample.py --sample TUMOUR --filter_homref --verbose < {source}/{sample}.caveman-1.1.2.flagged.muts.vcf > in/{sample}.caveman.vcf'.format(source=source, sample=sample))


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
  if 'config' in INCLUDE:
    with open('cfg/config.yaml', 'w') as fh:
      for line in open('cfg/config.yaml.template', 'r'):
        line = line.replace('SAMPLES', samples_list)
        line = line.replace('TUMOURS', tumours_list)
        fh.write(line)
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Copy input files')
  parser.add_argument('--source', required=True, help='root directory for input vcfs')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.source)

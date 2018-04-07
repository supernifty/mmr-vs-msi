#!/usr/bin/env python

import os
import sys

def run(cmd):
  sys.stderr.write('  executing "{}"...'.format(cmd))
  os.system(cmd)
  sys.stderr.write(' done\n')

sys.stderr.write('generating symlinks...')
samples = set()
for line in open('cfg/samples', 'r'):
  sample = line.strip('\n')
  run('ln -s /scratch/VR0211/pan-prostate/out/{sample}.varscan.vcf in/'.format(sample=sample))
  samples.add(sample)

sys.stderr.write('making config.yaml file...')
samples_list = ',\n  '.join(["'{}'".format(sample) for sample in sorted(samples)])
with open('cfg/config.yaml', 'w') as fh:
  for line in open('cfg/config.yaml.template', 'r'):
    line = line.replace('SAMPLES', samples_list)
    fh.write(line)
sys.stderr.write(' done\n')

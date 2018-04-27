#!/usr/bin/env python

import collections
import logging
import os
import sys

def main(vcfs):
  logging.info('processing %i vcfs...', len(vcfs))

  samples = {}
  callers = set()
  
  for idx, vcf in enumerate(vcfs):
    sample = os.path.basename(vcf).split('.')[0]
    caller = os.path.basename(vcf).split('.')[1]
    callers.add(caller)
    if sample not in samples:
      samples[sample] = collections.defaultdict(int)
    for line in open(vcf, 'r'):
      if not line.startswith('#'):
        samples[sample][caller] += 1
    logging.debug('processed %s: %i of %i vcfs', vcf, idx + 1, len(vcfs))

  logging.info('processing %i vcfs: done', len(vcfs))
    
  # write results
  sys.stdout.write("Name\t{}\n".format('\t'.join(sorted(list(callers)))))
  for sample in sorted(samples):
    sys.stdout.write('{sample}\t{callers}\n'.format(sample=sample, callers='\t'.join([str(samples[sample][caller]) for caller in sorted(list(callers))])))

if __name__ == '__main__':
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  main(sys.argv[1:])

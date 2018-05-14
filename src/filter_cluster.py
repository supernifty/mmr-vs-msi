#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys


def main(fh, out):
  logging.info('starting...')

  #"Sample","Patient","Lab ID","Tissue","Sequenced","Status","Tags","Notes"
  #"CMHS319","CMHP125","PP74_RARP","Primary","Library","ContEst","","Inferred barcode"

  primaries = set()
  for line in open('cfg/htsdb.csv', 'r'):
    fields = line.strip('\n').split(',')
    if len(fields) < 3:
      continue
    if fields[3].strip('"') == 'Primary':
      primaries.add(fields[0].strip('"'))

  logging.info('%i primaries', len(primaries))

  first = True
  columns = set()
  # Chrom   Begin   End     Annotation      Len     Samples CMHS1   CMHS101 CMHS104 CMHS105
  for line in fh:
    fields = line.strip('\n').split('\t')
    if first:
      first = False
      logging.info('originally %i columns', len(fields))
      for idx in range(len(fields)):
        if idx <= 5 or fields[idx] in primaries:
          columns.add(idx)
      logging.info('post %i columns', len(columns))
    # write filtered row
    new_row = []
    for idx in range(len(fields)):
      if idx in columns:
        new_row.append(fields[idx])
    out.write('{}\n'.format('\t'.join(new_row)))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(sys.stdin, sys.stdout)

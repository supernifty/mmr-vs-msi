#!/usr/bin/env python

# filters an msi bed file

import argparse
import sys

def get_fields(values):
  result = {}
  for value in values:
    key, val = value.split('=')
    result[key] = val
  return result

def main(src, target, minlen, maxlen, minrepeat, maxrepeat, exons_only):
  allowed = 0
  denied = 0
  for idx, line in enumerate(src):
    params = get_fields(line.strip('\n').split('\t')[3].split(';'))
    if minlen <= int(params['length']) <= maxlen:
      # test repeat length
      if minrepeat <= len(params['repeat']) <= maxrepeat:
        if not exons_only or 'exon' in params:
          allowed += 1
          if line.startswith('chr'):
            line = line[3:]
          target.write(line)
        else:
          denied += 1
      else:
        denied += 1
    else:
      denied += 1

    if idx % 1000000 == 0:
      sys.stderr.write('written {} of {}.\n'.format(allowed, allowed + denied))

  sys.stderr.write('finished. wrote {} of {}.\n'.format(allowed, allowed + denied))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter MSI bed')
  parser.add_argument('--minlen', type=int, default=10)
  parser.add_argument('--maxlen', type=int, default=100)
  parser.add_argument('--minrepeat', type=int, default=1)
  parser.add_argument('--maxrepeat', type=int, default=4)
  parser.add_argument('--exons_only', action='store_true')
  args = parser.parse_args()
  main(sys.stdin, sys.stdout, minlen=args.minlen, maxlen=args.maxlen, minrepeat=args.minrepeat, maxrepeat=args.maxrepeat, exons_only=args.exons_only)

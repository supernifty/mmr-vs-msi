#!/usr/bin/env python

import collections
import sys

def main(fh_in, out, err):
  repeat = collections.defaultdict(int)
  length = collections.defaultdict(int)
  repeat_exon = collections.defaultdict(int)
  length_exon = collections.defaultdict(int)
  total = 0
  exon_total = 0
  for total, line in enumerate(fh_in):
    fields = line.strip('\n').split('\t')
    if 'exon' in fields[3]:
      exon_total += 1
    # chr1    10103   10113   repeat=AACCC;length=10
    #detail = {}
    for kv in fields[3].split(';'):
      k, v = kv.split('=')
      #detail[k] = v
      if k == 'repeat':
        repeat[v] += 1
        if 'exon' in fields[3]:
          repeat_exon[v] += 1
      if k == 'length':
        length[int(v)] += 1
        if 'exon' in fields[3]:
          length_exon[int(v)] += 1

    if total % 100000 == 0:
      sys.stderr.write('{} lines...\n'.format(total))
      
  out.write('Type\n')
  for key in sorted(repeat.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat[key], repeat[key] / (total + 1)))

  out.write('\nType in exons\n')
  for key in sorted(repeat_exon.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat_exon[key], repeat_exon[key] / exon_total))

  out.write('\nLength\n')
  for key in sorted(length.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, length[key], length[key] / total))

  out.write('\nLength in exons\n')
  for key in sorted(length_exon.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, length_exon[key], length_exon[key] / exon_total))

if __name__ == '__main__':
  main(sys.stdin, sys.stdout, sys.stderr)

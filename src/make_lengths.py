#!/usr/bin/env python
'''
  given a fasta file, generate a lengths file as required by bedtools

  usage: make_lengths < fasta > lengths
'''

import sys

chrom = None
length = 0
seen = set()
for idx, line in enumerate(sys.stdin):
  if line.startswith('>'):
    if chrom is not None and chrom not in seen:
      sys.stdout.write('{chrom}\t{length}\n'.format(chrom=chrom, length=length))
      seen.add(chrom)
    chrom = line.strip('\n').split(' ')[0][1:]
    length = 0
  else:
    length += len(line.strip('\n'))

  if idx < 1000000 and idx % 100000 == 0 or idx % 1000000 == 0:
    sys.stderr.write('{lines} lines processed. chromosome {chrom}:{length}...\n'.format(lines=idx, chrom=chrom, length=length))

if length > 0:
  sys.stdout.write('{chrom}\t{length}\n'.format(chrom=chrom, length=length))

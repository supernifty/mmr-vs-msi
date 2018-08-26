#!/usr/bin/env python

import collections
import sys

def rotate_repeat(r):
  '''
    normalize repeats to start from earliest base e.g. GCAA -> AAGC
  '''
  best = 0
  for idx, el in enumerate(r):
    if r[idx] < r[best]:
      best = idx

  return r[best:] + r[:best]

def main(fh_in, out, err):
  repeat = collections.defaultdict(int) # counts of repeat types i.e. AGGC
  repeat_rotated = collections.defaultdict(int) # normalized repeat types
  length = collections.defaultdict(int) # repeat lengths
  repeat_exon = collections.defaultdict(int) # limited to exons
  repeat_exon_rotated = collections.defaultdict(int) # limited to exons
  length_exon = collections.defaultdict(int) # limited to exons
  total = total_seq = 0 # total number of entries
  exon_total = 0 # total number of exon repeats
  five_utr = three_utr = five_utr_seq = three_utr_seq = exon_seq = ncrna = ncrna_seq = 0
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
        repeat_rotated[rotate_repeat(v)] += 1
        if 'exon' in fields[3]:
          repeat_exon[v] += 1
          repeat_exon_rotated[rotate_repeat(v)] += 1
      if k == 'length':
        length[int(v)] += 1
        total_seq += int(v)
        if 'exon' in fields[3]:
          length_exon[int(v)] += 1
          exon_seq += int(v)
        if '5UTR' in fields[3]:
          five_utr += 1
          five_utr_seq += int(v)
        if '3UTR' in fields[3]:
          three_utr += 1
          three_utr_seq += int(v)
        if 'ncRNA' in fields[3]:
          ncrna += 1
          ncrna_seq += int(v)

    if total % 100000 == 0:
      sys.stderr.write('{} lines...\n'.format(total))
      
  out.write('## Totals\n')
  out.write('Exon count\t{}\t{:.3f}\n'.format(exon_total, exon_total / (total + 1)))
  out.write('5UTR count\t{}\t{:.3f}\n'.format(five_utr, five_utr / (total + 1)))
  out.write('3UTR count\t{}\t{:.3f}\n'.format(three_utr, three_utr / (total + 1)))
  out.write('ncRNA count\t{}\t{:.3f}\n'.format(ncrna, ncrna / (total + 1)))
  out.write('Exon sequence\t{}\t{:.3f}\n'.format(exon_seq, exon_seq / total_seq))
  out.write('5UTR sequence\t{}\t{:.3f}\n'.format(five_utr_seq, five_utr_seq / total_seq))
  out.write('3UTR sequence\t{}\t{:.3f}\n'.format(three_utr_seq, three_utr_seq / total_seq))
  out.write('ncRNA sequence\t{}\t{:.3f}\n'.format(ncrna_seq, ncrna_seq / total_seq))

  out.write('## Type\n')
  for key in sorted(repeat.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat[key], repeat[key] / (total + 1)))

  out.write('\n## Type rotated\n')
  for key in sorted(repeat_rotated.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat_rotated[key], repeat_rotated[key] / (total + 1)))

  out.write('\n## Type in exons rotated\n')
  for key in sorted(repeat_exon_rotated.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat_exon_rotated[key], repeat_exon_rotated[key] / exon_total))

  out.write('\n## Type in exons\n')
  for key in sorted(repeat_exon.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, repeat_exon[key], repeat_exon[key] / exon_total))

  out.write('\n## Length\n')
  for key in sorted(length.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, length[key], length[key] / total))

  out.write('\n## Length in exons\n')
  for key in sorted(length_exon.keys()):
    out.write('{}\t{}\t{:.3f}\n'.format(key, length_exon[key], length_exon[key] / exon_total))

if __name__ == '__main__':
  main(sys.stdin, sys.stdout, sys.stderr)

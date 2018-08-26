#!/usr/bin/env python

import argparse
import logging
import sys

REMOVE_CHR=False
SKIP_CONTIGS=True

def main(filter, annotate):
  sys.stdin.readline() # header
  chroms = set()
  sequence = {'coding': 0, 'five': 0, 'three': 0, 'ncRNA': 0}
  for line in sys.stdin:
    #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
    #0       NM_001308203.1  chr1    +       66999251        67216822        67000041        67208778        22      66999251,66999928,67091529,67098752,67105459,67108492,67109226,67136677,67137626,67138963,67142686,67145360,67154830,67155872,67160121,67184976,67194946,67199430,67205017,67206340,67206954,67208755,  66999355,67000051,67091593,67098777,67105516,67108547,67109402,67136702,67137678,67139049,67142779,67145435,67154958,67155999,67160187,67185088,67195102,67199563,67205220,67206405,67207119,67216822,  0       SGIP1   cmpl    cmpl    -1,0,1,2,0,0,1,0,1,2,1,1,1,0,1,1,2,2,0,2,1,1,
    #1       NM_001323574.1  chr1    -       48998526        50489626        48999844        50489468        14      48998526,49000561,49005313,49052675,49056504,49100164,49119008,49128823,49332862,49511255,49711441,50162948,50317067,50489434,  48999965,49000588,49005410,49052838,49056657,49100276,49119123,49128913,49332902,49511472,49711536,50163109,50317190,50489626,  0       AGBL4   cmpl    cmpl    2,2,1,0,0,2,1,1,0,2,0,1,1,0,
    fields = line.strip('\n').split('\t')
    codingStart = int(fields[6])
    codingEnd = int(fields[7])
    exonStarts = fields[9].split(',')
    exonEnds = fields[10].split(',')
    for rank, (start, end) in enumerate(zip(exonStarts, exonEnds)):
      if start == '':
        continue

      chrom = fields[2]
      if SKIP_CONTIGS and '_' in chrom:
        continue

      if REMOVE_CHR and chrom.startswith('chr'):
        chrom = chrom[3:]

      if fields[3] == '+':
        strand = 'f'
      else:
        strand = 'r'

      start = int(start)
      end = int(end)

      # update stats
      if codingStart > start and codingStart < end or codingEnd < end and codingEnd > start: # one point in the coding region
        cstart = max(start, codingStart)
        cend = min(end, codingEnd)
        sequence['coding'] += cend - cstart

      if filter == '3UTR':
        if fields[1].startswith('NR_'):
          continue
        if strand == 'f':
          if end <= codingEnd: # nothing to include
            continue
          start = max(start, codingEnd) # potential truncation
        else:
          if start >= codingStart: # nothing to include
            continue
          end = min(end, codingStart) # potential truncation

      if filter == '5UTR':
        if fields[1].startswith('NR_'):
          continue
        if strand == 'f':
          if start >= codingStart: # nothing to include
            continue
          end = min(end, codingStart) # potential truncation
        else:
          if end <= codingEnd: # nothing to include
            continue
          start = max(start, codingEnd) # potential truncation

      if filter == 'ncRNA':
        if fields[1].startswith('NM_'):
          continue
        if strand == 'f':
          if start >= codingStart: # nothing to include
            continue
          end = min(end, codingStart) # potential truncation
        else:
          if end <= codingEnd: # nothing to include
            continue
          start = max(start, codingEnd) # potential truncation

      if annotate:
        annotation = "{gene}_exon_{rank}_{count}_{vcf_pos}_{strand}".format(gene=fields[12], rank=rank, count=fields[8], vcf_pos=int(start) + 1, strand=strand)
        sys.stdout.write('{chrom}\t{start}\t{end}\t{annotation}\n'.format(chrom=chrom, start=start, end=end, annotation=annotation))
      else:
        sys.stdout.write('{chrom}\t{start}\t{end}\n'.format(chrom=chrom, start=start, end=end))
      chroms.add(chrom)

  logging.info('chroms:\n  %s', '\n  '.join(sorted(list(chroms))))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate a bed file from refseq')
  parser.add_argument('--filter', required=False, help='3UTR or 5UTR or ncRNA, omit for all')
  parser.add_argument('--annotate', action='store_true', help='include annotation')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)
  main(args.filter, args.annotate)

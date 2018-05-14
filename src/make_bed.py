#!/usr/bin/env python
'''
  given ucsc output, generate transcript locations for specified genes
 
  Usage:
    make_bed genes.txt < ucsc.out > bed
'''

REMOVE_CHR=True

import sys

genes = set()
for gene_list in sys.argv[1:]:
  for line in open(gene_list):
    genes.add(line.strip('\n'))

found = set()

# new
#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
#0       NM_001308203.1  chr1    +       66999251        67216822        67000041        67208778        22      66999251,66999928,67091529,67098752,67105459,67108492,67109226,67136677,67137626,67138963,67142686,67145360,67154830,67155872,67160121,67184976,67194946,67199430,67205017,67206340,67206954,67208755,      66999355,67000051,67091593,67098777,67105516,67108547,67109402,67136702,67137678,67139049,67142779,67145435,67154958,67155999,67160187,67185088,67195102,67199563,67205220,67206405,67207119,67216822,      0       SGIP1   cmpl    cmpl    -1,0,1,2,0,0,1,0,1,2,1,1,1,0,1,1,2,2,0,2,1,1,

written = 0
read = 0
for read, line in enumerate(sys.stdin):
  fields = line.strip('\n').split('\t')
  gene = fields[12]
  if gene in genes:
    found.add(gene)
    chrom = fields[2]
    if REMOVE_CHR and chrom.startswith('chr'):
      chrom = chrom[3:]
    sys.stdout.write('{chrom}\t{start}\t{end}\t{gene}\t{transcript}\n'.format(chrom=chrom, start=fields[4], end=fields[5], gene=gene, transcript=fields[1]))
    written += 1

sys.stderr.write('Read {}, wrote {}\n'.format(read, written))
sys.stderr.write('Genes not found: {}\n'.format(', '.join(sorted(list(genes.difference(found))))))

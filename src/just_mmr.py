#!/usr/bin/env python

import sys

def main(summary_fh, genes_fn, out):
  genes = set()
  for line in open(genes_fn, 'r'):
    genes.add(line.strip('\n'))

  first = True
  cols = set()
  for line in summary_fh:
    fields = line.strip('\n').split('\t')
    if first:
      first = False
      for col in range(len(fields)):
        if '|' in fields[col]:
          gene = fields[col].split('|')[0]
          if gene in genes:
            cols.add(col)
        elif fields[col].startswith('rt_'):
          pass
        else:
          cols.add(col)

    # print allowed cols
    outline = '\t'.join([fields[x] for x in range(len(fields)) if x in cols])
    sys.stdout.write('{}\n'.format(outline))

if __name__ == '__main__':
  main(sys.stdin, sys.argv[1], sys.stdout)

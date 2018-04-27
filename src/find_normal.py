#!/usr/bin/env python

import sys

def main(metadata, sample):
  patient = None
  normals = {}
  # CMHS1,CMHP1,299bT,BT,N,,,
  for line in metadata:
    fields = line.strip('\n').split(',')
    if fields[0] == sample:
      patient = fields[1]
    if fields[4] == 'Y':
      normals[fields[1]] = fields[0]

  if patient in normals:
    sys.stdout.write('{}\n'.format(normals[patient]))
  else:
    sys.stdout.write('NOTFOUND\n')

if __name__ == '__main__':
  main(sys.stdin, sys.argv[1])

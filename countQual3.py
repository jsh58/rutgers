#!/usr/bin/python

# JMG 6/24/16
# Retrieving a subset of reads from a SAM.
# Counting good/bad qual methylation data
# Version 2: querying a fastq file directly
# Version 3: using numpy to calculate avg/stddev

import sys
import gzip
import numpy

def parseFASTQ(f, lis):
  '''
  Parse the FASTQ file. Record quality scores for each base.
  '''
  total = 0
  qual = [[] for i in range(5)]
  head = f.readline()
  while head:
    seq = f.readline().rstrip()
    plus = f.readline()
    qul = f.readline().rstrip()

    # record qual scores
    for i in range(len(seq)):
      idx = lis.index(seq[i])
      qual[idx].append( ord(qul[i]) - 33 )

    total += 1
    head = f.readline()

  # calculate average quality
  avgQual = []
  for i in range(len(qual)):
    avgQual.append(( numpy.mean( qual[i] ), numpy.std( qual[i] ), len( qual[i] ) ))
  return avgQual, total


def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 1:
    print 'Usage: python %s  <FASTQfile>' % sys.argv[0]
    print '  file must be gzip-compressed'
    sys.exit(-1)

  try:
    f = gzip.open(args[0], 'rb')
  except IOError:
    sys.stderr.write('Error! Cannot open %s\n' % args[0])
    sys.exit(-1)

  lis = ['A', 'C', 'G', 'T', 'N']
  avgQual, total = parseFASTQ(f, lis)
  f.close()

  for i in range(len(lis)):
    print lis[i],
    for j in range(len(avgQual[i])):
      print '\t', avgQual[i][j],
    print

  sys.stderr.write('Reads analyzed: %d\n' % total)

if __name__ == '__main__':
  main()

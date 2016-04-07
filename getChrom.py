#!/usr/bin/python

# JMG 3/17/16
# Retrieve a segment from a genome.

import sys

def main():
  args = sys.argv[1:]
  if len(args) < 4:
    print 'Usage: python getChrom.py  <fastaGenome>  <chrom>  <pos>  <len>'
    print '  note: <pos> is 0-based position'
    sys.exit(-1)

  try:
    f = open(args[0], 'rU')
  except IOError:
    print 'Error! Cannot open genome file %s' % args[0]
    sys.exit(-1)

  chrom = '>' + args[1]
  try:
    pos = int(args[2])
    length = int(args[3])
  except ValueError:
    print 'Error! <pos> and <len> must be integers'
    sys.exit(-1)

  chr = ''
  for line in f:
    if line.rstrip() == chrom:
      # save chromosome
      for line in f:
        if line[0] == '>':
          break
        chr += line.rstrip()
      break
  f.close()

  # print output
  if chr:
    print '%s %d %d\n%s' % (chrom, pos, length, chr[pos:pos+length])
  else:
    print 'Error! Chromosome %s not found' % chrom[1:]

if __name__ == '__main__':
  main()

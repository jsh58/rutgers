#!/usr/bin/python

# JMG 3/17/16
# Retrieve a segment from a genome.

import sys

def revComp(seq):
  rc = ''
  comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  for nuc in seq[::-1]:
    if not nuc.upper() in comp:
      print 'Error! Unknown nucleotide:', nuc
      sys.exit(-1)
    rc += comp[nuc.upper()]
  return rc

def main():
  args = sys.argv[1:]
  if len(args) < 4:
    print 'Usage: python %s  <genome> ' % sys.argv[0] + \
      ' <chrom>  <pos>  <len>  [rc]'
    print '  <genome>  FASTA format'
    print '  <pos>     0-based position'
    print '  [rc]      option to print reverse-complement'
    sys.exit(-1)

  try:
    f = open(args[0], 'rU')
  except IOError:
    print 'Error! Cannot open genome file %s' % args[0]
    sys.exit(-1)

  # save CL args
  chrom = '>' + args[1]
  try:
    pos = int(args[2])
    length = int(args[3])
  except ValueError:
    print 'Error! <pos> and <len> must be integers'
    sys.exit(-1)
  rc = 0
  if len(args) > 4 and args[4] == 'rc':
    rc = 1

  # find chromosome
  chr = ''
  for line in f:
    if line.rstrip() == chrom:
      # save sequence
      for line in f:
        if line[0] == '>':
          break
        chr += line.rstrip()
      break
  f.close()

  # print output
  if chr:
    if pos + length >= len(chr):
      print 'Warning! Chromosome %s length is %d' % \
        (chrom[1:], len(chr))
    print '%s %d %d' % (chrom, pos, length),
    if rc:
      print 'rc\n%s' % revComp(chr[pos:pos+length])
    else:
      print '\n%s' % chr[pos:pos+length].upper()
  else:
    print 'Error! Chromosome %s not found' % chrom[1:]

if __name__ == '__main__':
  main()

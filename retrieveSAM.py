#!/usr/bin/python

# JMG 6/13/16
# Retrieving a subset of reads from a SAM.

import sys
import re

def parseCigar(cigar):
  '''
  Determine indel offset of alignment from CIGAR.
  '''
  offset = 0
  ins = re.findall(r'(\d+)I', cigar)
  for i in ins:
    offset -= int(i)
  de = re.findall(r'(\d+)D', cigar)
  for d in de:
    offset += int(d)
  return offset

def parseSAM(f, fOut, chrom, start, end):
  '''
  Parse the SAM file.
  '''
  count = total = 0
  for line in f:
    if line[0] == '@':
      fOut.write(line)
      continue
    spl = line.rstrip().split('\t')
    count += 1

    # skip if wrong chromosome
    chr = spl[2]
    if chr != chrom:
      continue

    # write if within bounds
    loc = int(spl[3])
    loc3 = loc + len(spl[9]) + parseCigar(spl[5])
    if (loc >= start and loc <= end) or \
       (loc3 >= start and loc3 <= end):
      fOut.write(line)
      total += 1

  return count, total

def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 5:
    print 'Usage: python %s  <SAMfile>  <out>  <chr>  <start>  <end>' % sys.argv[0]
    print '  Use \'-\' for stdin/stdout'
    sys.exit(-1)
  try:
    f = open(args[0], 'rU')
  except IOError:
    # use '-' for stdin
    if args[0] == '-':
      f = sys.stdin
    else:
      sys.stderr.write('Error! Cannot open %s\n' % args[0])
      sys.exit(-1)

  # use '-' for stdout
  if args[1] == '-':
    fOut = sys.stdout
  else:
    fOut = open(args[1], 'w')

  chr = args[2]
  try:
    start = int(args[3])
    end = int(args[4])
  except ValueError:
    sys.stderr.write('Error! Cannot convert to int')
    sys.exit(-1)

  count, total = parseSAM(f, fOut, chr, start, end)
  f.close()
  fOut.close()

  sys.stderr.write('Reads analyzed: %d\n' % count)
  sys.stderr.write('Reads written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()

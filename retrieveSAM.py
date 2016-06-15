#!/usr/bin/python

# JMG 6/13/16
# Retrieving a subset of reads from a SAM.

import sys
import re

def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to int\n' % arg)
    #usage()
  return val

def parseCigar(cigar):
  '''
  Return in/del offset, plus list of tuples for cigar.
  '''
  parts = re.findall(r'(\d+)([IDM])', cigar)
  offset = 0
  for part in parts:
    if part[1] == 'D':
      offset += int(part[0])
    elif part[1] == 'I':
      offset -= int(part[0])
  return offset, parts

def getXM(lis):
  '''
  Get "XM" methylation string.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == 'XM':
      return spl[-1]
  sys.stderr.write('Error! Cannot find XM in SAM record\n')
  sys.exit(-1)

def adjustCig(seq, cigar):
  '''
  Adjust the given string based on the cigar,
    inserting "-" for each D.
  '''
  ans = ''
  pos = 0
  for c in cigar:
    dist = int(c[0])
    if c[1] == 'M':
      ans += seq[pos:pos+dist]
      pos += dist
    elif c[1] == 'I':
      ans += seq[pos:pos+dist]
      pos += dist
    elif c[1] == 'D':
      ans += '-' * dist
  return ans

def loadGen(fGen, chr):
  '''
  Load a chromosome from the given genome file.
  '''
  try:
    f = open(fGen, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open genome file %s\n' % fGen)
    sys.exit(-1)

  # find chromosome
  seq = ''
  for line in f:
    if line.rstrip()[1:] == chr:
      # save sequence
      for line in f:
        if line[0] == '>':
          return seq
        seq += line.rstrip()
  if seq:
    return seq
  sys.stderr.write('Error! Cannot find chromosome %s\n' % chr)
  sys.exit(-1)

def parseSAM(f, fOut, chrom, start, end, meth):
  '''
  Parse the SAM file. Select reads overlapping the
    given coordinates.
  '''
  minLoc = 1000000000
  maxLoc = 0
  count = total = 0
  for line in f:
    if line[0] == '@':
      #fOut.write(line)
      continue
    spl = line.rstrip().split('\t')
    count += 1

    # skip if wrong chromosome
    chr = spl[2]
    if chr != chrom:
      continue

    # write if within bounds
    loc = int(spl[3])
    offset, cigar = parseCigar(spl[5])
    loc3 = loc + len(spl[9]) + offset
    if (loc >= start and loc <= end) or \
        (loc3 >= start and loc3 <= end) or \
        (loc <= start and loc3 >= end):
      meth[spl[0]] = (loc, adjustCig(getXM(spl[11:]), cigar))
      #fOut.write(line)
      total += 1
      if loc3 > maxLoc:
        maxLoc = loc3
      if loc < minLoc:
        minLoc = loc

  return count, total, minLoc, maxLoc

def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 5:
    print 'Usage: python %s  <SAMfile>  <out> ' % sys.argv[0],
    print '<chr>  <start>  <end>  [<genome>]'
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

  # save region of interest
  chr = args[2]
  start = getInt(args[3])
  end = getInt(args[4])

  # process file
  meth = {}
  count, total, minLoc, maxLoc = parseSAM(f, fOut, chr, start, end, meth)
  f.close()

  # get genomic segment
  if len(args) > 5:
    seq = loadGen(args[5], chr)
    fOut.write('>%s %d %d\n' % (chr, minLoc, maxLoc))
    fOut.write(seq[minLoc-1:maxLoc-1].upper() + '\n')
  def firstVal(tup):
    return tup[1][0]
  for read in sorted(meth.items(), key=firstVal):
    #fOut.write(' ' * (meth[read[0]][0] - minLoc) + meth[read[0]][1] + '\n')
    fOut.write(' ' * (meth[read[0]][0] - minLoc) + meth[read[0]][1] + \
      '  ' + read[0] + '\n')
  fOut.close()

  sys.stderr.write('Reads analyzed: %d\n' % count)
  sys.stderr.write('Reads written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()

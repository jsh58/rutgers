#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: determining wtf is wrong with Picard MarkDuplicates

import sys
import gzip
import re

class Read():
  '''
  A class for representing the position(s) of
    a read's alignment(s).
  '''
  def __init__(self, spl):
    self.r1 = []    # alignments of first read
    self.r2 = []    # alignments of second read
    self.pair = []  # paired alignments
    self.header = spl[0]
    self.addInfo(spl)  # save initial alignment

  def addInfo(self, spl):
    chrom, pos, rc, first, hi, ascore, mateChr, matePos \
      = self.loadInfo(spl)
    if first:
      #self.r1.append((chrom, pos))
      self.r1.append('%s:%d:%s:%s' % (chrom, pos, hi, ascore))
      #self.r1.append('%s:%d:%s:%d' % (chrom, pos, mateChr, matePos))
    else:
      #self.r2.append((chrom, pos))
      self.r2.append('%s:%d:%s:%s' % (chrom, pos, hi, ascore))
      #self.r2.append('%s:%d:%s:%d' % (chrom, pos, mateChr, matePos))

  def loadInfo(self, spl):
    '''
    Retrieve the RNAME, POS, and flag information from
      a SAM record.
    '''
    flag = int(spl[1])
    chrom = spl[2]
    pos = int(spl[3])
    rc = 0
    if flag & 0x10:
      # reverse-complement alignment: position is 3' end
      pos += self.parseCigar(spl[5])
      rc = 1
    first = 0   # first segment in template
    if flag & 0x80:
      first = 1  # last segment in template
    mateChr = spl[6]
    matePos = int(spl[7])
    hi = self.getTag(spl[11:], 'HI')
    ascore = self.getTag(spl[11:], 'AS')
    return chrom, pos, rc, first, hi, ascore, mateChr, matePos

  def parseCigar(self, cigar):
    '''
    Determine 3' position of alignment from CIGAR.
    '''
    ops = re.findall(r'(\d+)(\D)', cigar)
    offset = 0
    for op in ops:
      if op[1] in ['M', 'D', 'N', 'S', '=', 'X']:
        offset += int(op[0])
      elif op[1] in ['I', 'P', 'H']:
        pass
      else:
        sys.stderr.write('Error! Unknown op %s in cigar %s\n' \
          % (op[1], cigar))
        sys.exit(-1)
    return offset

  def getTag(self, lis, tag):
    '''
    Get optional tag from a SAM record.
    '''
    for t in lis:
      spl = t.split(':')
      if spl[0] == tag:
        return spl[-1]
    return -1
    #sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
    #sys.exit(-1)

  def getHeader(self):
    return self.header

  def getInfo(self):
    return self.r1, self.r2

  def __str__(self):
    res = self.header + '\n'
    res += 'first: ' + ' '.join(self.r1) + '\n'
    res += 'second: ' + ' '.join(self.r2)
    print res

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def main():
  args = sys.argv[1:]
  if len(args) < 1:
    print 'Usage: python %s  <SAMfile>  [<out>]' % sys.argv[0]
    print '  Use \'-\' for stdin/stdout'
    sys.exit(-1)
  f = openRead(args[0])
  fOut = None
  out = 0
  if len(args) > 1:
    fOut = openWrite(args[1])
    out = 1

  reads = []    # list of reads
  headers = {}  # list of read headers
  count = dups = uniq = notSeq = 0
  for line in f:
    if line[0] == '@':
      #if out: fOut.write(line)
      continue
    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n'
        + line)
      sys.exit(-1)

    if spl[0] in headers:
      for read in reads:
        if read.getHeader() == spl[0]:
          read.addInfo(spl)
          #print read.getInfo()
          #sys.exit(0)
          break
    else:
      reads.append(Read(spl))
      headers[spl[0]] = 1

    #for read in reads:
    #  print read.getInfo()
    #raw_input()
    count += 1
    if count % 100000 == 0:
      print count
    if count % 1000000 == 0:
      break

  # check alignment scores
  for read in reads:
    r1, r2 = read.getInfo()
    best = 1000
    for aln in r1:
      spl = aln.split(':')
      if best == 1000:
        best = spl[-1]
      elif spl[-1] != best:
        print read.getHeader(), r1

    best = 1000
    for aln in r2:
      spl = aln.split(':')
      if best == 1000:
        best = spl[-1]
      elif spl[-1] != best:
        print read.getHeader(), r2

  sys.exit(0)

  f.close()
  if out: fOut.close()

  sys.stderr.write('Reads: %10d\n' % count)
  sys.stderr.write('Unique: %9d\n' % uniq)
  sys.stderr.write('Dups: %11d\n' % dups)
  sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

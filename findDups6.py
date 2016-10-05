#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: determining wtf is wrong with Picard MarkDuplicates

import sys
import re

class Read():
  '''
  A class for representing the position(s) of
    a read's alignment(s).
  '''
  def __init__(self, spl):
    self.r1 = []  # alignments of first read
    self.r2 = []  # alignments of second read
    self.header = spl[0]
    self.addInfo(spl)  # save initial alignment

  def addInfo(self, spl):
    chr, pos, rc, first = self.loadInfo(spl)
    if first:
      #self.r1.append((chr, pos))
      self.r1.append('%s-%d' % (chr, pos))
    else:
      #self.r2.append((chr, pos))
      self.r2.append('%s-%d' % (chr, pos))

  def loadInfo(self, spl):
    '''
    Retrieve the RNAME, POS, and flag information from
      a SAM record.
    '''
    flag = int(spl[1])
    chr = spl[2]
    pos = int(spl[3])
    rc = 0
    if flag & 0x10:
      pos += len(spl[9]) + parseCigar(spl[5])
      rc = 1
    first = 0
    if flag & 0x80:
      first = 1
    return chr, pos, rc, first

  def getTag(lis, tag):
    '''
    Get optional tag from a SAM record.
    '''
    for t in lis:
      spl = t.split(':')
      if spl[0] == tag:
        return spl[-1]
    sys.stderr.write('Error! Cannot find %s in SAM record\n' % tag)
    sys.exit(-1)

  def getHeader(self):
    return self.header

  def getInfo(self):
    res = 'first: ' + ' '.join(self.r1) + '\n'
    res = 'second: ' + ' '.join(self.r2)
    return res

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
          print read.getHeader() + read.getInfo()
          break
      headers[spl[0]] = 1
    else:
      reads.append(Read(spl))

    #for read in reads:
    #  print read.getInfo()
    #raw_input()
    continue

    # determine location -- use 3' end if RC
    if chr in pos:
      if loc in pos[chr]:
        if rc in pos[chr][loc]:
          if spl[9] in pos[chr][loc][rc]:
            # location, strand, *and* seq match
            dups += 1
          else:
            # everything matches except seq
            pos[chr][loc][rc].append(spl[9])
            if out: fOut.write(line)
            notSeq += 1
        else:
          pos[chr][loc][rc] = [spl[9]]
          if out: fOut.write(line)
          uniq += 1
      else:
        pos[chr][loc] = {}
        pos[chr][loc][rc] = [spl[9]]
        if out: fOut.write(line)
        uniq += 1
    else:
      pos[chr] = {}
      pos[chr][loc] = {}
      pos[chr][loc][rc] = [spl[9]]
      if out: fOut.write(line)
      uniq += 1
    count += 1
  f.close()
  if out: fOut.close()

  sys.stderr.write('Reads: %10d\n' % count)
  sys.stderr.write('Unique: %9d\n' % uniq)
  sys.stderr.write('Dups: %11d\n' % dups)
  sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

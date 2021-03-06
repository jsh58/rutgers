#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: producing a list of PCR duplicates

import sys
import gzip
import re

class Align():
  '''
  A class for representing an alignment.
  '''
  def __init__(self, chrom, pos, revComp, offset, secondary, \
      hi, ascore):
    self.chrom = chrom
    self.pos = pos
    self.revComp = revComp  # may be redundant, since offset will be 0 if revComp == 0
    self.offset = offset  # dist. to 3' end for rc alignments
    self.secondary = secondary
    self.hi = hi
    self.ascore = ascore

  def __str__(self):
    return ' '.join([self.chrom, str(self.pos), str(self.revComp)])

class PEalign(Align):
  '''
  A class for representing an alignment for a
    paired-end read.
  '''
  def __init__(self, chrom, pos, revComp, offset, secondary, \
      hi, ascore, mateChr, matePos):
    Align.__init__(self, chrom, pos, revComp, offset, secondary, hi, ascore)
    self.mateChr = mateChr
    self.matePos = matePos

#class PEalign():
#  def __init__(self, first, aln):
#    if first:
#      self.r2 = aln
#      self.r1 = None
#    else:
#      self.r1 = aln
#      self.r2 = None
#    self.paired = 0  # alignment is full

class Read():
  '''
  A class for representing the position(s) of
    a read's alignment(s).
  '''
  def __init__(self):
    #self.aln = []   # all alignments
    #self.pair = []  # paired alignments

    self.r1 = []    # alignments of first read
    self.r2 = []    # alignments of second read
    #self.header = spl[0]
    #self.addInfo(spl)  # save initial alignment


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

################################################################

def parseCigar(cigar):
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
      sys.stderr.write('Error! Unknown Op "%s" in cigar %s\n' \
        % (op[1], cigar))
      sys.exit(-1)
  return offset

def getTag(lis, tag):
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

def loadInfo(spl):
  '''
  Create an Align object from a SAM record.
  '''
  flag = int(spl[1])
  chrom = spl[2]
  pos = int(spl[3])

  # deal with alignments to the reverse strand
  revComp = 0
  offset = 0  # distance to 3' end of alignment
  if flag & 0x10:
    revComp = 1  # may be redundant since offset will be > 0
    offset = parseCigar(spl[5])

  first = 0   # 0 -> first segment in template; 1 -> last segment
  if flag & 0x80:
    first = 1
    if flag & 0x40:
      sys.stderr.write('Error! Linear templates with > 2 reads '
        + 'are not allowed\n')
      sys.exit(-1)
  secondary = flag & 0x900

  hi = getTag(spl[11:], 'HI')  # query hit index
  ascore = getTag(spl[11:], 'AS')  # alignment score

  # return Align object
  if flag & 0x1 and not flag & 0x8:
    # paired alignment: save mate information
    matePos = int(spl[7])
    mateChr = ''
    if spl[6] != '*' and matePos != 0:
      mateChr = spl[2]
      if spl[6] != '=':
        mateChr = spl[6]
    return PEalign(chrom, pos, revComp, offset, secondary, hi, ascore, mateChr, matePos), first

  return Align(chrom, pos, revComp, offset, secondary, hi, ascore), first

#def addInfo(spl, read):
#  '''
#  Add alignment info to existing read.
#  '''
#  aln, first = loadInfo(spl)
#  if first:
#    read.r2.append(aln)
#  else:
#    read.r1.append(aln)
#
#    flag = 1
#    for r in self.r2:
#      # determine if paired alignment
#      if r.getPos() == aln.getMatePos() \
#          and aln.getPos() == r.getMatePos():
#        self.pair.append((r, aln))
#        self.r2.remove(r)
#        flag = 0
#        break
#    if flag:
#      self.r1.append(aln)



def processSAM(f, fOut):
  '''
  Process the SAM file.
  '''
  reads = {}    # dict of Reads
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

    # skip if unmapped
    if int(spl[1]) & 0x4:
      continue

    # save alignment info
    if spl[0] not in reads:
      reads[spl[0]] = Read()
    aln, first = loadInfo(spl)
    if first:
      reads[spl[0]].r2.append(aln)
    else:
      reads[spl[0]].r1.append(aln)

    #for read in reads:
    #  print read.getInfo()
    #raw_input()
    count += 1
    if count % 100000 == 0:
      print count
    #if count % 100000 == 0:
    #  break

  print count
  for r in reads:
    print 'first:', len(reads[r].r1)
    print 'second:', len(reads[r].r2)
    for aln in reads[r].r1:
      print r, aln
    #for aln in reads[r].r2:
    #  print r, aln

  sys.exit(0)
  # check alignment scores
  for header in reads:
    r1, r2 = reads[header].getInfo()
    best = 1000
    for aln in r1:
      spl = aln.split(':')
      if best == 1000:
        best = spl[-1]
      elif spl[-1] != best:
        print header, r1

    best = 1000
    for aln in r2:
      spl = aln.split(':')
      if best == 1000:
        best = spl[-1]
      elif spl[-1] != best:
        print header, r2

  sys.exit(0)

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

  processSAM(f, fOut)

  f.close()
  if out: fOut.close()

  sys.stderr.write('Reads: %10d\n' % count)
  sys.stderr.write('Unique: %9d\n' % uniq)
  sys.stderr.write('Dups: %11d\n' % dups)
  sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

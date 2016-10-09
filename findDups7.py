#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: producing a list of PCR duplicates
# Version 7: specific for TopHat's output file

import sys
import gzip
import re

class Align():
  '''
  A class for representing an alignment.
  '''
  def __init__(self, chrom, pos, revComp, offset, secondary, \
      hi, ascore):
    self.chrom5 = chrom
    self.chrom3 = chrom

    self.pos5 = pos
    self.pos3 = pos
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
    self.aln = {}   # all alignments


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
  if flag & 0x10:
    revComp = 1  # may be redundant since offset will be > 0
    pos += parseCigar(spl[5])

  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1

  return (chrom, pos), first, revComp




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
    aln, first, revComp = loadInfo(spl)
    if spl[0] not in reads:
      reads[spl[0]] = [{} for i in range(2)]

    hi = getTag(spl[11:], 'HI')  # query hit index
    if hi not in reads[spl[0]][first]:
      reads[spl[0]][first][hi] = [[] for i in range(2)]
    if reads[spl[0]][first][hi][revComp]:
      sys.stderr.write('Error! already a result for %s, %d, %d, %d\n' % (spl[0], first, hi, revComp))
      sys.exit(-1)
    reads[spl[0]][first][hi][revComp] = aln

    count += 1
    if count % 100000 == 0:
      print count

  print 'Total alignments:', count
  print 'Unique reads:', len(reads)

  # conglomerate alignments
  reads5 = {}
  reads3 = {}
  readsPE = {}
  for r in reads:
    for i in range(2):
      for h in reads[r][i]:
        if reads[r][i][h][0] and reads[r][i][h][1]:
          if r not in readsPE:
            readsPE[r] = []
          readsPE[r].append(reads[r][i][h][0] + reads[r][i][h][1])
          if r in reads5:
            reads5.remove(r)
          if r in reads3:
            reads3.remove(r)
        elif r not in readsPE:
          if not reads[r][i][h][1]:
            if r not in reads5:
              reads5[r] = []
            reads5[r].append(reads[r][i][h][0])
          else:
            if r not in reads3:
              reads3[r] = []
            reads3[r].append(reads[r][i][h][1])

  dups5 = []  # saving duplicate read headers
  dups3 = []  # ditto
  dupsPE = [] # ditto
  rep = {}    # saving alignments, to compare against
  for r in reads5:
    for aln in reads5[r]:
      if aln in rep:
        dups5.append(r)
      else:
        rep[aln] = 1
  rep = {}
  for r in reads3:
    for aln in reads3[r]:
      if aln in rep:
        dups3.append(r)
      else:
        rep[aln] = 1
  rep = {}
  for r in readsPE:
    for aln in readsPE[r]:
      if aln in rep:
        dupsPE.append(r)
      else:
        rep[aln] = 1

  print 'Reads PE:', len(readsPE)
  print '    Dups:', len(dupsPE)
  print 'Reads 5\':', len(reads5)
  print '    Dups:', len(dups5)
  print 'Reads 3\':', len(reads3)
  print '    Dups:', len(dups3)

  for r in dupsPE:
    fOut.write(r + '\n')
  for r in dups5:
    fOut.write(r + '\n')
  for r in dups3:
    fOut.write(r + '\n')


def main():
  args = sys.argv[1:]
  if len(args) < 2:
    print 'Usage: python %s  <SAMfile>  <out>' % sys.argv[0]
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

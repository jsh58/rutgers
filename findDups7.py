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
  offset = 0
  if flag & 0x10:
    revComp = 1  # may be redundant since offset will be > 0
    offset = parseCigar(spl[5])

  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1

  paired = 1
  if flag & 0x8:
    paired = 0

  return (chrom, pos, offset), first, paired




def processSAM(f, fOut):
  '''
  Process the SAM file.
  '''
  readsPE = {}    # dict of paired-end reads
  readsSE = {}    # dict of single-end reads
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
    raw, first, paired = loadInfo(spl)
    # reconfigure alignment
    if raw[2] > 0:
      aln = (raw[0], raw[1] + raw[2], '-')
    else:
      aln = (raw[0], raw[1], '+')

    # only one end aligned: save to readsSE
    if not paired:
      if spl[0] not in readsSE:
        readsSE[spl[0]] = []
      readsSE[spl[0]].append(aln)

    # paired-end (even if not properly paired): save to readsPE
    else:
      if spl[0] not in readsPE:
        readsPE[spl[0]] = {}
      hi = getTag(spl[11:], 'HI')  # query hit index
      if hi not in readsPE[spl[0]]:
        readsPE[spl[0]][hi] = [() for i in range(2)]
      if aln[2] == '+':
        readsPE[spl[0]][hi][0] = aln  # put '+' alignment 1st
      else:
        readsPE[spl[0]][hi][1] = aln  # '-' alignment 2nd

    #if reads[spl[0]][first][hi][revComp]:
    #  sys.stderr.write('Error! already a result for %s, %d, %d, %d\n' % (spl[0], first, hi, revComp))
    #  sys.exit(-1)
    #reads[spl[0]][first][hi][revComp] = aln

    count += 1
    if count % 1000000 == 0:
      print count

  print 'Total alignments:', count

  # find dups in SE alignments
  print 'SE reads:', len(readsSE)
  rep = {}    # saving alignments, to compare against
  dups = 0
  for r in readsSE:
    printed = 1
    for aln in readsSE[r]:
      if aln in rep:
        if printed:
          fOut.write(r + '\n')
          dups += 1
          printed = 0   # don't write dups twice!
      else:
        rep[aln] = 1
  print '  dups:', dups


  # find dups in PE alignments
  print 'PE reads:', len(readsPE)
  dups = 0
  rep = {}    # saving alignments, to compare against
  for r in readsPE:
    printed = 1
    for h in readsPE[r]:
      aln0, aln1 = readsPE[r][h]
      aln = aln0 + aln1
      if aln in rep:
        if printed:
          fOut.write(r + '\n')
          dups += 1
          printed = 0   # don't write dups twice!
      else:
        rep[aln] = 1
  print '  dups:', dups


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

  #sys.stderr.write('Reads: %10d\n' % count)
  #sys.stderr.write('Unique: %9d\n' % uniq)
  #sys.stderr.write('Dups: %11d\n' % dups)
  #sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

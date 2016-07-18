#!/usr/bin/python

# JMG 6/13/16
# Retrieving a subset of reads from a SAM.
# Version 3: condensing just CpG methylation data.

import sys
import re

class Read():
  '''
  A class for representing a read's position,
    seq, qual, and methylation string.
  '''
  def __init__(self, header, chr, pos, cigar, seq, qual, meth):
    self.header = header
    self.chr = chr
    self.pos = pos
    self.seq = seq
    self.qual = qual
    self.meth = meth

    # convert cigar to 'M/I/D' string
    cig = ''
    for c in cigar:
      dist = int(c[0])
      cig += c[1] * dist
    self.cigar = cig

  def getPos(self):
    return self.pos

  def getCigar(self):
    return self.cigar

  def setCigar(self, cigar):
    self.cigar = cigar

  def setMeth(self, meth):
    self.meth = meth

  def adjustMeth(self):
    '''
    Remove non-CpG methylation data from methylation string.
    #Convert CpG methylation data to '/' if quality is low.
    '''
    meth = self.meth
    for c in ['h', 'H', 'x', 'X']:
      meth = meth.replace(c, '.')

    # remove low quality CpG bases -- min qual 30
    for c in ['z', 'Z']:
      idx = meth.find(c)
      while idx != -1:
        # convert low qual bases (< 30) to '/'
        #if ord(self.qual[idx]) - 33 < 30:
        #  meth = meth[:idx] + '/' + meth[idx+1:]
        #print c, self.qual[idx], ord(self.qual[idx]), '\n'
        idx = meth.find(c, idx + 1)
    self.meth = meth

  def adjustCig(self, minLoc):
    '''
    Adjust the cigar -- add gaps to 5' end.
    '''
    # add gaps to 5' end
    offset = self.pos - minLoc
    self.cigar = 'G' * offset + self.cigar

  def getMeth(self):
    '''
    Construct the methylation string,
      including gaps.
    '''
    ans = ''
    loc = 0  # location in self.meth
    for i in range(len(self.cigar)):
      # add bases from methylation string
      c = self.cigar[i]
      if c == 'D':
        ans += '-'
      elif c == 'G':
        ans += ' '
      elif c in ['M', 'I']:
        ans += self.meth[loc]
        loc += 1
      else:
        sys.stderr.write('Error! Unknown CIGAR base: %s\n' % c)
        sys.exit(-1)
    return ans


def addIns(reads):
  '''
  Add inserted bases to the cigars.
  '''
  # create matrix of cigars
  mat = []
  for read in reads:
    mat.append(read.getCigar())

  # check each position for 'I'
  i = 0
  while i < len(mat[-1]):
    for cig in mat:
      if i < len(cig) and cig[i] == 'I':
        for j in range(len(mat)):
          if i < len(mat[j]) and mat[j][i] != 'I':
            if i == 0 or mat[j][i-1] == 'G':
              mat[j] = mat[j][:i] + 'G' + mat[j][i:]
            else:
              mat[j] = mat[j][:i] + 'D' + mat[j][i:]
        break
    i += 1

  # save adjusted cigars
  for i in range(len(reads)):
    reads[i].setCigar(mat[i])


def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to int\n' % arg)
    sys.exit(-1)
  return val

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
  '''
  if filename == '-':
    return sys.stdin
  try:
    f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
  '''
  if filename == '-':
    return sys.stdout
  try:
    f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

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

def loadGen(fGen, chrom):
  '''
  Load a chromosome from the given genome file.
  '''
  # find chromosome
  f = openRead(fGen)
  seq = ''
  for line in f:
    if line.rstrip()[1:] == chrom:
      # save sequence
      for line in f:
        if line[0] == '>':
          return seq
        seq += line.rstrip()
  if seq:
    return seq
  sys.stderr.write('Error! Cannot find chromosome %s\n' % chrom)
  sys.exit(-1)

def parseSAM(f, chrom, start, end, reads):
  '''
  Parse the SAM file. Select reads overlapping the
    given coordinates.
  '''
  minLoc = 1000000000
  maxLoc = 0
  count = total = 0
  for line in f:
    if line[0] == '@':
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
      reads.append(Read(spl[0], chr, loc, cigar, \
        spl[9], spl[10], getXM(spl[11:])))
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
  if len(args) < 6:
    print 'Usage: python %s  <SAMfile>  <out> ' % sys.argv[0],
    print '<chr>  <start>  <end>  <genome>'
    print '  Use \'-\' for stdin/stdout'
    sys.exit(-1)

  # open files
  f = openRead(args[0])
  fOut = openWrite(args[1])

  # save region of interest
  chr = args[2]
  start = getInt(args[3])
  end = getInt(args[4])

  # process file
  reads = []
  count, total, minLoc, maxLoc = parseSAM(f, chr, \
    start, end, reads)
  reads = sorted(reads, key=Read.getPos)  # sort by position
  f.close()

  # get genomic segment
  if len(args) > 5:
    seq = loadGen(args[5], chr)
    # prepend to reads list
    reads.insert(0, Read('genome', chr, minLoc-1, [(str(maxLoc - minLoc), 'M')], \
      seq[minLoc-1:maxLoc-1].upper(), '', seq[minLoc-1:maxLoc-1].upper()))

  # adjust reads for 5' ends, insertions
  for read in reads:
    read.adjustCig(minLoc)
    #read.adjustMeth()
  addIns(reads)

  # find methylation loci
  loci = []
  gap = []  # running tally of deleted bases
  tally = 0
  gMeth = reads[0].getMeth()
  for i in range(len(gMeth)):
    for read in reads[1:]:
      meth = read.getMeth()
      if i < len(meth) and meth[i] in ['z', 'Z']:
        loci.append(i)
        break
    if gMeth[i] == '-':
      tally += 1
    gap.append(tally)

  # print output
  #fOut.write('%s\t1-based positions:\n' % chr)
  summary = {}
  for read in reads[1:]:
    meth = read.getMeth()
    res = ''  # summarized methylation string
    idx = 0
    while idx < len(loci):
      pos = (loci[idx],)
      chrPos = str(pos[0] + minLoc - gap[pos[0]])  # actual position in genome
      if idx + 1 < len(loci) and loci[idx] + 1 == loci[idx + 1]:
        pos += (loci[idx + 1],)
        chrPos += '+1'
        idx += 1
      #if read == reads[1]:
      #  fOut.write('%d' % (chrPos))
      #  if len(pos) > 1:
      #    fOut.write('+1')
      #  fOut.write('\t')
      if not chrPos in summary:
        summary[chrPos] = [0, 0]
      if pos[0] >= len(meth) or \
          (len(pos) == 1 and meth[pos[0]] == ' ') or \
          (len(pos) > 1 and pos[1] < len(meth) and meth[pos[1]] == ' '):
        # no data: labeled '-'
        res += '-'
      else:
        flag = 1
        for i in pos:
          if i >= len(meth):
            # no data: labeled '-'
            res += '-'
            flag = 0
            break
          if meth[i] == 'z':
            # unmethylated: labeled '0'
            res += '0'
            flag = 0
            summary[chrPos][0] += 1
            break
          if meth[i] == 'Z':
            # methylated: labeled '1'
            res += '1'
            flag = 0
            summary[chrPos][1] += 1
            break
        if flag:
          # other labeled: 'x'
          res += 'x'
      idx += 1
    #if read == reads[1]:
    #  fOut.write('\n')
    fOut.write('%s\n' % res)

  fOut.write('\n' + '=' * 50 + '\n\n')
  for i in sorted(summary):
    # bismark output lists methylated before unmethylated
    fOut.write('%s\t%s\t%d\t%d\n' % (chr, i, summary[i][1], summary[i][0]))
  fOut.close()

  sys.stderr.write('Reads analyzed: %d\n' % count)
  sys.stderr.write('Reads written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()

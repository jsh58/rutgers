#!/usr/bin/python

# JMG 7/21/16
# Producing methylation summary information
#   directly from a SAM file.

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

def getTag(lis, tag):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  sys.stderr.write('Error! Cannot find ' + tag + ' in SAM record\n')
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

def parseSAM(f, d):
  '''
  Parse the SAM file. Select reads overlapping the
    given coordinates.
  '''
  genome = []  # ordered chromosome names
  minLoc = 1000000000
  maxLoc = 0
  count = total = 0
  for line in f:

    # save chromosome names
    if line[0] == '@':
      spl = line.rstrip().split('\t')[1].split(':')
      if line[1:3] == 'SQ' and spl[0] == 'SN':
        d[spl[1]] = {}
        genome.append(spl[1])
      continue

    spl = line.rstrip().split('\t')
    chrom = spl[2]
    pos = getInt(spl[3])
    offset, cigar = parseCigar(spl[5])
    meth = getTag(spl[11:], 'XM')

    rc = 0
    if getInt(spl[1]) & 0x10:
      rc = 1
    if chrom not in d:
      sys.stderr.write('Error! Cannot find chromosome %s in genome\n' % chrom)
      sys.exit(-1)

    for i in range(len(meth)):
      if meth[i] in ['z', 'Z']:
        loc = pos + i - rc  # need to adjust based on cigar I/D
        if loc not in d[chrom]:
          d[chrom][loc] = [0, 0]
        if meth[i] == 'z':
          d[chrom][loc][1] += 1
        else:
          d[chrom][loc][0] += 1

    count += 1

  sys.stderr.write('Reads analyzed: ' + str(count) + '\n')
  return genome

def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 2:
    print 'Usage: python %s  <SAMfile>  <out> ' % sys.argv[0]#,
    #print '<chr>  <start>  <end>  [<genome>]'
    print '  Use \'-\' for stdin'
    sys.exit(-1)

  # open files
  f = openRead(args[0])
  fOut = openWrite(args[1])

  # save region of interest
  #chr = args[2]
  #start = getInt(args[3])
  #end = getInt(args[4])

  # process file
  d = {}
  genome = parseSAM(f, d)
  f.close()

  # print output
  for chrom in genome:
    for loc in sorted(d[chrom]):
      print chrom, loc, d[chrom][loc]
    raw_input()



  for chrom in genome:
    fOut.write('>%s %d %d\n' % (chr, minLoc, maxLoc))
    fOut.write(read.getMeth() + '\n')
    #fOut.write(read.getCigar() + '\n')
  fOut.close()

  sys.stderr.write('Reads analyzed: %d\n' % count)
  sys.stderr.write('Reads written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()

#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: producing a list of PCR duplicates
# Version 7: specific for TopHat's output file
# Version 8: using quality scores to rank reads

import sys
import gzip
import re

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

def getScore(qual):
  '''
  Return sum of quality scores.
  '''
  qualSum = 0
  for c in qual:
    qualSum += ord(c) - 33
  return qualSum

def parseCigarFwd(cigar):
  '''
  Determine splice sites from CIGAR.
  '''
  ops = re.findall(r'(\d+)(\D)', cigar)
  offset = 0
  ntup = ()  # tuple of spliced lengths
  for op in ops:
    if op[1] in ['M', 'D', 'N', 'S', '=', 'X']:
      if op[1] == 'N':
        ntup += (str(offset) + '-' + op[0],)
      offset += int(op[0])
    elif op[1] in ['I', 'P', 'H']:
      pass
    else:
      sys.stderr.write('Error! Unknown Op "%s" in cigar %s\n' \
        % (op[1], cigar))
      sys.exit(-1)
  return ntup

def parseCigarRev(cigar):
  '''
  Determine splice sites and 3' position of
    alignment from CIGAR.
  '''
  ops = re.findall(r'(\d+)(\D)', cigar)
  offset = 0
  ntup = ()  # tuple of spliced lengths
  # parse cigar backwards
  for op in ops[::-1]:
    if op[1] in ['M', 'D', 'N', 'S', '=', 'X']:
      if op[1] == 'N':
        ntup += (str(offset) + '-' + op[0],)
      offset += int(op[0])
    elif op[1] in ['I', 'P', 'H']:
      pass
    else:
      sys.stderr.write('Error! Unknown Op "%s" in cigar %s\n' \
        % (op[1], cigar))
      sys.exit(-1)
  return offset, ntup

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
  Load alignment info from a SAM record.
  '''
  flag = int(spl[1])
  chrom = spl[2]
  pos = int(spl[3])

  # get tuple of splice sites (pos-len pairs) and
  #   offset (distance to 3' end) for rev-comp alignments
  if flag & 0x10:
    revComp = 1
    offset, ntup = parseCigarRev(spl[5])
  else:
    revComp = 0  # may be redundant since offset will be 0
    offset = 0
    ntup = parseCigarFwd(spl[5])

  # first or last segment?
  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1

  # paired alignments?
  paired = 1
  if flag & 0x8:
    paired = 0

  # primary alignment?
  primary = 0   # 0 -> secondary/supplementary alignment
                # 1 -> primary alignment
  if not flag & 0x900:
    primary = 1

  # ignoring 'first'
  return (chrom, pos, offset, ntup), paired, primary


def findDupsSE(fOut, readsSE, repSE, scoreSE):
  '''
  Find duplicates in single-end reads.
  '''
  # repSE dict already has single-end versions of PE alignments
  print 'SE reads:%10d' % len(readsSE)
  dups = 0
  for r in sorted(scoreSE, key=scoreSE.get, reverse=True):
    printed = 0  # boolean: 1 -> read classified a duplicate
    for aln in readsSE[r]:
      if aln in repSE:
        if not printed and repSE[aln] != r:
          fOut.write('\t'.join([r, 'aln:' + str(aln), repSE[aln], 'single-end\n']))
          dups += 1
          printed = 1   # don't write dups twice!
      else:
        repSE[aln] = r
  print '  dups:%12d' % dups


def findDups(fOut, readsSE, readsPE, scoreSE, scorePE):
  '''
  Find duplicates in paired-end alignments.
  '''
  print 'PE reads:%10d' % len(readsPE)
  dups = 0
  rep = {}    # for saving PE alignments (to compare against)
  repSE = {}  # for saving each SE alignment
  for r in sorted(scorePE, key=scorePE.get, reverse=True):
    printed = 0  # boolean: 1 -> read classified a duplicate
    for h in readsPE[r]:
      if not readsPE[r][h][0] or not readsPE[r][h][1]:
        sys.stderr.write('Error! Incomplete PE alignment for read %s, HI %d\n' % (r, h))
        sys.exit(-1)
      aln0, aln1 = readsPE[r][h]
      if (aln0, aln1) in rep:
        if not printed and rep[(aln0, aln1)] != r:
          fOut.write('\t'.join([r, 'aln:' + str((aln0, aln1)), rep[(aln0, aln1)], 'paired-end\n']))
          dups += 1
          printed = 1   # don't write dups twice!
      elif (aln1, aln0) in rep:
        if not printed and rep[(aln1, aln0)] != r:
          fOut.write('\t'.join([r, 'aln:' + str((aln1, aln0)), rep[(aln1, aln0)], 'paired-end\n']))
          dups += 1
          printed = 1   # don't write dups twice!
      else:
        rep[(aln0, aln1)] = r
        if aln0 not in repSE:
          repSE[aln0] = r
        if aln1 not in repSE:
          repSE[aln1] = r
  print '  dups:%12d' % dups

  findDupsSE(fOut, readsSE, repSE, scoreSE)


def processSAM(f, fOut):
  '''
  Process the SAM file.
  '''
  readsPE = {}    # dict of paired-end alignments
  readsSE = {}    # dict of single-end alignments
  scorePE = {}    # dict of quality score sum for each pair of reads
  scoreSE = {}    # dict of quality score sum for each read
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
    raw, paired, primary = loadInfo(spl)
    # reconfigure alignment -- for revComp alignments, move
    #   pos to 3' end (on fwd strand)
    if raw[2] > 0:
      aln = (raw[0], raw[1] + raw[2], raw[3], '-')
    else:
      aln = (raw[0], raw[1], raw[3], '+')

    # only one end aligned: save to readsSE
    if not paired:
      if spl[0] not in readsSE:
        readsSE[spl[0]] = []
      readsSE[spl[0]].append(aln)
      if primary:
        scoreSE[spl[0]] = getScore(spl[10])

    # paired-end (even if not properly paired): save to readsPE
    else:
      if spl[0] not in readsPE:
        readsPE[spl[0]] = {}
      hi = getTag(spl[11:], 'HI')  # query hit index
      # NOTE: if 'hi' is repeated -- e.g. if it is == -1 (not present), then
      #   can consider using RNEXT and PNEXT as the surrogate 'hi' value
      #   -> but then, must use native alignment info, i.e. for '-' alignments,
      #     5' POS (raw[1]) instead of aln[1]; these alignments may be repeated,
      #     but at least primary alignments can be matched (! 0x100)

      if hi not in readsPE[spl[0]]:
        readsPE[spl[0]][hi] = [() for i in range(2)]
        readsPE[spl[0]][hi][0] = aln  # just put this alignment first --
                                      #   can put '+' alignment 1st by checking aln[3]
      elif readsPE[spl[0]][hi][1]:
        sys.stderr.write('Error! More than two alignments for ' \
          + 'read %s, HI %s\n' % (spl[0], hi))
        sys.stderr.write('  alns: ' + str(readsPE[spl[0]][hi]) \
          + ',' + str(aln) + '\n')
        sys.exit(-1)
      else:
        readsPE[spl[0]][hi][1] = aln  # 2nd alignment

      if primary:
        scorePE[spl[0]] = scorePE.get(spl[0], 0) + getScore(spl[10])

    count += 1
    #if count % 10000000 == 0:
    #  print count

  print 'Total alignments:', count

  findDups(fOut, readsSE, readsPE, scoreSE, scorePE)


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

  if f != sys.stdin:
    f.close()
  if out and fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()

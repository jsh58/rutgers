#!/usr/bin/python

# JMG 5/29/16 (updated 10/26/16)

# Given a SAM file, this script classifies reads that are
#   potential PCR duplicates.
# The determination is based on each alignment's reference
#   name (RNAME), position (POS), position and length of
#   gaps (e.g. splice sites, defined by the CIGAR 'N' Op),
#   orientation ('+' or '-'), and status (paired or not).
# A paired alignment can be a duplicate only of other
#   paired alignments, with both alignments matching.
#   A single alignment can be a duplicate of a paired
#   alignment (if it matches either end) or another single
#   alignment.  An alignment is considered 'single' if the
#   0x1 bit of the FLAG is not set, or if the 0x8 bit is
#   set.
# Secondary alignments are analyzed the same as primary
#   alignments.  A read is considered a duplicate if any
#   of its multiple alignments is classified a duplicate.
# The output file lists, for each duplicate, the read
#   header (QNAME), the alignment that matched another,
#   the 'parent' read header, and whether the match was
#   single or paired, all tab-delimited.
# Note that this script does not remove duplicates; it
#   just finds them.  The script getReads.py will filter
#   the duplicate reads from a SAM file.  So, a complete
#   find-and-remove process would go like this:
# $ samtools view <BAM> | python findDups.py - dups.txt;
#   samtools view -h <BAM> | python getReads.py - no dups.txt - | \
#   samtools view -bS - > <outBAM>

# SE version: treat all reads like SE

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
  Determine start position and splice sites from CIGAR.
  '''
  ops = re.findall(r'(\d+)(\D)', cigar)

  # determine actual start position of sequence
  #   (adjusting backwards for 'S', forwards for 'D')
  start = 0
  if ops[0][1] == 'S':
    start = - int(ops[0][0])
  elif ops[0][1] == 'D':
    start = int(ops[0][0])

  # parse cigar
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
  return start, ntup

def parseCigarRev(cigar):
  '''
  Determine start position, splice sites, and
    3' position of alignment from CIGAR.
  '''
  ops = re.findall(r'(\d+)(\D)', cigar)

  # determine actual start position of sequence
  #   (adjusting backwards for 'S', forwards for 'D')
  start = 0
  if ops[0][1] == 'S':
    start = - int(ops[0][0])
  elif ops[0][1] == 'D':
    start = int(ops[0][0])

  # parse cigar backwards
  offset = 0
  ntup = ()  # tuple of spliced lengths
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
  return start, offset, ntup

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

  # get start site, tuple of splice sites (pos-len pairs),
  #   and (for rev-comp alignments) distance to 3' end
  if flag & 0x10:
    revComp = 1
    start, offset, ntup = parseCigarRev(spl[5])
  else:
    revComp = 0  # may be redundant since offset will be 0
    offset = 0
    start, ntup = parseCigarFwd(spl[5])

  # first or last segment?
  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1
    if flag & 0x40:
      sys.stderr.write('Error! Linear templates with > 2 reads ' \
        + 'are not allowed\n')
      sys.exit(-1)

  # paired alignments?
  paired = 1
  if not flag & 0x1 or flag & 0x8:
    paired = 0

  # primary alignment?
  primary = 0   # 0 -> secondary/supplementary alignment
                # 1 -> primary alignment
  if not flag & 0x900:
    primary = 1

  # return alignment info
  return (chrom, pos, offset, ntup), start, revComp, \
    first, paired, primary

def processNoHI(header, first, readsPE, noHI, aln, mate):
  '''
  In the absence of an 'HI' tag, try to find a matching PE alignment.
    'aln' is this alignment (RNAME, POS);
    'mate' is the mate's alignment (RNEXT, PNEXT).
  '''
  if mate[0] == '=':
    mate[0] = aln[0]
  nextAln = (mate[0], int(mate[1]))  # tuple of mate's alignment
  if nextAln[0] == '*' or nextAln[1] == 0:
    sys.stderr.write('Error! Missing paired alignment ' \
      + 'for read %s\n' % header)
    sys.exit(-1)
  other = (first + 1) % 2  # mate's 'first' designation

  # try to find a matching alignment
  hi = 1  # match index (0 reserved for primary alignments)
  match = -1  # alternative match index
  while hi in readsPE[header]:
    if header in noHI[other] \
        and hi in noHI[other][header] \
        and (nextAln, aln) == noHI[other][header][hi]:
      # alignments match!
      if readsPE[header][hi][first]:
        match = hi  # if already has alignment, save index and continue
      else:
        return hi
    hi += 1

  # for alternative match, copy readsPE alignment for new index
  if match != -1:
    readsPE[header][hi] = [() for i in range(2)]
    readsPE[header][hi][other] = readsPE[header][match][other]

  # no matches -- save these alignments to noHI
  if header not in noHI[first]:
    noHI[first][header] = {}
  noHI[first][header][hi] = (aln, nextAln)
  return hi

def addMissing(readsPE, noHI):
  '''
  For PE alignments without 'HI' tags that are not matched,
    fill in each missing alignment with the corresponding
    primary alignment.
  '''
  for first in range(2):
    other = (first + 1) % 2
    for header in noHI[first]:
      for hi in noHI[first][header]:
        if not readsPE[header][hi][other]:
          # missing alignment: copy primary (hi=0)
          readsPE[header][hi][other] = readsPE[header][0][other]

def findDupsSE(fOut, readsSE, repSE, scoreSE):
  '''
  Find duplicates in single-end reads.
  '''
  # repSE dict already has single-end versions of PE alignments
  print 'SE-aligned reads:%10d' % len(readsSE)
  dups = 0
  for r in sorted(scoreSE, key=scoreSE.get, reverse=True):
    printed = 0  # boolean: 1 -> read classified a duplicate
    for aln in readsSE[r]:
      if aln in repSE:
        # duplicate
        fOut.write('\t'.join([r, 'aln:' + str(aln), \
          repSE[aln], 'single-end\n']))
        dups += 1
        printed = 1
        break

    # save non-dup alignments
    if not printed:
      for aln in readsSE[r]:
        repSE[aln] = r

  print '  duplicates:%14d' % dups

def findDups(fOut, readsSE, readsPE, scoreSE, scorePE):
  '''
  Find duplicates in paired-end alignments.
  '''
  print 'PE-aligned reads:%10d' % len(readsPE)
  dups = 0
  rep = {}    # for saving PE alignments (to compare against)
  repSE = {}  # for saving each SE alignment

  # check for duplicates, in descending order by quality scores
  for r in sorted(scorePE, key=scorePE.get, reverse=True):
    printed = 0  # boolean: 1 -> read classified a duplicate
    for h in readsPE[r]:
      if not readsPE[r][h][0] or not readsPE[r][h][1]:
        sys.stderr.write('Error! Incomplete PE alignment for ' \
          + 'read %s, HI %d\n' % (r, h))
        sys.exit(-1)
      aln0, aln1 = readsPE[r][h]
      if (aln0, aln1) in rep:
        # duplicate
        fOut.write('\t'.join([r, 'aln:' + str((aln0, aln1)), \
          rep[(aln0, aln1)], 'paired-end\n']))
        dups += 1
        printed = 1
        break
      elif (aln1, aln0) in rep:
        # duplicate, in reverse order
        fOut.write('\t'.join([r, 'aln:' + str((aln1, aln0)), \
          rep[(aln1, aln0)], 'paired-end\n']))
        dups += 1
        printed = 1
        break

    # save non-dup alignments
    if not printed:
      for h in readsPE[r]:
        aln0, aln1 = readsPE[r][h]
        rep[(aln0, aln1)] = r
        if aln0 not in repSE:
          repSE[aln0] = r
        if aln1 not in repSE:
          repSE[aln1] = r

  print '  duplicates:%14d' % dups

  # find single-end duplicates
  findDupsSE(fOut, readsSE, repSE, scoreSE)

def processSAM(f, fOut):
  '''
  Process the SAM file.
  '''
  readsPE = {}    # dict of paired-end alignments
  readsSE = {}    # dict of single-end alignments
  scorePE = {}    # dict of quality score sum for each pair of reads
  scoreSE = {}    # dict of quality score sum for each read
  noHI = [{} for i in range(2)] # dicts for matching PE alignments
                                #   lacking 'HI' tags
  count = unmap = 0
  for line in f:
    if line[0] == '@':
      continue

    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n'
        + line)
      sys.exit(-1)
    count += 1

    # skip if unmapped
    if int(spl[1]) & 0x4:
      unmap += 1
      continue

    # save alignment info
    raw, start, revComp, first, paired, primary = loadInfo(spl)
    paired = 0  # treat everything as SE
    # reconfigure alignment -- for revComp alignments, move
    #   pos to 3' end (on fwd strand)
    if revComp:
      aln = (raw[0], raw[1] + start + raw[2], raw[3], '-')
    else:
      aln = (raw[0], raw[1] + start, raw[3], '+')

    # only one end aligned: save to readsSE
    if not paired:
      if spl[0] not in readsSE:
        readsSE[spl[0]] = []
      readsSE[spl[0]].append(aln)
      # save sum of quality scores
      if primary:
        scoreSE[spl[0]] = getScore(spl[10])

    # paired-end (even if not properly paired): save to readsPE
    else:
      if spl[0] not in readsPE:
        readsPE[spl[0]] = {}

      # get query hit index (to match paired alignments)
      hi = getTag(spl[11:], 'HI')
      # if 'hi' unavailable (-1), with multiple alignments,
      #   try to find matching alignments (using noHI dict)
      if hi == -1:
        if primary:
          hi = 0
        else:
          hi = processNoHI(spl[0], first, readsPE, \
            noHI, raw[0:2], spl[6:8])

      # save alignment
      if hi not in readsPE[spl[0]]:
        readsPE[spl[0]][hi] = [() for i in range(2)]
      if readsPE[spl[0]][hi][first]:
        sys.stderr.write('Error! Multiple alignments for ' \
          + 'read %s, R%d, HI %s\n' % (spl[0], first + 1, str(hi)) \
          + '  alns: ' + str(readsPE[spl[0]][hi]) \
          + ', ' + str(aln) + '\n')
        sys.exit(-1)
      readsPE[spl[0]][hi][first] = aln

      # save sum of quality scores
      if primary:
        scorePE[spl[0]] = scorePE.get(spl[0], 0) + getScore(spl[10])

  print 'Total alignments:%10d' % count
  print '  Unmapped:%16d' % unmap

  # add in missing PE alignments
  addMissing(readsPE, noHI)
  # find duplicates
  findDups(fOut, readsSE, readsPE, scoreSE, scorePE)

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <SAMfile>  <out>\n' % sys.argv[0] \
      + '  Use \'-\' for stdin/stdout\n')
    sys.exit(-1)
  fIn = openRead(args[0])
  fOut = openWrite(args[1])

  processSAM(fIn, fOut)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()

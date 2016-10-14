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
  Create an Align object from a SAM record.
  '''
  flag = int(spl[1])
  chrom = spl[2]
  pos = int(spl[3])

  # get offset (distance to 3' end) and tuple of splice sites
  if flag & 0x10:
    revComp = 1
    offset, ntup = parseCigarRev(spl[5])
  else:
    revComp = 0  # may be redundant since offset will be > 0
    offset = 0
    ntup = parseCigarFwd(spl[5])

  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1

  # paired alignments?
  paired = 1
  if flag & 0x8:
    paired = 0

  # ignoring 'first'
  return (chrom, pos, offset, ntup), paired


def findDupsSE(fOut, readsSE, repSE):
  '''
  Find duplicates in single-end reads.
  '''
  # repSE dict already has single-end versions of PE alignments
  print 'SE reads:%10d' % len(readsSE)
  dups = 0
  for r in readsSE:
    printed = 1  # boolean: has read been classified a dup? (don't count twice)
    for aln in readsSE[r]:
      if aln in repSE:
        if printed:
          fOut.write('\t'.join([r, 'aln:' + str(aln), repSE[aln], 'single-end\n']))
          dups += 1
          printed = 0   # don't write dups twice!
      else:
        repSE[aln] = r
  print '  dups:%12d' % dups


def findDups(fOut, readsSE, readsPE):
  '''
  Find duplicates in paired-end alignments.
  '''
  print 'PE reads:%10d' % len(readsPE)
  dups = 0
  rep = {}    # for saving PE alignments (to compare against)
  repSE = {}  # for saving each SE alignment
  for r in readsPE:
    printed = 1  # boolean: has read been classified a dup? (don't count twice)
    for h in readsPE[r]:
      if not readsPE[r][h][0] or not readsPE[r][h][1]:
        sys.stderr.write('Error! Incomplete PE alignment for read %s, HI %d\n' % (r, h))
        sys.exit(-1)
      aln0, aln1 = readsPE[r][h]
      if (aln0, aln1) in rep:
        if printed:
          fOut.write('\t'.join([r, 'aln:' + str((aln0, aln1)), rep[(aln0, aln1)], 'paired-end\n']))
          dups += 1
          printed = 0
      elif (aln1, aln0) in rep:
        if printed:
          fOut.write('\t'.join([r, 'aln:' + str((aln1, aln0)), rep[(aln1, aln0)], 'paired-end\n']))
          dups += 1
          printed = 0
      else:
        rep[(aln0, aln1)] = r
        if aln0 not in repSE:
          repSE[aln0] = r
        if aln1 not in repSE:
          repSE[aln1] = r
  print '  dups:%12d' % dups

  findDupsSE(fOut, readsSE, repSE)


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
    raw, paired = loadInfo(spl)
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
      if aln in readsSE[spl[0]]:
        sys.stderr.write('Error! Repeated alignment for read %s\n' % spl[0])
        sys.stderr.write('  aln: ' + str(aln) + '\n')
        sys.stderr.write(readsSE[spl[0]] + '\n')
        sys.exit(-1)
      readsSE[spl[0]].append(aln)  # only save if aln isn't already there


    # paired-end (even if not properly paired): save to readsPE
    else:
      if spl[0] not in readsPE:
        readsPE[spl[0]] = {}
      hi = getTag(spl[11:], 'HI')  # query hit index
      # if 'hi' is repeated -- e.g. if it is == -1 (not present), then can
      #   consider using RNEXT and PNEXT as the surrogate 'hi' value
      #   -- but then must use native alignment info, not 5' end for '-' alignments

      if hi not in readsPE[spl[0]]:
        readsPE[spl[0]][hi] = [() for i in range(2)]
        readsPE[spl[0]][hi][0] = aln  # just put this alignment first --
                                      #   can put '+' alignment 1st by checking if aln[3] == '+'
      elif readsPE[spl[0]][hi][1]:
        sys.stderr.write('Error! Repeated alignment for read %s, HI %s\n' % (spl[0], hi))
        sys.exit(-1)
      else:
        readsPE[spl[0]][hi][1] = aln  # 2nd alignment

    count += 1
    if count % 10000000 == 0:
      print count

  print 'Total alignments:', count

  findDups(fOut, readsSE, readsPE)


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

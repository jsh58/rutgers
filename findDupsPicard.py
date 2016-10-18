#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 6: producing a list of PCR duplicates
# Version 8: mimicking Picard -- counting duplicates

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
  #strange = 0  # are both alignments in same orientation?
  if flag & 0x10:
    revComp = 1  # may be redundant since offset will be > 0
    offset = parseCigar(spl[5])
    #if flag & 0x20:
    #  strange = 1  # both alignments are '-'
  #elif not flag & 0x20:
  #  strange = 1  # both alignments are '+'

  first = 0   # 0 -> first segment
              # 1 -> last segment
  if flag & 0x80:
    first = 1

  paired = 1
  if flag & 0x8:
    paired = 0

  return (chrom, pos, offset), paired




def processSAM(f, fOut):
  '''
  Process the SAM file.
  '''
  readsPE = {}    # dict of paired-end reads
  readsSE = {}    # dict of single-end reads
  count = second = unmap = dups = uniq = notSeq = 0
  printed = 0
  for line in f:
    if line[0] == '@':
      #if out: fOut.write(line)
      continue
    count += 1

    spl = line.rstrip().split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n'
        + line)
      sys.exit(-1)

    # skip if unmapped
    if int(spl[1]) & 0x4:
      unmap += 1
      continue

    ### SKIP if secondary
    if int(spl[1]) & 0x100:
      second += 1
      printed += 1
      continue

    # save alignment info
    raw, paired = loadInfo(spl)
    # reconfigure alignment
    if raw[2] > 0:
      aln = (raw[0], raw[1] + raw[2], '-')
    else:
      aln = (raw[0], raw[1], '+')

    # only one end aligned: save to readsSE
    if not paired:
      if spl[0] in readsSE:
        sys.stderr.write('Error! Multiple alignments for SE read %s\n' % spl[0])
        sys.exit(-1)
      readsSE[spl[0]] = aln

    # paired-end (even if not properly paired): save to readsPE
    else:
      if spl[0] not in readsPE:
        readsPE[spl[0]] = [() for i in range(2)]
        readsPE[spl[0]][0] = aln
      else:
        readsPE[spl[0]][1] = aln

    count += 1
    if count % 10000000 == 0:
      print count

  print 'Total alignments:', count
  print '  Secondary alignments:', second
  print '  Unmapped alignments:', unmap

  # find dups in PE alignments
  print 'PE reads:%10d' % len(readsPE)
  dups = 0
  rep = {}    # saving PE alignments, to compare against
  repSE = {}  # for single-end alignments
  for r in readsPE:
    if not readsPE[r][0] or not readsPE[r][1]:
      sys.stderr.write('Error! Incomplete PE alignment for read %s\n' % r)
      sys.exit(-1)
    aln0, aln1 = readsPE[r]
    if (aln0, aln1) in rep:
      fOut.write('\t'.join([r, 'aln:' + str((aln0, aln1)), rep[(aln0, aln1)], 'paired-end\n']))
      dups += 1
    elif (aln1, aln0) in rep:
      fOut.write('\t'.join([r, 'aln:' + str((aln1, aln0)), rep[(aln1, aln0)], 'paired-end\n']))
      dups += 1
    else:
      rep[(aln0, aln1)] = r
      printed += 2
      repSE[aln0] = r
      repSE[aln1] = r
  print '  dups:%12d' % dups

  # find dups in SE alignments
  print 'SE reads:%10d' % len(readsSE)
  dups = 0
  for r in readsSE:
    flag = 0
    if readsSE[r] in repSE:
      fOut.write('\t'.join([r, 'aln:' + str(readsSE[r]), repSE[readsSE[r]], 'SE->PE\n']))
      dups += 1
    else:
      printed += 1
      repSE[readsSE[r]] = r

  print '  dups:%12d' % dups

  print 'SAM records printed:', printed


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

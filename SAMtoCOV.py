#!/usr/bin/python

# JMG 7/21/16
# Producing methylation summary information
#   directly from a SAM file made by Bismark.

import sys
import re
import gzip

def usage():
  print '''Usage: python SAMtoCOV.py  -i <SAM>  -o <OUT>  [options]
    -i <SAM>  SAM alignment file produced by Bismark.
                Can use '-' for stdin, but must specify '-h'
                with samtools view, e.g.:
              samtools view -h <BAM> | python SAMtoCOV.py -i - -o <OUT>
    -o <OUT>  Output file listing counts of methylated and
                unmethylated CpGs, merged and sorted. Each
                line of the output lists the following six
                values for a single genomic CpG, tab-delimited:
              <chrom>          the chromosome name;
              <start>,<end>    the 1-based position of the 'C' in
                                 the CpG;
              <strand>         '+' (invariable for this merged file);
              <meth>,<unmeth>  Counts of methylated and unmethylated
                                 bases.
  Options:
    -m <int>  Minimum coverage (methylation counts) to report a CpG
                (def. 1 [all sites reported])
    -pct      Replace <strand> with methylation percent in the
                output file (fourth column)'''
  sys.exit(-1)

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

def printOutput(fOut, genome, meth, minCov, pct):
  '''
  Print the sorted output -- location and methylation counts
    for each CpG with at least minCov data values.
  '''
  printed = 0
  for chrom in genome:
    for loc in sorted(meth[chrom]):
      total = meth[chrom][loc][0] + meth[chrom][loc][1]
      if total < minCov:
        continue  # fails to meet minimum coverage
      if pct:
        # 4th column is percent methylated
        fOut.write('%s\t%d\t%d\t%.6f\t%d\t%d\n' % (chrom, loc, loc+1,
          100.0 * meth[chrom][loc][0] / total,
          meth[chrom][loc][0], meth[chrom][loc][1]))
      else:
        # 4th column is strand ('+')
        fOut.write('%s\t%d\t%d\t+\t%d\t%d\n' % (chrom, loc, loc+1,
          meth[chrom][loc][0], meth[chrom][loc][1]))
      printed += 1
  return printed

def parseCigar(cigar):
  '''
  Return in/del offset, plus string representation of cigar.
  '''
  parts = re.findall(r'(\d+)([IDM])', cigar)
  offset = 0
  cigar = ''
  for part in parts:
    if part[1] == 'D':
      offset += int(part[0])
    elif part[1] == 'I':
      offset -= int(part[0])
    cigar += int(part[0]) * part[1]
  return offset, cigar

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

def saveMeth(d, chrom, loc, meth):
  '''
  Save the methylation info for a given genomic position
    to the given dict (d).
  '''
  if chrom not in d:
    d[chrom] = {}
  if loc not in d[chrom]:
    d[chrom][loc] = [0, 0]
  if meth == 'z':
    d[chrom][loc][1] += 1  # unmethylated: index 1
  else:
    d[chrom][loc][0] += 1  # methylated: index 0

def loadMeth(cigar, strXM, chrom, pos, rc, meth, ins, dup, peMeth):
  '''
  Load methylation info using a methylation string (strXM).
  '''
  methCount = count = 0  # counting variables
  offset = 0  # in/del offset
  cigPos = 0  # position in cigar -- mirrors i (position in meth)
  for i in range(len(strXM)):

    # CpG methylation calls are 'z' (unmethylated) or 'Z' (methylated)
    while strXM[i] in ['z', 'Z']:

      # determine genomic location (C of 'CG' on the forward strand)
      loc = pos + i - rc + offset

      # for "novel" CpG sites created by a deletion,
      #   adjust location to the 5' end of the 'D's
      if rc and cigar[cigPos-1] == 'D':
        j = cigPos - 1
        while j > -1 and cigar[j] != 'M':
          j -= 1
        loc -= cigPos - j - 1

      # skip if position has been counted in a previous alignment
      if dup == 2 and chrom in peMeth and loc in peMeth[chrom]:
        break

      # for "novel" CpG sites created by an insertion,
      #   save to 'ins' dictionary
      if (not rc and cigar[cigPos] == 'I') or \
          (rc and cigPos > 0 and cigar[cigPos-1] == 'I'):
        saveMeth(ins, chrom, loc, strXM[i])
      # otherwise, save to regular 'meth' dict
      else:
        saveMeth(meth, chrom, loc, strXM[i])

      # if p-e alignment, also save to 'peMeth' dict
      if dup == 1:
        saveMeth(peMeth, chrom, loc, strXM[i])

      # update counts
      if strXM[i] == 'Z':
        methCount += 1
      count += 1
      break

    # change in/del offset
    cigPos += 1
    while cigPos < len(cigar) and cigar[cigPos] == 'D':
      offset += 1
      cigPos += 1
    if cigPos < len(cigar) and cigar[cigPos] == 'I':
      offset -= 1

  return methCount, count

def parseSAM(f, meth):
  '''
  Parse the SAM file. Save methylation data.
  '''
  genome = []  # for chromosome names, ordered in SAM header
  ins = {}     # for novel CpGs caused by insertion
  peMeth = {}  # for checking overlapping of paired-end alignments
  total = mapped = 0     # counting variables for reads
  methCount = count = 0  # counting variables for methylation data
  for line in f:

    # save chromosome names
    if line[0] == '@':
      spl = line.rstrip().split('\t')[1].split(':')
      if line[1:3] == 'SQ' and spl[0] == 'SN':
        meth[spl[1]] = {}
        genome.append(spl[1])
      continue

    # load SAM data
    spl = line.rstrip().split('\t')
    if len(spl) < 12:
      sys.stderr.write('Error! Poorly formatted SAM record\n' + line)
      sys.exit(-1)
    total += 1

    # get alignment info from flag
    flag = getInt(spl[1])
    if flag & 0x4:
      continue  # skip unmapped
    mapped += 1
    rc = 0
    if getTag(spl[11:], 'XG') == 'GA':
      rc = 1  # alignment is to G->A converted genome,
              #   so methylation data is on G of 'CG'

    # determine if read has multiple segments
    dup = 0  # 0 -> single-end alignment
    if flag & 0x1:
      if spl[0] in peMeth:
        dup = 2  # read seen before
      else:
        # first segment -- initialize dict
        peMeth[spl[0]] = {}
        dup = 1

    # load chrom, position
    chrom = spl[2]
    if chrom not in meth:
      sys.stderr.write('Error! Cannot find chromosome ' \
        + '%s in genome\n' % chrom)
      sys.exit(-1)
    pos = getInt(spl[3])

    # load CIGAR, methylation string
    offset, cigar = parseCigar(spl[5])
    strXM = getTag(spl[11:], 'XM')  # methylation string from Bismark

    # load CpG methylation info
    count1, count2 = loadMeth(cigar, strXM, chrom, pos, rc, meth, \
      ins, dup, peMeth.get(spl[0], None))
    methCount += count1
    count += count2

  # warn about novel inserted CpGs
  if ins:
    sys.stderr.write('Warning! Novel CpG(s) caused by ' \
      + 'insertion(s) -- will be ignored.\n')

  return genome, total, mapped, methCount, count

def main():
  '''
  Main.
  '''
  # Default parameters
  minCov = 1   # min. coverage to report a CpG site
  pct = 0      # report methylation percents in output
  fIn = None
  fOut = None

  # get command-line args
  args = sys.argv[1:]
  if len(args) < 4: usage()
  i = 0
  while i < len(args):
    if args[i] == '-i':
      fIn = openRead(args[i+1])
    elif args[i] == '-o':
      fOut = openWrite(args[i+1])
    elif args[i] == '-m':
      minCov = getInt(args[i+1])
    elif args[i] == '-pct':
      pct = 1
    elif args[i] == '-h':
      usage()
    else:
      sys.stderr.write('Error! Unknown parameter: %s\n' % args[i])
      usage()

    if len(args[i]) == 2:
      i += 2
    else:
      i += 1

  # check for errors
  if fIn == None or fOut == None:
    sys.stderr.write('Error! Must specify input and output files\n')
    usage()

  # process file
  meth = {}  # dict for methylation counts
  genome, total, mapped, methCount, count = parseSAM(fIn, meth)
  if fIn != sys.stdin:
    fIn.close()

  # print output
  printed = printOutput(fOut, genome, meth, minCov, pct)
  if fOut != sys.stdout:
    fOut.close()

  # print summary counts
  sys.stderr.write('Reads analyzed: %d\n' % total)
  sys.stderr.write('  Mapped: %d\n' % mapped)
  sys.stderr.write('Total CpGs analyzed: %d\n' % count)
  if count:
    sys.stderr.write('  Percent methylated: %.1f%%\n' % \
      (100.0 * methCount / count))
  sys.stderr.write('Total CpGs printed: %d\n' % printed)

if __name__ == '__main__':
  main()

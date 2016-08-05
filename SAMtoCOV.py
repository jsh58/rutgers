#!/usr/bin/python

# JMG 7/21/16
# Producing methylation summary information
#   directly from a SAM file made by Bismark.

import sys
import re
import gzip

def usage():
  print "Usage: python SAMtoCOV.py  <SAM>  <OUT>                 \n\
                                                                 \n\
    <SAM>   SAM alignment file produced by Bismark.              \n\
              Can use '-' for stdin, but must specify '-h'       \n\
              with samtools view, e.g.:                          \n\
            samtools view -h <BAM> | python SAMtoCOV.py - <OUT>  \n\
                                                                 \n\
    <OUT>   Output file listing counts of methylated and         \n\
              unmethylated CpGs (similar to that produced by     \n\
              Bismark's coverage2cytosine with --merge_CpG)       "
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
  sys.stderr.write('Error! Cannot find ' + tag + ' in SAM record\n')
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

def loadMeth(cigar, strXM, chrom, pos, rc, meth, ins):
  '''
  Load methylation info using a methylation string (strXM).
  '''
  methCount = count = 0  # counting variables
  offset = 0  # in/del offset
  cigPos = 0  # position in cigar -- mirrors i (position in meth)
  for i in range(len(strXM)):

    # CpG methylation calls are 'z' (unmethylated) or 'Z' (methylated)
    if strXM[i] in ['z', 'Z']:

      # location is C of 'CG' on the forward strand
      loc = pos + i - rc + offset

      # for "novel" CpG sites created by a deletion,
      #   adjust location to the 5' end of the 'D's
      if rc and cigar[cigPos-1] == 'D':
        j = cigPos - 1
        while j > -1 and cigar[j] != 'M':
          j -= 1
        loc -= cigPos - j - 1

      # for "novel" CpG sites created by an insertion,
      #   save to 'ins' dictionary
      if (not rc and cigar[cigPos] == 'I') or \
          (rc and cigPos > 0 and cigar[cigPos-1] == 'I'):
        saveMeth(ins, chrom, loc, strXM[i])

      # otherwise, save to regular 'meth' dict
      else:
        saveMeth(meth, chrom, loc, strXM[i])

      # update counts
      if strXM[i] == 'Z':
        methCount += 1
      count += 1

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
    #if flag & 0x10:
    if getTag(spl[11:], 'XG') == 'GA':
      rc = 1  # alignment is to G->A converted genome

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

    # check methylation string for CpG sites
    count1, count2 = loadMeth(cigar, strXM, chrom, pos, rc, meth, ins)
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
  args = sys.argv[1:]
  if len(args) < 2:
    usage()

  # open files
  f = openRead(args[0])
  fOut = openWrite(args[1])

  # process file
  meth = {}  # dict for methylation counts
  genome, total, mapped, methCount, count = parseSAM(f, meth)
  if f != sys.stdin:
    f.close()

  # print output
  for chrom in genome:
    for loc in sorted(meth[chrom]):
      #fOut.write('%s\t%d\t%d\t+\t%d\t%d\n' % (chrom, loc, loc+1, \
      #  meth[chrom][loc][0], meth[chrom][loc][1]))
      fOut.write('%s\t%d\t%d\t%.6f\t%d\t%d\n' % (chrom, loc, loc+1, \
        100.0*meth[chrom][loc][0] / (meth[chrom][loc][0]+meth[chrom][loc][1]), \
        meth[chrom][loc][0], meth[chrom][loc][1]))
  if fOut != sys.stdout:
    fOut.close()

  # print summary counts
  sys.stderr.write('Reads analyzed: %d\n' % total)
  sys.stderr.write('  Mapped: %d\n' % mapped)
  sys.stderr.write('Total CpGs analyzed: %d\n' % count)
  if count:
    sys.stderr.write('  Percent methylated: %.1f%%\n' % \
      (100.0 * methCount / count))

if __name__ == '__main__':
  main()

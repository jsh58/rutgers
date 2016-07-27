#!/usr/bin/python

# JMG 7/21/16
# Producing methylation summary information
#   directly from a SAM file.

import sys
import re

def usage():
  print "Usage: python SAMtoCOV.py  <SAM>  <OUT>                \n\
                                                                \n\
    <SAM>   SAM alignment file produced by Bismark --           \n\
              can use '-' for stdin, but must specify '-h'      \n\
              with samtools view, e.g.:                         \n\
            samtools view -h <BAM> | python SAMtoCOV.py - <OUT> \n\
                                                                \n\
    <OUT>   Output file listing counts of methylated and        \n\
              unmethylated CpGs (similar to that produced by    \n\
              Bismark's coverage2cytosine with --merge_CpG)      "
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

def parseSAM(f, d):
  '''
  Parse the SAM file. Save methylation data.
  '''
  genome = []  # ordered chromosome names
  ins = {}     # for novel CpGs caused by insertion
  count = 0
  for line in f:

    # save chromosome names
    if line[0] == '@':
      spl = line.rstrip().split('\t')[1].split(':')
      if line[1:3] == 'SQ' and spl[0] == 'SN':
        d[spl[1]] = {}
        genome.append(spl[1])
      continue

    # load SAM data
    spl = line.rstrip().split('\t')
    if len(spl) < 12:
      sys.stderr.write('Error! Poorly formatted SAM record\n' + line)
      sys.exit(-1)
    chrom = spl[2]
    pos = getInt(spl[3])
    offset, cigar = parseCigar(spl[5])
    meth = getTag(spl[11:], 'XM')  # methylation string from Bismark

    rc = 0
    if getInt(spl[1]) & 0x10:
      rc = 1  # alignment is to reverse strand
    if chrom not in d:
      sys.stderr.write('Error! Cannot find chromosome %s in genome\n' % chrom)
      sys.exit(-1)

    # check methylation string for CpG sites
    offset = 0  # reset in/del offset -- do not use net offset from parseCigar()
    cigPos = 0  # position in cigar -- mirrors i
    for i in range(len(meth)):

      if meth[i] in ['z', 'Z']:

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
          if (chrom, loc) not in ins:
            ins[(chrom, loc)] = [0, 0]
          if meth[i] == 'z':
            ins[(chrom, loc)][1] += 1
          else:
            ins[(chrom, loc)][0] += 1

        else:

          # save methylation data
          if loc not in d[chrom]:
            d[chrom][loc] = [0, 0]
          if meth[i] == 'z':
            d[chrom][loc][1] += 1
          else:
            d[chrom][loc][0] += 1

      # change in/del offset
      cigPos += 1
      while cigPos < len(cigar) and cigar[cigPos] == 'D':
        offset += 1
        cigPos += 1
      if cigPos < len(cigar) and cigar[cigPos] == 'I':
        offset -= 1

    count += 1

  # warn about novel inserted CpGs
  if ins:
    sys.stderr.write('Warning! Novel CpG(s) caused by ' + \
      'insertion(s) -- will be ignored.\n')

  return genome, count

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
  d = {}
  genome, count = parseSAM(f, d)
  if f != sys.stdin:
    f.close()

  # print output
  for chrom in genome:
    for loc in sorted(d[chrom]):
      #fOut.write('%s\t%d\t%d\t+\t%d\t%d\n' % (chrom, loc, loc+1, \
      #  d[chrom][loc][0], d[chrom][loc][1]))
      fOut.write('%s\t%d\t%d\t%.6f\t%d\t%d\n' % (chrom, loc, loc+1, \
        100.0 * d[chrom][loc][0] / ( d[chrom][loc][0] + d[chrom][loc][1] ), \
        d[chrom][loc][0], d[chrom][loc][1]))
  if fOut != sys.stdout:
    fOut.close()

  sys.stderr.write('Reads analyzed: %d\n' % count)
  #sys.stderr.write('Reads written to %s: %d\n' % (args[1], total))

if __name__ == '__main__':
  main()

#!/usr/bin/python

# JMG 4/22/16
# Combining regions for a set of bismark cov files.

import sys
import os.path
import gzip
import math

def usage():
  print '''Usage: python combineRegions5.py  [options]  -o <outfile>  <infile(s)>
    <outfile>     Output file listing combined regions (for each: genomic
                    position [chrom, start, end], number of CpGs, and
                    methylation data for each sample, tab-delimited)
    <infile(s)>   One or more files listing methylation counts at each
                    genomic position (produced by SAMtoCOV.py, or
                    coverage2cytosine from Bismark)
  Options:
    To consider a particular CpG:
      -r <int>    Minimum number of reads at a position in a sample (def. 1)
      -s <int>    Minimum number of samples with the minimum number of reads
                    to consider a position (def. 1)
    To analyze a region of CpGs:
      -d <int>    Maximum distance between CpGs to combine into the same
                    region (def. 100)
      -c <int>    Minimum number of CpGs in a region to report (def. 1)
      -x <int>    Maximum length of a region of CpGs -- regions longer than
                    this will be split into smaller regions, as long as each
                    smaller region meets the -c threshold (def. 1e9)
    To report a particular result:
      -m <int>    Minimum total reads in a region for a sample (def. 1)
    Other:
      -f          Report methylation fraction for each sample, rather than
                    methylated-unmethylated counts
      -v          Run in verbose mode'''
  sys.exit(-1)

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
  '''
  if filename == '-':
    return sys.stdout
  try:
    f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

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

def splitRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fOut):
  '''
  Split a CpG region that is too large and process
    each subregion via processRegion().
  '''
  # determine number of subregions and length
  subReg = math.ceil((reg[-1] - reg[0]) / float(maxLen))
  while len(reg) / subReg < minCpG:
    subReg -= 1  # too few CpGs: decrease number
  lengthReg = (reg[-1] - reg[0]) / subReg
  subReg = int(subReg)

  # create subregions based on length
  start = 0     # index of beginning of subregion
  prev = reg[0] # genomic position of beginning of subregion
  ends = []
  for j in range(len(reg)):
    if reg[j] > prev + lengthReg and j - start >= minCpG:
      ends.append(j)
      if len(ends) == subReg - 1: break
      start = j
      prev += lengthReg
  while len(ends) < subReg:
    ends.append(len(reg))

  # make sure each region has at least minCpG
  j = len(ends) - 1
  while j and ends[j] - ends[j - 1] < minCpG:
    ends[j - 1] = ends[j] - minCpG
    j -= 1

  # process subregions
  start = 0     # index of beginning of subregion
  total = 0     # number of regions printed
  for end in ends:
    # pass to processRegion()
    total += processRegion(chrom, reg[start:end], count, minCpG,
      minReg, float('inf'), samples, fOut)
    start = end

  return total

def processRegion(chrom, reg, count, minCpG, minReg, \
    maxLen, samples, fraction, fOut):
  '''
  Produce output for a given region of CpGs: a line
    containing chromosome name, start and end
    coordinates, number of CpGs, and methylation
    data for each sample, all tab-delimited.
  To print a line, the region must have at least
    <minCpG> sites, and at least one sample must have
    at least <minReg> counts.
  Any region longer than <maxLen> will be split via
    splitRegion(), as long as <minCpG> is still
    maintained by the subregions.
  Any sample that does not have <minReg> counts gets
    an 'NA' designation.
  '''
  if len(reg) < minCpG:
    return 0

  # split region larger than maxLen
  if reg[-1] - reg[0] > maxLen:
    return splitRegion(chrom, reg, count, minCpG, minReg, \
      maxLen, samples, fOut)

  flag = 0  # boolean for printing line
  res = '%s\t%d\t%d\t%d' % (chrom, reg[0], reg[-1], len(reg))
  for sample in samples:
    meth = unmeth = 0
    # sum methylated/unmeth bases at each position in region
    for r in reg:
      pos = str(r)
      if sample in count[chrom][pos]:
        meth += count[chrom][pos][sample][0]
        unmeth += count[chrom][pos][sample][1]
    if meth + unmeth < minReg:
      # less than minimum number of counts
      res += '\tNA'
    else:
      if fraction:
        # compute methylated fraction
        res += '\t%f' % (meth / float(meth + unmeth))
      else:
        res += '\t%d-%d' % (meth, unmeth)  # save actual counts
      flag = 1
  if flag:
    fOut.write(res + '\n')
    return 1
  return 0

def combineRegions(count, total, order, minSamples, maxDist, \
    minCpG, minReg, maxLen, samples, fraction, fOut):
  '''
  Combine data from CpG positions that are close to each
    other. Process combined regions on the fly (via
    processRegion() function).
  '''
  printed = 0  # count of printed regions
  for chrom in order:
    reg = []  # for saving connected positions
    pos3 = 0
    for pos in sorted(total[chrom], key=int):
      # require a min. number of samples
      if total[chrom][pos] >= minSamples:
        loc = int(pos)
        # if next position is more than maxDist away,
        #   process previous genomic region
        if pos3 and loc - pos3 > maxDist:
          printed += processRegion(chrom, reg, count, minCpG, \
            minReg, maxLen, samples, fraction, fOut)
          reg = []  # reset list
        reg.append(loc)
        pos3 = loc
    # process last genomic region for this chromosome
    printed += processRegion(chrom, reg, count, minCpG, \
      minReg, maxLen, samples, fraction, fOut)
  return printed

def processFile(fname, minReads, count, total, order, \
    samples):
  '''
  Load the methylation/unmethylation counts for a file.
  '''
  f = openRead(fname)

  # save sample name
  sample = fname.split('/')[-1].split('.')[0]
  while sample in samples:
    sample += '-'
  samples.append(sample)

  # load counts from file
  for line in f:
    try:
      chrom, pos, end, strand, meth, unmeth \
        = line.rstrip().split('\t')
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record' \
        + ' in %s:\n%s' % (fname, line))
      sys.exit(-1)
    meth = getInt(meth)
    unmeth = getInt(unmeth)

    # save counts and total
    if not chrom in count:
      count[chrom] = {}
      total[chrom] = {}
      order.append(chrom)
    if not pos in count[chrom]:
      count[chrom][pos] = {}
    count[chrom][pos][sample] = [meth, unmeth]
    # save to 'total' dict. only if sufficient coverage
    if meth + unmeth >= minReads:
      total[chrom][pos] = total[chrom].get(pos, 0) + 1

  f.close()

def main():
  '''
  Main.
  '''
  # Default parameters
  minReads = 1        # min. reads in a sample at a position
  minSamples = 1      # min. samples with min. reads at a position
  maxDist = 100       # max. distance between CpGs
  minCpG = 1          # min. CpGs in a region
  minReg = 1          # min. reads in a sample for a region
  maxLen = 1000000000 # max. length of a combined region
  outfile = ''        # output file
  fIn = []            # list of input files
  fraction = 0        # report methylated fractions option
  verbose = 0         # verbose option

  # Get command-line args
  args = sys.argv[1:]
  if len(args) < 2: usage()
  i = 0
  while i < len(args):
    if args[i][0] == '-':
      if args[i] == '-r':
        minReads = getInt(args[i+1])
      elif args[i] == '-s':
        minSamples = getInt(args[i+1])
      elif args[i] == '-d':
        maxDist = getInt(args[i+1])
      elif args[i] == '-c':
        minCpG = getInt(args[i+1])
      elif args[i] == '-m':
        minReg = getInt(args[i+1])
      elif args[i] == '-x':
        maxLen = getInt(args[i+1])
      elif args[i] == '-o':
        outfile = args[i+1]
      elif args[i] == '-v':
        verbose = 1
        i -= 1
      elif args[i] == '-f':
        fraction = 1
        i -= 1
      elif args[i] == '-h':
        usage()
      else:
        sys.stderr.write('Error! Unknown argument: %s\n' % args[i])
        usage()
      i += 1
    else:
      fIn.append(args[i])
    i += 1

  # check for I/O errors
  if outfile == '':
    sys.stderr.write('Error! Must supply an output file\n')
    usage()
  if len(fIn) == 0:
    sys.stderr.write('Error! Must supply one or more input files\n')
    usage()
  for fname in fIn:
    if not os.path.isfile(fname):
      sys.stderr.write('Error! Cannot open input file %s\n' % fname)
      usage()
  fOut = openWrite(outfile)

  # load methylation information for each sample
  if verbose:
    sys.stderr.write('Loading methylation information\n')
  count = {}    # for methylated, unmethylated counts
  total = {}    # for number of samples with min. coverage
  order = []    # for ordered chromosome names
  samples = []  # list of sample names
  for fname in fIn:
    if verbose:
      sys.stderr.write('  file: %s\n' % fname)
    processFile(fname, minReads, count, total, order, samples)

  # produce output
  if verbose:
    sys.stderr.write('Combining regions and producing output\n')
  fOut.write('\t'.join(['chr', 'start', 'end', 'CpG'] \
    + samples) + '\n')
  printed = combineRegions(count, total, order, minSamples, \
    maxDist, minCpG, minReg, maxLen, samples, fraction, fOut)
  if verbose:
    sys.stderr.write('Regions printed: %d\n' % printed)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()

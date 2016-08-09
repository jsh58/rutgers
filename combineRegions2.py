#!/usr/bin/python

# JMG 4/22/16
# Combining regions for a set of bismark cov files.

import sys
import gzip

def usage():
  print "Usage: python combineRegions2.py  [options]  -o <outfile>  <infile(s)> \n\
    <outfile>     Output file listing combined regions                          \n\
    <infile(s)>   One or more files produced by coverage2cytosine (Bismark)     \n\
  Options:                                                                      \n\
    To consider a particular CpG:                                               \n\
      -r <int>    Minimum number of reads at a position in a sample (def. 1)    \n\
      -s <int>    Minimum number of samples with the minimum number of reads    \n\
                    to consider a position (def. 1)                             \n\
    To analyze a region of CpGs:                                                \n\
      -d <int>    Maximum distance between CpG's to combine into the same       \n\
                    region (def. 1000)                                          \n\
      -c <int>    Minimum number of CpG's in a region to report (def. 1)        \n\
    To report a particular result:                                              \n\
      -m <int>    Minimum total reads in a region for a sample (def. 1)         \n\
    Other:                                                                      \n\
      -v          Run in verbose mode                                            "
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

def valChr(k):
  '''
  For sorting chromosomes (mouse genome).
  '''
  if k == 'chrM':
    return 0
  elif k == 'chrX':
    return 20
  elif k == 'chrY':
    return 21
  return int(k[3:])

def processRegion(chrom, reg, count, minCpG, minReg, \
    samples, fOut):
  '''
  Produce output for a given region of CpGs.
  '''
  if len(reg) < minCpG:
    return 0

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
      res += '\tNA'
    else:
      # compute methylated fraction
      res += '\t%f' % (meth / float(meth+unmeth))
      flag = 1
  if flag:
    fOut.write(res + '\n')
    return 1
  return 0

def combineRegions(count, total, minSamples, maxDist, minCpG, \
    minReg, samples, fOut):
  '''
  Combine data from CpG positions that are close to each
    other. Process combined regions on the fly.
  '''
  printed = 0  # count of printed regions
  for chrom in sorted(total, key=valChr):
    reg = []  # for saving connected positions
    pos3 = 0
    for pos in sorted(total[chrom], key=int):
      # require a min. number of samples
      if total[chrom][pos] >= minSamples:
        loc = int(pos)
        # if next position is more than maxDist away,
        #   process previous genomic region
        if pos3 and loc - pos3 > maxDist:
          printed += processRegion(chrom, reg, count, \
            minCpG, minReg, samples, fOut)
          reg = []  # reset list
        reg.append(loc)
        pos3 = loc
    # process last genomic region
    printed += processRegion(chrom, reg, count, minCpG, \
      minReg, samples, fOut)
  return printed

def processFile(fname, minReads, count, total, samples):
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
  minReads = 1    # min. reads in a sample at a position
  minSamples = 1  # min. samples with min. reads at a position
  maxDist = 1000  # max. distance between CpGs
  minCpG = 1      # min. CpGs in a region
  minReg = 1      # min. reads in a sample for a region
  fOut = None     # output file
  fIn = []        # list of input files
  verbose = 0     # verbose option

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
      elif args[i] == '-o':
        fOut = openWrite(args[i+1])
      elif args[i] == '-v':
        verbose = 1
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
  if fOut == None:
    sys.stderr.write('Error! Must supply an output file\n')
    usage()
  if len(fIn) == 0:
    sys.stderr.write('Error! Must supply one or more input files\n')
    usage()

  # load methylation information for each sample
  if verbose:
    sys.stderr.write('Loading methylation information\n')
  count = {}    # for methylated, unmethylated counts
  total = {}    # for number of samples with min. coverage
  samples = []  # list of sample names
  for fname in fIn:
    if verbose:
      sys.stderr.write('  file: %s\n' % fname)
    processFile(fname, minReads, count, total, samples)

  # produce output
  if verbose:
    sys.stderr.write('Combining regions and producing output\n')
  fOut.write('\t'.join(['chr', 'start', 'end', 'CpG'] \
    + samples) + '\n')
  printed = combineRegions(count, total, minSamples, maxDist, \
    minCpG, minReg, samples, fOut)
  if verbose:
    sys.stderr.write('Regions printed: %d\n' % printed)
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()

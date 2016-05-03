#!/usr/bin/python

# JMG 4/22/16
# Combining regions for a set of bismark cov files.

import sys

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

def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    print 'Error! Cannot convert %s to int' % arg
    usage()
  return val

def valChr(k):
  '''
  For sorting chromosomes.
  '''
  if k == 'chrM':
    return 0
  elif k == 'chrX':
    return 20
  elif k == 'chrY':
    return 21
  return int(k[3:])

def val(k):
  '''
  For sorting positions.
  '''
  return int(k)

def processRegion(chr, reg, d, minCpG, minReg, samples, fOut):
  '''
  Produce output for a given region of CpGs.
  '''
  if len(reg) < minCpG:
    return

  flag = 0  # boolean for printing line
  res = '%s\t%d\t%d\t%d' % (chr, reg[0], reg[-1], len(reg))
  for sample in samples:
    meth = unmeth = 0
    for r in reg:
      pos = str(r)
      if sample in d[chr][pos]:
        meth += d[chr][pos][sample][0]
        unmeth += d[chr][pos][sample][1]
    if meth + unmeth < minReg:
      res += '\tNA'
    else:
      res += '\t%f' % (meth / float(meth+unmeth))
      flag = 1
  if flag:
    fOut.write(res + '\n')

def combineRegions(d, tot, minSamples, maxDist, minCpG, minReg, samples, fOut):
  '''
  Combine data from CpG positions that are close to each other.
  Process combined regions on the fly.
  '''
  for chr in sorted(tot, key=valChr):
    reg = []  # for saving connected positions
    pos3 = 0
    for pos in sorted(tot[chr], key=val):
      # require a min. number of samples
      if tot[chr][pos] >= minSamples:
        loc = int(pos)
        if pos3 and loc - pos3 > maxDist:
          processRegion(chr, reg, d, minCpG, minReg, samples, fOut)
          reg = []  # reset list
        reg.append(loc)
        pos3 = loc
    processRegion(chr, reg, d, minCpG, minReg, samples, fOut)

def processFile(fname, minReads, d, tot, samples):
  '''
  Load the methylation/unmethylation counts for a file.
  '''
  try:
    f = open(fname, 'rU')
  except IOError:
    print 'Error! Cannot open', fname
    usage()

  # save sample name
  sample = fname.split('.')[0]
  samples.append(sample)

  # load counts from file
  for line in f:
    try:
      chr, pos, end, pct, meth, unmeth = line.rstrip().split('\t')
    except ValueError:
      print 'Error! Poorly formatted cov record in %s:\n' % fname, line
      sys.exit(-1)
    meth = getInt(meth)
    unmeth = getInt(unmeth)

    # save counts and total
    if not chr in d:
      d[chr] = {}
      tot[chr] = {}
    if not pos in d[chr]:
      d[chr][pos] = {}
    d[chr][pos][sample] = [meth, unmeth]
    if meth + unmeth >= minReads:
      tot[chr][pos] = tot[chr].get(pos, 0) + 1

  f.close()

def main():
  '''
  Main.
  '''
  # Default parameters
  minReads = 1     # min. reads in a sample at a position
  minSamples = 1   # min. samples with min. reads at a position
  maxDist = 1000   # max. distance between CpGs
  minCpG = 1       # min. CpGs in a region
  minReg = 1       # min. reads in a sample for a region
  fOut = None      # output file
  fIn = []         # list of input files
  verbose = 0      # verbose option

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
        fOut = open(args[i+1], 'w')
      elif args[i] == '-v':
        verbose = 1
        i -= 1
      elif args[i] == '-h':
        usage()
      else:
        print 'Error! Unknown argument:', args[i]
        usage()
      i += 1
    else:
      fIn.append(args[i])
    i += 1

  # check for I/O errors
  if fOut == None:
    print 'Error! Must supply an output file'
    usage()
  if len(fIn) == 0:
    print 'Error! Must supply one or more input files'
    usage()

  # load methylation information for each sample
  if verbose:
    print 'Loading methylation information'
  d = {}    # for methylated, unmethylated counts
  tot = {}  # for number of samples with min. coverage
  samples = []  # list of sample names
  for fname in fIn:
    if verbose:
      print '  ' + fname
    processFile(fname, minReads, d, tot, samples)

  # produce output
  if verbose:
    print 'Combining regions and producing output'
  fOut.write('\t'.join(['chr', 'start', 'end', 'CpG'] + samples) + '\n')
  combineRegions(d, tot, minSamples, maxDist, minCpG, minReg, samples, fOut)
  fOut.close()

if __name__ == '__main__':
  main()

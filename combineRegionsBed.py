#!/usr/bin/python

# JMG 4/28/16
# Analyzing bismark cov files using a BED file.

import sys

def usage():
  print "Usage: python combineRegionsBed.py  [options]  -b <bedfile>  -o <outfile>  <infile(s)> \n\
    <bedfile>     BED file listing regions of interest                          \n\
    <outfile>     Output file listing combined regions                          \n\
    <infile(s)>   One or more files produced by coverage2cytosine (Bismark)     \n\
  Options:                                                                      \n\
    -g <genome>   FASTA file of reference genome (to count CpGs)                \n\
                    (warning: the genome is loaded into memory)                  "
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

def cpgcount(seq):
  '''
  Count the number of 'CG' dinucleotides in a sequence.
  '''
  seq = seq.upper()
  count = 0
  idx = seq.find('CG')
  while idx != -1:
    count += 1
    idx = seq.find('CG', idx + 1)
  return count

def openFile(fname):
  try:
    f = open(fname, 'rU')
  except IOError:
    print 'Error! Cannot open', fname
    usage()
  return f

def loadGenome(f):
  '''
  Load a fasta genome to memory.
  '''
  if f == None:
    return {}
  gen = {}
  seq = ''
  for line in f:
    if line[0] == '>':
      if seq:
        gen[chr] = seq
        seq = ''
      chr = line[1:].rstrip()
    else:
      seq += line.rstrip()
  gen[chr] = seq
  f.close()
  return gen

def processFile(fname, d, samples):
  '''
  Load the methylation/unmethylation counts for a file.
  '''
  f = openFile(fname)

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
    pos = str(int(pos)-1)  # convert pos to 0-based

    # save counts and total
    if not chr in d:
      d[chr] = {}
    if not pos in d[chr]:
      d[chr][pos] = {}
    d[chr][pos][sample] = [meth, unmeth]

  f.close()

def processRegion(line, d, samples, gen, fOut):
  '''
  Analyze a region (defined by the BED file).
  Total methylated and unmethylated counts for each sample.
  '''
  chr, st, end = line.split('\t')
  st = getInt(st)
  end = getInt(end)
  cpg = 'NA'
  if chr in gen:
    cpg = cpgcount( gen[ chr ][ st : end ] )

  res = '%s\t%d\t%d\t%s' % (chr, st, end, str(cpg))
  for sample in samples:
    meth = unmeth = 0
    if not chr in d:
      res += '\tNA'
      continue
    for r in range(st, end):
      pos = str(r)
      if not pos in d[chr]: continue
      if sample in d[chr][pos]:
        meth += d[chr][pos][sample][0]
        unmeth += d[chr][pos][sample][1]
    if meth + unmeth == 0:
      res += '\tNA'
    else:
      res += '\t%f' % (meth / float(meth+unmeth))
  fOut.write(res + '\n')


def main():
  '''
  Main.
  '''
  # Default parameters
  fOut = None      # output file
  fIn = []         # list of input files
  fBed = None      # BED file
  fGen = None      # genome file

  # Get command-line args
  args = sys.argv[1:]
  if len(args) < 2: usage()
  i = 0
  while i < len(args):
    if args[i][0] == '-':
      if args[i] == '-b':
        fBed = openFile(args[i+1])
      elif args[i] == '-g':
        fGen = openFile(args[i+1])
      elif args[i] == '-o':
        fOut = open(args[i+1], 'w')
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
  if fBed == None:
    print 'Error! Must supply a BED file'
    usage()
  if fOut == None:
    print 'Error! Must supply an output file'
    usage()
  if len(fIn) == 0:
    print 'Error! Must supply one or more input files'
    usage()

  # load genome to memory
  gen = loadGenome(fGen)

  # load methylation information for each sample
  d = {}    # for methylated, unmethylated counts
  samples = []  # list of sample names
  for fname in fIn:
    processFile(fname, d, samples)

  # produce output
  fOut.write('\t'.join(['chr', 'start', 'end', 'CpG'] + samples) + '\n')
  for line in fBed:
    processRegion(line.rstrip(), d, samples, gen, fOut)
  fBed.close()
  fOut.close()

if __name__ == '__main__':
  main()

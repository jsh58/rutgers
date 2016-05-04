#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.

import sys, math
try:
  from scipy import stats
except ImportError:
  print 'Error! Must have the scipy module installed --'
  print '  info can be found at www.scipy.org/install.html'
  sys.exit(-1)

def usage():
  print "Usage: python diffMeth.py  -1 <group1>  -2 <group2>  \         \n\
      -i <input>  -o <output>  [options]                                \n\
    <group>     Comma-separated list of sample names (as found in       \n\
                  the header of <input>)                                \n\
    <input>     File listing genomic regions and methylation results    \n\
                  (output from combineRegions2.py)                      \n\
  Options (whether or not to report a region):                          \n\
    -d <float>  Minimum methylation difference between sample groups    \n\
                  ([0-1]; def. 0 [all results reported])                \n\
    -up         Report only regions hypermethylated in group2           \n\
    -down       Report only regions hypomethylated in group2            \n\
    -p <float>  Maximum p-value ([0-1]; def. 1 [all results reported])  \n\
    -dna        Report regions whose diff is 'NA' (occurs when one      \n\
                  group has no methylation data)                        \n\
    -pna        Report regions whose p-value is 'NA' (occurs when one   \n\
                  group does not have multiple data points)              "
  sys.exit(-1)

def openFile(fname):
  '''
  Open a file for reading.
  '''
  try:
    f = open(fname, 'rU')
  except IOError:
    print 'Error! Cannot open', fname
    usage()
  return f

def getFloat(arg):
  '''
  Convert given argument to float.
  '''
  try:
    val = float(arg)
  except ValueError:
    print 'Error! Cannot convert %s to int' % arg
    usage()
  if val > 1 or val < 0:
    print 'Error! Value %s is not in [0,1]' % arg
    usage()
  return val

def getSample(csv):
  '''
  Return a list of samples.
  '''
  arr = []
  for tok in csv.split(','):
    arr.append(tok)
  return arr

def saveIndexes(fIn, sample1, sample2):
  '''
  Find indexes for samples in input file.
  '''
  idx1 = []  # for indexes
  idx2 = []
  res1 = []  # for ordered sample names
  res2 = []
  idxExtra = []  # for gene, distance, location
  resExtra = []

  # get indexes
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  for i in range(len(spl)):
    if spl[i] in sample1:
      idx1.append(i)
      res1.append(spl[i])
    elif spl[i] in sample2:
      idx2.append(i)
      res2.append(spl[i])
    elif spl[i] in ['gene', 'distance', 'location']:
      idxExtra.append(i)
      resExtra.append(spl[i])
  if len(idx1) != len(sample1) or len(idx2) != len(sample2):
    print 'Error! Cannot find all sample names in input file'
    sys.exit(-1)

  # construct header for output file
  res = spl[:4] + res1 + res2 + resExtra
  return idx1, idx2, idxExtra, res

def calcAvg(spl, idxs):
  '''
  Calculate average for a subset of values in a list.
  '''
  avg = 0
  val = []  # for saving re-ordered values (incl. NA)
  sample = []  # for saving only numerical values
  for idx in idxs:
    val.append(spl[idx])
    if spl[idx] == 'NA':
      continue
    try:
      avg += float(spl[idx])
    except ValueError:
      print 'Error! Poorly formatted record:\n', line
    sample.append(float(spl[idx]))
  if sample:
    avg /= len(sample)
  return avg, val, sample

def processLine(line, idx1, idx2, idxExtra):
  '''
  Calculate mean difference and p-value.
  '''
  spl = line.split('\t')
  if len(spl) < max(max(idx1), max(idx2)):
    print 'Error! Poorly formatted record:\n', line

  # calculate averages, difference
  avg1, val1, sample1 = calcAvg(spl, idx1)
  avg2, val2, sample2 = calcAvg(spl, idx2)
  if not sample1 and not sample2:
    return 0, 0, ''  # no data for these groups
  diff = 'NA'
  if sample1 and sample2:
    diff = avg2 - avg1

  # calculate p-value (Welch's t-test)
  if len(sample1) < 2 or len(sample2) < 2:
    pval = 'NA'
  elif diff == 0:
    pval = 1
  elif max(sample1) - min(sample1) == 0 and \
      max(sample2) - min(sample2) == 0:
    pval = 0
  else:
    pval = stats.ttest_ind(sample1, sample2, equal_var=False)[1]

  # construct record for output file
  res = spl[:4] + val1 + val2
  for idx in idxExtra:
    res.append(spl[idx])
  return diff, pval, res

def main():
  '''
  Main.
  '''
  # Default parameters
  d = 0      # min. methylation difference
  p = 1      # max. p-value
  up = 0     # report only hypermethylated (boolean)
  down = 0   # report only hypomethylated (boolean)
  dna = 0    # report results with diff of 'NA'
  pna = 0    # report results with pval of 'NA'

  # get command-line args
  args = sys.argv[1:]
  if len(args) < 2: usage()
  fIn = None
  fOut = None
  sample1 = []
  sample2 = []
  i = 0
  while i < len(args):
    if args[i] == '-i':
      fIn = openFile(args[i+1])
    elif args[i] == '-o':
      fOut = open(args[i+1], 'w')
    elif args[i] == '-1':
      sample1 = getSample(args[i+1])
    elif args[i] == '-2':
      sample2 = getSample(args[i+1])
    elif args[i] == '-d':
      d = getFloat(args[i+1])
    elif args[i] == '-p':
      p = getFloat(args[i+1])
    elif args[i] == '-up':
      up = 1
    elif args[i] == '-down':
      down = 1
    elif args[i] == '-dna':
      dna = 1
    elif args[i] == '-pna':
      pna = 1
    elif args[i] == '-h':
      usage()
    else:
      print 'Error! Unknown argument:', args[i]
      usage()

    if args[i] in ['-up', '-down', '-dna', '-pna']:
      i += 1
    else:
      i += 2

  # check for errors
  if fIn is None or fOut is None:
    print 'Error! Must specify input and output files'
    usage()
  if not sample1 or not sample2:
    print 'Error! Must have two sets of samples'
    usage()

  # save indexes of samples from header
  idx1, idx2, idxExtra, res = saveIndexes(fIn, sample1, sample2)
  fOut.write('\t'.join(res + ['diff', 'p-value']) + '\n')

  # process file
  for line in fIn:
    diff, pval, res = processLine(line.rstrip(), idx1, idx2, idxExtra)

    # determine if result should be reported
    if res:
      if diff == 'NA':
        if not dna: continue
      elif (up and diff <= 0) or (down and diff >= 0) or \
          abs(diff) < d: continue
      if pval == 'NA':
        if not pna: continue
      elif pval > p: continue
      fOut.write('\t'.join(res + [str(diff), str(pval)]) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

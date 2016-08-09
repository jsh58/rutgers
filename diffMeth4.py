#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.
# Version 2: comparing sets of samples
# Version 4: using t-test from original diffMeth.py

import sys
scip = 0  # boolean for scipy import
try:
  from scipy import stats
except ImportError:
  sys.stderr.write('Warning! The scipy module is not installed\n'
    + '  (see www.scipy.org/install.html).\n'
    + '  All p-values will be reported as "NA".\n')
  scip = 1

def usage():
  print "Usage: python diffMeth4.py  [options]  -i <input>  -o <output> \  \n\
          <groupList1>  <groupList2>  [...]                                \n\
    <groupList>  Comma-separated list of sample names (as found in         \n\
                   the header of <input>)                                  \n\
    <input>      File listing genomic regions and methylation results      \n\
                   (output from combineRegions2.py)                        \n\
  Options (whether or not to report a region):                             \n\
    -c <int>     Minimum number of CpGs in a region (def. 1)               \n\
    -d <float>   Minimum methylation difference between sample groups      \n\
                   ([0-1]; def. 0 [all results reported])                  \n\
    -p <float>   Maximum p-value ([0-1]; def. 1 [all results reported])    \n\
    -up          Report only regions hypermethylated in later group        \n\
    -down        Report only regions hypomethylated in later group         \n\
    -dna         Report regions whose diff is 'NA' (occurs when one        \n\
                   group has no methylation data)                          \n\
    -pna         Report regions whose p-value is 'NA' (occurs when one     \n\
                   group does not have multiple data points)                "
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

def getFloat(arg, minVal = None, maxVal = None):
  '''
  Convert given argument to float. Ensure it is within a
    supplied range, if applicable.
  '''
  try:
    val = float(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to float\n' % arg)
    sys.exit(-1)
  if (minVal != None and val < minVal) or \
      (maxVal != None and val > maxVal):
    sys.stderr.write('Error! Value %f is outside of range [%f,%f]\n' \
      % (val, minVal, maxVal))
    sys.exit(-1)
  return val

def getSample(csv):
  '''
  Return a list of samples.
  '''
  arr = []
  for tok in csv.split(','):
    arr.append(tok)
  return arr

def saveIndexes(fIn, samples):
  '''
  Find indexes for samples in header of input file.
  '''
  idxs = []  # for indexes
  res = []   # for ordered sample names
  for i in range(len(samples)):
    idxs.append([])
    res.append([])

  # extra fields to keep
  extra = ['gene', 'distance', 'location']
  idxExtra = []
  resExtra = []

  # get indexes
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  for i in range(len(spl)):
    for j in range(len(samples)):
      if spl[i] in samples[j]:
        idxs[j].append(i)
        res[j].append(spl[i])
    if spl[i] in extra:
      idxExtra.append(i)
      resExtra.append(spl[i])

  # make sure all samples were found
  for i in range(len(idxs)):
    if len(idxs[i]) != len(samples[i]):
      sys.stderr.write('Error! Cannot find all sample names in input file\n')
      sys.exit(-1)

  # construct header for output file
  head = spl[:4]
  for i in range(len(res)):
    head += [','.join(res[i])]
  head += resExtra
  for i in range(len(res)-1):
    for j in range(i+1, len(res)):
      base = ','.join(res[i]) + '->' + ','.join(res[j])
      head += [base + '_diff', base + '_pval']
  return idxs, idxExtra, head

def calcDiff(avg1, avg2, meth1, meth2):
  '''
  Calculate methylation difference and p-value for
    two samples.
  '''
  if avg1 == 'NA' or avg2 == 'NA':
    return 'NA', 'NA'

  # calculate methylation difference
  diff = avg2 - avg1

  # calculate p-value (Welch's t-test)
  if len(meth1) < 2 or len(meth2) < 2 or scip:
    pval = 'NA'
  elif diff == 0:
    pval = 1
  elif max(meth1) - min(meth1) == 0 and \
      max(meth2) - min(meth2) == 0:
    pval = 0
  else:
    pval = stats.ttest_ind(meth1, meth2, equal_var=False)[1]

  return diff, pval

def calcAvg(spl, idxs):
  '''
  Calculate average for a subset of values in a list.
  '''
  avg = 0.0
  vals = []  # for saving numerical values
  for idx in idxs:
    if spl[idx] == 'NA':
      continue
    val = getFloat(spl[idx])
    avg += val
    vals.append(val)
  if vals:
    avg /= len(vals)
  return avg, vals

def isPrintable(diff, pval, minDiff, maxPval, up, down, dna, pna):
  '''
  Determine if result should be printed, based
    on diff, pval, and CL parameters.
  '''
  if diff == 'NA':
    if not dna: return 0
  elif (up and diff <= 0) or (down and diff >= 0) or \
      abs(diff) < minDiff: return 0
  if pval == 'NA':
    if not pna: return 0
  elif pval > maxPval: return 0
  return 1

def processLine(line, idxs, idxExtra, minCpG, minDiff, maxPval,
    up, down, dna, pna):
  '''
  Process a line containing methylation data for a set of samples.
  '''
  spl = line.split('\t')
  if len(spl) < max([max(idx) for idx in idxs]):
    sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
    return []
  if getInt(spl[3]) < minCpG:
    return []  # fewer than min. CpGs

  # calculate averages
  avg = []   # for average methylation fractions
  meth = []  # for methylation values
  data = 0   # valid methylation data boolean
  for idx in idxs:
    val, vals = calcAvg(spl, idx)
    if vals:
      avg.append(val)
      meth.append(vals)
      data = 1
    else:
      avg.append('NA')
      meth.append('NA')
  if not data:
    return []

  # calculate differences and p-values for each possible
  #   combination of two samples
  ans = []  # for saving diffs and p-vals
  pr = 0  # boolean for printing line
  for i in range(len(avg)-1):
    for j in range(i+1, len(avg)):
      diff, pval = calcDiff(avg[i], avg[j], meth[i], meth[j])
      ans.extend((diff, pval))

      # determine if result should be printed
      if pr == 0:
        pr = isPrintable(diff, pval, minDiff, maxPval, \
          up, down, dna, pna)

  # return if no (significant) diffs for any group comparisons
  if not pr:
    return []

  # construct record for output file
  res = spl[:4]  # chrom, pos, CpGs
  for val in avg:
    res.append(str(val))  # avg. methylation fraction
  for idx in idxExtra:
    res.append(spl[idx])  # extra fields
  for val in ans:
    res.append(str(val))  # diffs and p-values
  return res

def main():
  '''
  Main.
  '''
  # Default parameters
  minCpG = 0   # min. number of CpGs
  minDiff = 0  # min. methylation difference
  maxPval = 1  # max. p-value
  up = 0       # report only hypermethylated (boolean)
  down = 0     # report only hypomethylated (boolean)
  dna = 0      # report results with diff of 'NA'
  pna = 0      # report results with pval of 'NA'
  if scip:
    pna = 1    # if no scipy, all pvals are 'NA'

  # get command-line args
  args = sys.argv[1:]
  if len(args) < 2: usage()
  fIn = None
  fOut = None
  samples = []
  i = 0
  while i < len(args):
    if args[i] == '-i':
      fIn = openRead(args[i+1])
    elif args[i] == '-o':
      fOut = openWrite(args[i+1])
    elif args[i] == '-c':
      minCpG = getInt(args[i+1])
    elif args[i] == '-d':
      minDiff = getFloat(args[i+1], 0, 1)
    elif args[i] == '-p':
      maxPval = getFloat(args[i+1], 0, 1)
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
      samples.append(getSample(args[i]))

    if args[i][0] == '-' and len(args[i]) == 2:
      i += 2
    else:
      i += 1

  # check for errors
  if fIn == None or fOut == None:
    sys.stderr.write('Error! Must specify input and output files\n')
    usage()
  if len(samples) < 2:
    sys.stderr.write('Error! Must have at least two sets of samples\n')
    usage()

  # save indexes of samples from header
  idxs, idxExtra, res = saveIndexes(fIn, samples)
  fOut.write('\t'.join(res) + '\n')

  # process file
  count = 0
  for line in fIn:
    res = processLine(line.rstrip(), idxs, idxExtra, \
      minCpG, minDiff, maxPval, up, down, dna, pna)
    if res:
      fOut.write('\t'.join(res) + '\n')
      count += 1
  sys.stderr.write('Regions printed: %d\n' % count)

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

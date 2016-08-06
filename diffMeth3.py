#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.
# Version 2: comparing sets of samples

import sys, math
scip = 0  # boolean for scipy import
try:
  from scipy import stats
except ImportError:
  sys.stderr.write('Warning! The scipy module is not installed\n'
    + '  (see www.scipy.org/install.html).\n'
    + '  All p-values will be reported as "NA".\n')
  scip = 1

def usage():
  print "Usage: python diffMeth3.py  [options]  -i <input>  -o <output> \  \n\
          <groupList1>  <groupList2>  [...]                                \n\
    <groupList>  Comma-separated list of sample names (as found in         \n\
                   the header of <input>)                                  \n\
    <input>      File listing genomic regions and methylation results      \n\
                   (output from combineRegions2.py)                        \n\
  Options (whether or not to report a region):                             \n\
    -c <int>     Minimum number of CpGs in a region (def. 1)               \n\
    -d <float>   Minimum methylation difference between sample groups      \n\
                   ([0-1]; def. 0 [all results reported])                  \n\
    -p <float>   Maximum p-value ([0-1]; def. 1 [all results reported])     "
  sys.exit(-1)

def openFile(fname):
  '''
  Open a file for reading.
  '''
  try:
    f = open(fname, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % fname)
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
  Find indexes for samples in input file.
  '''
  idxs = []  # for indexes
  res = []   # for ordered sample names
  for i in range(len(samples)):
    idxs.append([])
    res.append([])

  idxExtra = []  # for gene, distance, location
  resExtra = []

  # get indexes
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  for i in range(len(spl)):
    for j in range(len(samples)):
      if spl[i] in samples[j]:
        idxs[j].append(i)
        res[j].append(spl[i])
    if spl[i] in ['gene', 'distance', 'location']:
      idxExtra.append(i)
      resExtra.append(spl[i])

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
      head += ['->'.join([str(i), str(j)]) + '_diff', \
        '->'.join([str(i), str(j)]) + '_pval']
  return idxs, idxExtra, head

def countBases(spl, idxs):
  '''
  Sum methylated and unmethylated bases from a list.
  '''
  meth = unmeth = 0
  sample = 0
  for idx in idxs:
    if spl[idx] == 'NA':
      continue
    div = spl[idx].split('-')
    meth += getInt(div[0])
    unmeth += getInt(div[1])
  return meth, unmeth

def calcDiff(count1, count2):
  '''
  Calculate methylation difference and p-value for
    two samples.
  '''
  if count1[0] == 'NA' or count2[0] == 'NA':
    return 'NA', 'NA'

  # calculate methylation difference
  n1 = count1[0] + count1[1]
  n2 = count2[0] + count2[1]
  p1 = float(count1[0]) / n1
  p2 = float(count2[0]) / n2
  diff = p2 - p1

  # calculate Z-value
  phat = float(count1[0] + count2[0]) / (n1 + n2)
  if phat in [0.0, 1.0]:
    zval = 0
  else:
    zval = diff / (phat * (1-phat) * (1.0 / n1 + 1.0 / n2))**0.5

  # calculate two-tailed p-value
  pval = stats.norm.sf(abs(zval))*2
  return diff, pval

def processLine(line, idxs, idxExtra, minCpG, minDiff, maxPval):
  '''
  Process a line containing methylation data for a set of samples.
  '''
  spl = line.split('\t')
  if len(spl) < max([max(idx) for idx in idxs]):
    sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
  if int(spl[3]) < minCpG:
    return []  # fewer than min. CpGs

  # save counts of meth/unmeth bases
  count = []
  for idx in idxs:
    meth, unmeth = countBases(spl, idx)
    if meth or unmeth:
      count.append([meth, unmeth])
    else:
      count.append(['NA'])

  # calculate differences and p-values for each possible
  #   combination of two samples
  ans = []  # for saving diffs and p-vals
  pr = 0  # boolean for printing line
  for i in range(len(count)-1):
    for j in range(i+1, len(count)):
      diff, pval = calcDiff(count[i], count[j])
      ans.extend((diff, pval))
      if diff != 'NA' and abs(diff) >= minDiff \
          and pval <= maxPval:
        pr = 1

  # return if no (significant) diffs for any groups
  if not pr:
    return []

  # construct record for output file
  res = spl[:4]
  for c in count:
    if c[0] == 'NA':
      res.append('NA')
    else:
      frac = float(c[0]) / (c[0] + c[1])
      res.append(str(frac)) # methylation fraction
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

  # get command-line args
  args = sys.argv[1:]
  if len(args) < 2: usage()
  fIn = None
  fOut = None
  samples = []
  i = 0
  while i < len(args):
    if args[i] == '-i':
      fIn = openFile(args[i+1])
    elif args[i] == '-o':
      fOut = open(args[i+1], 'w')
    elif args[i] == '-c':
      minCpG = getInt(args[i+1])
    elif args[i] == '-d':
      minDiff = getFloat(args[i+1], 0, 1)
    elif args[i] == '-p':
      maxPval = getFloat(args[i+1], 0, 1)
    elif args[i] == '-h':
      usage()
    else:
      sample = getSample(args[i])
      samples.append(sample)

    if args[i][0] != '-':
      i += 1
    else:
      i += 2

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
  for line in fIn:
    res = processLine(line.rstrip(), idxs, idxExtra, \
      minCpG, minDiff, maxPval)
    if res:
      fOut.write('\t'.join(res) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

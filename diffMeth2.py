#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.
# Version 2: comparing sets of samples

import sys, math
#scip = 0  # boolean for scipy import
#try:
#  from scipy import stats
#except ImportError:
#  print 'Warning! The scipy module is not installed'
#  print '  (see www.scipy.org/install.html).'
#  print '  All p-values will be reported as "NA".'
#  scip = 1

def usage():
  print "Usage: python diffMeth2.py  [options]  -i <input>  -o <output> \ \n\
      <groupList>                                                       \n\
    <groupList> Comma-separated list of sample names (as found in       \n\
                  the header of <input>)                                \n\
    <input>     File listing genomic regions and methylation results    \n\
                  (output from combineRegions2.py)                      \n\
  Options (whether or not to report a region):                          \n\
    -c <int>    Minimum number of CpGs in a region (def. 1)             \n\
    -d <float>  Minimum methylation difference between sample groups    \n\
                  ([0-1]; def. 0 [all results reported])                 "
    #-p <float>  Maximum p-value ([0-1]; def. 1 [all results reported])   "
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

def getFloat(arg):
  '''
  Convert given argument to float.
  '''
  try:
    val = float(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to float\n' % arg)
    sys.exit(-1)
  if val > 1 or val < 0:
    sys.stderr.write('Error! Value %s is not in [0,1]\n' % arg)
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
      head += ['->'.join([str(i), str(j)])]
  return idxs, idxExtra, head

def calcAvg(spl, idxs):
  '''
  Calculate average for a subset of values in a list.
  '''
  avg = 0
  sample = 0
  for idx in idxs:
    if spl[idx] == 'NA':
      continue
    try:
      avg += float(spl[idx])
    except ValueError:
      sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
    sample += 1
  if sample:
    avg /= sample
  return avg, sample

def processLine(line, cpg, idxs, idxExtra):
  '''
  Calculate mean difference and p-value.
  '''
  spl = line.split('\t')
  if len(spl) < max([max(idx) for idx in idxs]):
    sys.stderr.write('Error! Poorly formatted record:\n%s' % line)
  if int(spl[3]) < cpg:
    return [], 0  # fewer than min. CpGs

  # calculate averages
  avg = []
  for idx in idxs:
    av, sample = calcAvg(spl, idx)
    if sample:
      avg.append(av)
    else:
      avg.append('NA')

  # calculate differences
  diff = []
  maxDiff = -1
  pr = 0  # boolean for printing line
  for i in range(len(avg)-1):
    for j in range(i+1, len(avg)):
      if avg[i] == 'NA' or avg[j] == 'NA':
        diff.append('NA')
      else:
        dif = avg[j] - avg[i]
        diff.append(dif)
        if abs(dif) > maxDiff:
          maxDiff = abs(dif)
        pr = 1
  if not pr:
    return [], 0  # no diffs for these groups

  # calculate p-value (Welch's t-test)
  #if len(sample1) < 2 or len(sample2) < 2 or scip:
  #  pval = 'NA'
  #elif diff == 0:
  #  pval = 1
  #elif max(sample1) - min(sample1) == 0 and \
  #    max(sample2) - min(sample2) == 0:
  #  pval = 0
  #else:
  #  pval = stats.ttest_ind(sample1, sample2, equal_var=False)[1]

  # construct record for output file
  res = spl[:4]
  for av in avg:
    res.append(str(av))
  for idx in idxExtra:
    res.append(spl[idx])
  for dif in diff:
    res.append(str(dif))
  return res, maxDiff

def main():
  '''
  Main.
  '''
  # Default parameters
  cpg = 0    # min. number of CpGs
  d = 0      # min. methylation difference
  #p = 1      # max. p-value

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
      cpg = getInt(args[i+1])
    elif args[i] == '-d':
      d = getFloat(args[i+1])
    #elif args[i] == '-p':
      #p = getFloat(args[i+1])
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
    res, maxDiff = processLine(line.rstrip(), cpg, \
      idxs, idxExtra)

    # print result
    if res and maxDiff >= d:
      fOut.write('\t'.join(res) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

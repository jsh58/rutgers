#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.

import sys, math
try:
  from scipy import stats
except ImportError:
  print 'Error! Must have the scipy module installed'
  sys.exit(-1)

def usage():
  print "Usage: python diffMeth.py  -1 <sample1List>  -2 <sample2List>  \ \n\
    -i <input>  -o <output>                                               \n\
  <sampleList>  Comma-separated list of samples                           \n\
  <input>       Output from combineRegions2.py                             "
  sys.exit(-1)

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
    return 0, 0, ''
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
      try:
        fIn = open(args[i+1], 'rU')
      except IOError:
        print 'Error! Cannot open', args[i+1]
        usage()
    elif args[i] == '-o':
      fOut = open(args[i+1], 'w')
    elif args[i] == '-1':
      sample1 = getSample(args[i+1])
    elif args[i] == '-2':
      sample2 = getSample(args[i+1])
    elif args[i] == '-h':
      usage()
    else:
      print 'Error! Unknown argument:', args[i]
      usage()
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
    if res:
      fOut.write('\t'.join(res + [str(diff), str(pval)]) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

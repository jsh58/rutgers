#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.

import sys, math
from scipy import stats

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

def processLine(line, idx1, idx2):
  '''
  Calculate mean difference and p-value.
  '''
  spl = line.split('\t')
  if len(spl) < max(max(idx1), max(idx2)):
    print 'Error! Poorly formatted record:\n', line

  # calculate avg for sample1
  avg1 = 0
  sample1 = []
  for idx in idx1:
    if spl[idx] == 'NA': continue
    try:
      avg1 += float(spl[idx])
    except ValueError:
      print 'Error! Poorly formatted record:\n', line
    sample1.append(float(spl[idx]))
  if sample1:
    avg1 /= len(sample1)

  # calculate avg for sample2
  avg2 = 0
  sample2 = []
  for idx in idx2:
    if spl[idx] == 'NA': continue
    try:
      avg2 += float(spl[idx])
    except ValueError:
      print 'Error! Poorly formatted record:\n', line
    sample2.append(float(spl[idx]))
  if sample2:
    avg2 /= len(sample2)
  diff = abs(avg1 - avg2)

  # calculate p-value (Welch's t-test)
  if len(sample1) < 2 or len(sample2) < 2 or \
      (max(sample1) - min(sample1) == 0 and \
      max(sample2) - min(sample2) == 0):
    pval = 'NA'
  elif diff == 0:
    pval = 1
  else:
    pval = stats.ttest_ind(sample1, sample2, equal_var=False)[1]
  return diff, pval

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
  idx1 = []
  idx2 = []
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  for i in range(len(spl)):
    if spl[i] in sample1:
      idx1.append(i)
    elif spl[i] in sample2:
      idx2.append(i)
  if len(idx1) != len(sample1) or len(idx2) != len(sample2):
    print 'Error! Cannot find all sample names in input file'
    sys.exit(-1)
  fOut.write('\t'.join([header, 'diff', 'p-value']) + '\n')

  # process file
  for line in fIn:
    diff, pval = processLine(line.rstrip(), idx1, idx2)
    fOut.write('\t'.join([line.rstrip(), str(diff), str(pval)]) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

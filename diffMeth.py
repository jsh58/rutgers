#!/usr/bin/python

# JMG 4/26/16
# Testing regions for differential methylation.

import sys, math
from scipy import stats

def usage():
  print "Usage: python diffMeth.py  <input>  <output>"
  print '  <input>   Output from combineRegions2.py'
  sys.exit(-1)

def processLine(line, idx1, idx2):
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
  args = sys.argv[1:]
  if len(args) < 2: usage()
  try:
    fIn = open(args[0], 'rU')
  except IOError:
    print 'Error! Cannot open', args[0]
    usage()
  fOut = open(args[1], 'w')

  # save indexes of samples from header
  sam1 = ['C14', 'C20']
  sam2 = ['C34', 'C40']
  idx1 = []
  idx2 = []
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  for i in range(len(spl)):
    if spl[i] in sam1:
      idx1.append(i)
    elif spl[i] in sam2:
      idx2.append(i)

  # process file
  for line in fIn:
    diff, pval = processLine(line.rstrip(), idx1, idx2)
    fOut.write('\t'.join([line.rstrip(), str(diff), str(pval)]) + '\n')

  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

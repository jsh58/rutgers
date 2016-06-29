#!/usr/bin/python

# JMG 6/24/16
# Retrieving a subset of reads from a SAM.
# Counting good/bad qual methylation data

import sys

def getXM(lis):
  '''
  Get "XM" methylation string.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == 'XM':
      return spl[-1]
  sys.stderr.write('Error! Cannot find XM in SAM record\n')
  sys.exit(-1)

def parseSAM(f, lis, minQual):
  '''
  Parse the SAM file. Select reads overlapping the
    given coordinates.
  '''
  res = [[0] * len(lis) for i in range(2)]
  count = 0
  for line in f:
    if line[0] == '@':
      continue
    spl = line.rstrip().split('\t')
    count += 1

    qual = spl[10]
    xm = getXM(spl[11:])
    for i in range(len(xm)):
      if xm[i] in lis:
        idx = lis.index(xm[i])
        if ord(qual[i]) - 33 < minQual:
          res[0][idx] += 1
        else:
          res[1][idx] += 1

  return count, res

def main():
  '''
  Main.
  '''
  args = sys.argv[1:]
  if len(args) < 1:
    print 'Usage: python %s  <SAMfile>  [<qual>]' % sys.argv[0]
    print '  Use \'-\' for stdin/stdout'
    print '  <qual>   Quality score for good/bad classification (def. 30)'
    sys.exit(-1)
  try:
    f = open(args[0], 'rU')
  except IOError:
    # use '-' for stdin
    if args[0] == '-':
      f = sys.stdin
    else:
      sys.stderr.write('Error! Cannot open %s\n' % args[0])
      sys.exit(-1)

  # get minimum quality value
  minQual = 30
  if len(args) > 1:
    try:
      minQual = int(args[1])
    except ValueError:
      sys.stderr.write('Error! Cannot convert %s to int\n' % args[1])

  # process file
  lis = ['h', 'H', 'x', 'X', 'z', 'Z']
  count, res = parseSAM(f, lis, minQual)
  f.close()

  # print results
  for c in lis:
    print '\t' + c,
  print
  for i in range(2):
    if i: print 'Good\t',
    else: print 'Bad\t',
    for j in range(len(lis)):
      print '%d\t' % res[i][j],
    print
  print 'BadFrac',
  for j in range(len(lis)):
    print '\t%.3f' % (res[0][j] / float(res[0][j] + res[1][j])),
  print

  sys.stderr.write('Reads analyzed: %d\n' % count)

if __name__ == '__main__':
  main()

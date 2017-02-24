#!/usr/bin/python

# JMG 9/21/16
# Convert featureCounts output for use with
#   edgeR or DESeq2.
# Version 4: keep transcript lengths

import sys

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

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <in>  <out>\n' % sys.argv[0])
    sys.exit(-1)
  f1 = openRead(args[0])
  f2 = openWrite(args[1])

  line = f1.readline()
  line = f1.readline()
  spl = line.rstrip().split('\t')
  f2.write('Symbol\tLength')
  for sample in spl[6:]:
    div = sample.split('/')
    f2.write('\t%s' % (div[-2]))
  f2.write('\n')
  for line in f1:
    spl = line.rstrip().split('\t')
    f2.write('%s' % (spl[0]))
    for val in spl[5:]:
      f2.write('\t%d' % int(round(float(val))))
    f2.write('\n')
  f1.close()
  f2.close()

if __name__ == '__main__':
  main()

#!/usr/bin/python

# JMG 10/6/16
# Exclude reads listed in a given headers.txt
#   file (e.g. produced by findDups.py) from
#   a SAM file.

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
  if len(args) < 4:
    sys.stderr.write('Usage: python getReads.py  <inSAM>  ' \
      + '(yes|no)  <headers.txt>  <outSAM>\n' \
      + '  Use \'-\' for stdin/stdout\n' \
      + '  (yes|no)  Keep reads in headers.txt? \'yes\' or \'no\'\n')
    sys.exit(-1)

  # load headers -- 1st column of file
  head = {}
  fHead = openRead(args[2])
  for line in fHead:
    head[line.rstrip().split('\t')[0]] = 1
  if fHead != sys.stdin:
    fHead.close()

  # open SAM files
  fIn = openRead(args[0])
  fOut = openWrite(args[3])
  yes = 0
  if args[1] == 'yes':
    yes = 1

  # process SAM file
  if yes:
    for line in fIn:
      if line[0] == '@' or line.split('\t')[0] in head:
        fOut.write(line)
  else:
    for line in fIn:
      if line[0] == '@' or line.split('\t')[0] not in head:
        fOut.write(line)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()

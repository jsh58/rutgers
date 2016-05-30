#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)

import sys
import re

def parseCigar(cigar):
  '''
  Determine indel offset of alignment from CIGAR.
  '''
  offset = 0
  ins = re.findall(r'(\d+)I', cigar)
  for i in ins:
    offset -= int(i)
  de = re.findall(r'(\d+)D', cigar)
  for d in de:
    offset += int(d)
  return offset

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    print 'Usage: python %s  <SAMfile>  <out>' % sys.argv[0]
    sys.exit(-1)
  try:
    f = open(args[0], 'rU')
  except IOError:
    # use '-' for stdin
    if args[0] == '-':
      f = sys.stdin
    else:
      print 'Error! Cannot open', args[0]
      sys.exit(-1)

  # use '-' for stdout
  if args[1] == '-':
    fOut = sys.stdout
  else:
    fOut = open(args[1], 'w')

  pos = {}  # dict of positions
  count = dups = uniq = 0
  for line in f:
    if line[0] == '@':
      fOut.write(line)
      continue
    spl = line.rstrip().split('\t')

    # determine location -- use 3' end if RC
    chr = spl[2]
    loc = spl[3]
    rc = 0
    if int(spl[1]) & 0x10:
      loc = str(int(spl[3]) + len(spl[9]) + parseCigar(spl[5]))
      rc = 1

    if chr in pos:
      if loc in pos[chr]:
        if rc in pos[chr][loc]:
          dups += 1
        else:
          pos[chr][loc][rc] = 1
          fOut.write(line)
          uniq += 1
      else:
        pos[chr][loc] = {}
        pos[chr][loc][rc] = 1
        fOut.write(line)
        uniq += 1
    else:
      pos[chr] = {}
      pos[chr][loc] = {}
      pos[chr][loc][rc] = 1
      fOut.write(line)
      uniq += 1
    count += 1
  f.close()
  fOut.close()

  #print 'Reads: ', count
  #print 'Dups:   ', dups
  #print 'Unique:', uniq

if __name__ == '__main__':
  main()

#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches

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
  count = dups = uniq = notSeq = 0
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
          if spl[9] in pos[chr][loc][rc]:
            # location, strand, *and* seq match
            dups += 1
          else:
            # everything matches except seq
            pos[chr][loc][rc].append(spl[9])
            fOut.write(line)
            notSeq += 1
        else:
          pos[chr][loc][rc] = [spl[9]]
          fOut.write(line)
          uniq += 1
      else:
        pos[chr][loc] = {}
        pos[chr][loc][rc] = [spl[9]]
        fOut.write(line)
        uniq += 1
    else:
      pos[chr] = {}
      pos[chr][loc] = {}
      pos[chr][loc][rc] = [spl[9]]
      fOut.write(line)
      uniq += 1
    count += 1
  f.close()
  fOut.close()

  sys.stderr.write('Reads: %10d\n' % count)
  sys.stderr.write('Unique: %9d\n' % uniq)
  sys.stderr.write('Dups: %11d\n' % dups)
  sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

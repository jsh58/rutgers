#!/usr/bin/python

# JMG 5/29/16
# Finding putative PCR duplicates.
# (mimicking bismark's deduplicate script)
# Version 3: counting results, including actual seq matches
# Version 5: streamlining dict, using consensus seq

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
  if len(args) < 1:
    print 'Usage: python %s  <SAMfile>  [<out>]' % sys.argv[0]
    print '  Use \'-\' for stdin/stdout'
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

  out = 0
  if len(args) > 1:
    # use '-' for stdout
    if args[1] == '-':
      fOut = sys.stdout
    else:
      fOut = open(args[1], 'w')
    out = 1

  pos = {}  # dict of positions
  first = {}
  count = dups = uniq = notSeq = 0
  for line in f:
    if line[0] == '@':
      if out: fOut.write(line)
      continue
    spl = line.rstrip().split('\t')

    # determine location -- use 3' end if RC
    chr = spl[2]
    loc = spl[3]
    rc = 0
    if int(spl[1]) & 0x10:
      loc = str(int(spl[3]) + len(spl[9]) + parseCigar(spl[5]))
      rc = 1

    if (chr, loc, rc) in pos:
      flag = 1
      for seq in pos[(chr, loc, rc)]:
        if seq[0] == spl[9]:
          seq[1] += 1
          dups += 1
          flag = 0
          break
      if flag:
        pos[(chr, loc, rc)].append([spl[9], 1])
        first[(chr, loc, rc, spl[9])] = line
        uniq += 1
      #pos[(chr, loc, rc)][spl[9]] = pos[(chr, loc, rc)].get(spl[9], 0) + 1

    else:
      pos[(chr, loc, rc)] = [[spl[9], 1]]
      #pos[(chr, loc, rc)][spl[9]] = 1
      first[(chr, loc, rc, spl[9])] = line
      uniq += 1

    count += 1
  f.close()

  # produce output
  if out:
    for p in pos:

      if len(pos[p]) > 1:
        fOut.write(first[p + (max(pos[p], key=lambda x: x[1])[0],)])
      else:
        fOut.write(first[p + (pos[p][0][0],)])

    fOut.close()


  sys.stderr.write('Reads: %10d\n' % count)
  sys.stderr.write('Unique: %9d\n' % uniq)
  sys.stderr.write('Dups: %11d\n' % dups)
  #sys.stderr.write('NotSeq: %9d\n' % notSeq)

if __name__ == '__main__':
  main()

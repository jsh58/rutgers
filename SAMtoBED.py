#!/usr/bin/python

# JMG 4/13/16

# Produce a BED file from a SAM (or piped-in BAM).
# Combine paired-end alignments.

import sys
import re

def usage():
  print "Usage: python SAMtoBED.py  [options]  <input>  <output>  \n\
    <input>     SAM file, or '-' for stdin                        \n\
    <output>    BED file                                          \n\
  Options for unpaired alignments:                                \n\
    -n        Do not print unpaired alignments (default)          \n\
    -y        Print unpaired alignments                           \n\
    -a <int>  Print unpaired alignments, with read length         \n\
                increased to specified value"
  sys.exit(-1)

def loadChrLen(line, chr):
  '''
  Load chromosome lengths from the SAM file header.
  '''
  match = re.search(r'@SQ\s+SN:(\S+)\s+LN:(\d+)', line)
  if match:
    chr[match.group(1)] = int(match.group(2))

def writeOut(fOut, ref, start, end, read, chr):
  '''
  Write BED output. Check chromosome lengths.
  '''
  fOut.write('%s\t%d\t%d\t%s\n' % \
    (ref, start, end, read))

def processPaired(spl, flag, start, pos, chr, fOut):
  '''
  Process a properly paired SAM record.
  '''
  # 2nd of PE reads
  if spl[0] in pos:
    if pos[spl[0]] == -1:
      print 'Error! Read %s already analyzed' % spl[0]
      sys.exit(-1)

    # save end position
    if flag & 0x10:
      start += len(spl[9])

    writeOut(fOut, spl[2], min(start, pos[spl[0]]), \
      max(start, pos[spl[0]]), spl[0], chr)

    pos[spl[0]] = -1

  # 1st of PE reads: save end position
  else:
    if flag & 0x10:
      pos[spl[0]] = start + len(spl[9])
    else:
      pos[spl[0]] = start


def parseSAM(fIn, fOut, add):
  '''
  Parse the input file, and produce the output file.
  '''
  chr = {}
  pos = {}
  line = fIn.readline().rstrip()
  while line:

    # skip header
    if line[0] == '@':
      loadChrLen(line, chr)  # load chromosome length
      line = fIn.readline().rstrip()
      continue

    # save flag and start position
    spl = line.split('\t')
    if len(spl) < 11:
      print 'Error! Poorly formatted SAM record:\n', line
      sys.exit(-1)
    try:
      flag = int(spl[1])
      start = int(spl[3]) - 1
    except ValueError:
      print 'Error parsing SAM record:\n', line
      sys.exit(-1)

    # skip unmapped
    if flag & 0x4:
      line = fIn.readline().rstrip()
      continue

    # properly paired alignment
    if flag & 0x2:
      processPaired(spl, flag, start, pos, chr, fOut)

    # unpaired alignment
    elif add != -1:
      end = start + len(spl[9])
      # adjust ends (parameter 'add')
      if add != 0:
        if flag & 0x10:
          start -= add - len(spl[9])
        else:
          end = start + add

      out.write('%s\t%d\t%d\t%s\n' % \
        (spl[2], start, end, spl[0]))

      pos[spl[0]] = -1

    line = fIn.readline().rstrip()

def main():
  '''
  Runs the program.
  '''
  # parse command-line args
  args = sys.argv[1:]
  if not args: usage()
  add = -1  # number of bp to add to unpaired reads
  i = 0
  if args[i] == '-y':
    add = 0
    i += 1
  elif args[i] == '-a':
    try:
      add = int(args[i+1])
    except ValueError:
      print 'Error! Cannot convert %s to int' % args[i+1]
      sys.exit(-1)
    i += 2
  elif args[i] == '-n':
    i += 1

  # open files
  if len(args[i:]) < 2: usage()
  try:
    fIn = open(args[i], 'rU')
  except IOError:
    if args[i] == '-':
      fIn = sys.stdin
    else:
      print 'Error! Cannot open', args[i]
      usage()
  fOut = open(args[i+1], 'w')

  parseSAM(fIn, fOut, add)
  fIn.close()
  fOut.close()

if __name__ == '__main__':
  main()

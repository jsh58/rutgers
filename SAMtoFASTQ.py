#!/usr/bin/python

# JMG 4/6/16

# Produce a FASTQ file from a SAM (or piped-in BAM).
# Avoids duplicated reads.

import sys

def revComp(dna):
  '''
  Reverse-complements the given DNA sequence.
  '''
  rc = ''
  for nuc in dna[::-1]:
    comp = ''
    if nuc == 'A': comp = 'T'
    elif nuc == 'C': comp = 'G'
    elif nuc == 'G': comp = 'C'
    elif nuc == 'T': comp = 'A'
    elif nuc == 'N': comp = 'N'
    else:
      print 'Error! Unknown nucleotide: %s' % nuc
    rc += comp
  return rc

def parseSAM(f, out):
  '''
  Parse the input file, and produce the output file.
  Avoid duplicating reads. Rev-comp reads if necessary.
  '''
  headers = {}
  line = f.readline().rstrip()
  while line:

    # skip header
    if line[0] == '@':
      line = f.readline().rstrip()
      continue

    spl = line.split('\t')
    if len(spl) < 11:
      print 'Error! Poorly formatted SAM record:\n', line
      sys.exit(-1)

    # skip duplicates
    if spl[0] in headers:
      line = f.readline().rstrip()
      continue

    # get FLAG
    try:
      spl[1] = int(spl[1])
    except ValueError:
      print 'Error parsing FLAG in SAM record:\n', line
      sys.exit(-1)

    # skip if secondary or supplementary
    if spl[1] & 0x900:
      line = f.readline().rstrip()
      continue

    # rev-comp if necessary
    if spl[1] & 0x10:
      spl[9] = revComp(spl[9])
      spl[10] = spl[10][::-1]
    out.write('@%s\n%s\n+\n%s\n' % (spl[0], spl[9], spl[10]))
    headers[spl[0]] = 1

    line = f.readline().rstrip()

def main():
  '''
  Runs the program.
  '''
  args = sys.argv[1:]
  if len(args) < 2:
    print 'Usage: python SAMtoFASTQ.py  <input>  <output>'
    print "    <input> should be a SAM file, or '-' for stdin"
    sys.exit(-1)
  try:
    f = open(args[0], 'rU')
  except IOError:
    if args[0] == '-':
      f = sys.stdin
    else:
      print 'Error! Cannot open', args[0]
      sys.exit(-1)
  out = open(args[1], 'w')

  parseSAM(f, out)
  f.close()
  out.close()

if __name__ == '__main__':
  main()

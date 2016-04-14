#!/usr/bin/python

# JMG 4/6/16

# Produce a FASTQ file from a SAM (or piped-in BAM).
# Avoids duplicated reads.
# Version 2: dealing with paired reads.

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

def parseSAM(f, out1, out2, out3):
  '''
  Parse the input file, and produce the output files.
  Avoid duplicating reads. Rev-comp reads if necessary.
  '''
  headers = {}
  paired = {}
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

    # if paired, check other record or save info
    if spl[1] & 0x2:
      first = ''
      if spl[0] in paired and \
        paired[spl[0]][1] == spl[7] and paired[spl[0]][2] == spl[3]:
        # making sure coordinates match
        first = paired[spl[0]][0]

      # print both records
      if first:
        # rev-comp if necessary
        if spl[1] & 0x10:
          spl[9] = revComp(spl[9])
          spl[10] = spl[10][::-1]
        record = '@%s\n%s\n+\n%s\n' % (spl[0] + ' /2', spl[9], spl[10])
        out1.write(first)
        out2.write(record)
        headers[spl[0]] = 1
        del paired[spl[0]]

      # other record not loaded yet: save info
      else:
        paired[spl[0]] = []
        # rev-comp if necessary
        if spl[1] & 0x10:
          spl[9] = revComp(spl[9])
          spl[10] = spl[10][::-1]
        record = '@%s\n%s\n+\n%s\n' % (spl[0] + ' /1', spl[9], spl[10])
        paired[spl[0]].append(record)
        paired[spl[0]].append(spl[3])
        paired[spl[0]].append(spl[7])

    # unpaired alignment: print it
    else:
      # rev-comp if necessary
      if spl[1] & 0x10:
        spl[9] = revComp(spl[9])
        spl[10] = spl[10][::-1]

      out3.write('@%s\n%s\n+\n%s\n' % (spl[0], spl[9], spl[10]))
      headers[spl[0]] = 1

    line = f.readline().rstrip()

  # print paired leftovers
  if paired:
    print 'Error! Reads labeled paired but without mate:'
    for read in paired:
      print read, paired[read]

def main():
  '''
  Runs the program.
  '''
  args = sys.argv[1:]
  if len(args) < 4:
    print 'Usage: python SAMtoFASTQ.py  <input>  <outputPE1>  <outputPE2>  <unpaired>'
    print "    <input> should be a SAM file, or '-' for stdin"
    sys.exit(-1)
  try:
    f = open(args[0], 'rU')
  except IOError:
    if args[0] is '-':
      f = sys.stdin
    else:
      print 'Error! Cannot open', args[0]
      sys.exit(-1)
  out1 = open(args[1], 'w')
  out2 = open(args[2], 'w')
  out3 = open(args[3], 'w')

  parseSAM(f, out1, out2, out3)
  f.close()
  out1.close()
  out2.close()
  out3.close()

if __name__ == '__main__':
  main()

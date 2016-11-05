#!/usr/bin/python

# JMG 11/5/16

# Simulate a ChIP experiment.

import sys
import random

def revComp(seq):
  rc = ''
  comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
  for nuc in seq[::-1]:
    if not nuc.upper() in comp:
      print 'Error! Unknown nucleotide:', nuc
      sys.exit(-1)
    rc += comp[nuc.upper()]
  return rc

def loadGenome(f1):
  '''
  Load fasta genome to dict.
  '''
  sys.stderr.write('Loading genome...\n')
  gen = {}
  length = {}
  chrom = ''
  seq = ''
  for line in f1:
    if line[0] == '>':
      if seq:
        gen[chrom] = seq
        length[chrom] = len(seq)
        seq = ''
      chrom = line[1:].rstrip().split(' ')[0]
      if chrom in gen:
        sys.stderr.write('Error! Duplicate chromosome names: %s\n' % chrom)
        sys.exit(-1)
    else:
      seq += line.rstrip()
  if seq:
    gen[chrom] = seq
    length[chrom] = len(seq)
  f1.close()
  return gen, length

def createRegions(length):
  '''
  Create regions in [0, 1) for given chromosome lengths.
  '''
  total = float(sum(length.values()))
  reg = {}
  start = 0
  for chrom in length:
    val = length[chrom] / total
    reg[chrom] = [start, start + val]
    start += val
  reg[chrom][1] = 1  # ensure last segment ends at 1
  return reg

def loadBED(fBed):
  bed = {}
  for line in fBed:
    spl = line.rstrip().split('\t')
    if len(spl) < 3:
      sys.stderr.write('Warning! Not enough info in BED record:\n%s' % line)
      continue
    if len(spl) > 3:
      if spl[3] in bed:
        sys.stderr.write('Warning! Duplicate BED region names: %s\n' % spl[3])
        continue
      bed[spl[3]] = (spl[0], int(spl[1]), int(spl[2]))
    else:
      i = 0
      while 1:
        if i not in bed:
          bed[i] = (spl[0], int(spl[1]), int(spl[2]))
          break
        i += 1
  fBed.close()
  return bed

def createReads(f, gen, bed, reg, length, lengthDNA, lengthRead, number):

  sys.stderr.write('Creating reads...\n')

  randomProb = 0.01  # probability a random DNA fragment will be sequenced
  regionProb = 0.1   # probability a fragment overlapping a region will be sequenced

  # simulated reads
  for i in range(number):

    # find a suitable fragment
    while 1:

      # randomly choose chromosome (weighted by length)
      rand = random.random()
      chrom = ''
      for r in reg:
        if reg[r][0] <= rand < reg[r][1]:
          chrom = r
          break
      if not chrom:
        chrom = r

      # randomly choose position
      last = length[chrom] - lengthDNA  # last valid 5' position
      pos = random.randint(0, last)

      # determine if fragment lies within a BED region
      prob = randomProb  # probability of sequencing a fragment
      for b in bed:
        if bed[b][0] == chrom \
            and (bed[b][1] <= pos < bed[b][2] \
            or (pos < bed[b][1] and pos + lengthDNA >= bed[b][1])):
          prob = regionProb  # increase prob.
          break

      # print read
      if random.random() < prob:
        # randomly choose which end of fragment to sequence
        if random.random() < 0.5:
          seq = gen[chrom][pos:pos+lengthRead]
          strand = 'fwd'
        else:
          seq = revComp(gen[chrom][pos+lengthDNA-lengthRead:pos+lengthDNA])
          strand = 'rev'

        # skip if sequence has an 'N'
        if seq.find('N') != -1:
          continue

        f.write('@read' + ' '.join([str(i), chrom, str(pos), strand]) \
          + '\n' + seq.upper() \
          + '\n+\n' + 'I'*lengthRead + '\n')
        break

      # if read isn't printed, start over again

  f.close()

def main():
  args = sys.argv[1:]
  if len(args) < 6:
    sys.stderr.write('Usage: python %s  <gen>  <bed>  <out>' % sys.argv[0] \
      + '  <len1>  <len2>  <cov>\n' \
      + '  <gen>   Fasta genome\n' \
      + '  <bed>   BED file listing peak regions\n' \
      + '  <out>   Fastq output\n' \
      + '  <len1>  Length of DNA fragments\n' \
      + '  <len2>  Length of reads\n' \
      + '  <cov>   Number of reads\n')
    sys.exit(-1)

  fIn = open(args[0], 'rU')
  fBed = open(args[1], 'rU')
  fOut = open(args[2], 'w')
  lengthDNA = int(args[3])
  lengthRead = int(args[4])
  coverage = int(args[5])

  # produce reads
  gen, length = loadGenome(fIn)
  reg = createRegions(length)
  bed = loadBED(fBed)
  createReads(fOut, gen, bed, reg, length, lengthDNA, lengthRead, coverage)

if __name__ == '__main__':
  main()

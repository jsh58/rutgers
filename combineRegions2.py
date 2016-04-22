#!/usr/bin/python

# JMG 4/22/16
# Combining regions for a set of bismark cov files.

import sys

def valChr(k):
  '''
  Method to sort chromosomes.
  '''
  if k == 'chrM':
    return 0
  elif k == 'chrX':
    return 20
  elif k == 'chrY':
    return 21
  return int(k[3:])

def val(k):
  return int(k)


def processRegion(chr, reg, d, samples, fOut):
  #print 'Calling procReg for %s' % chr
  #for r in reg:
  #  print r,
  #print
  #sys.exit(0)

  fOut.write('%s\t%d\t%d\t%d' % (chr, reg[0], reg[-1], len(reg)))
  for sample in samples:
    meth = unmeth = 0
    for r in reg:
      pos = str(r)
      if sample in d[chr][pos]:
        meth += d[chr][pos][sample][0]
        unmeth += d[chr][pos][sample][1]
    if meth + unmeth == 0:
      fOut.write('\t0')
    else:
      fOut.write('\t%f' % (meth / float(meth+unmeth)))
  fOut.write('\n')

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    print 'Usage: python %s  <outfile>  <infile1>  <infile2> ...' % sys.argv[0]
    sys.exit(-1)

  fOut = open(args[0], 'w')

  d = {}
  tot = {}
  samples = []
  num = 1  # must have more than this number of reads
  for fname in args[1:]:
    try:
      fIn = open(fname, 'rU')
    except IOError:
      print 'Error! Cannot open', fname
      sys.exit(-1)

    # save sample name
    dot = fname.split('.')
    sample = dot[0]
    samples.append(sample)

    for line in fIn:
      spl = line.rstrip().split('\t')
      if len(spl) < 6:
        print 'Error! Poorly formatted cov record:\n', line
        sys.exit(-1)
      try:
        total = int(spl[4]) + int(spl[5])
      except ValueError:
        print 'Error! Poorly formatted cov record:\n', line

      # save counts and total
      if not spl[0] in d:
        d[spl[0]] = {}
        tot[spl[0]] = {}
      if not spl[1] in d[spl[0]]:
        d[spl[0]][spl[1]] = {}
      d[spl[0]][spl[1]][sample] = [int(spl[4]), int(spl[5])]
      if total > num:
        tot[spl[0]][spl[1]] = tot[spl[0]].get(spl[1], 0) + 1
    fIn.close()

  fOut.write('\t'.join(['chr', 'pos5', 'pos3', 'CpG'] + samples) + '\n')

  # save connected positions
  for chr in sorted(tot, key=valChr):
    reg = []
    pos3 = 0
    for pos in sorted(tot[chr], key=val):
      # require 2 samples
      if tot[chr][pos] > 1:
        loc = int(pos)
        if pos3 and loc - pos3 > 1000:
          processRegion(chr, reg, d, samples, fOut)
          reg = []
        reg.append(loc)
        pos3 = loc
    processRegion(chr, reg, d, samples, fOut)
  fOut.close()

#        chunk = [[] for sample in samples]
 #       rec = '\t'.join([chr, pos])
  #      for sample in samples:
   #       data = d[chr][pos].get(sample, None)
    #      if sample in d[chr][pos]:
     #       rec += '\t'.join(d[chr][pos][sample])
      #    else:
       #     rec += '\t0' * 2
        #fOut.write(rec + '\n')


if __name__ == '__main__':
  main()

#!/usr/bin/python

# JMG 6/2/16

# Finding diffs of specific types from diffMeth2.py output

import sys

def main():
  args = sys.argv[1:]
  if len(args) < 3:
    print 'Usage: python %s  <CSVfile>  <output>  (up|down)' % sys.argv[0]
    print '  <CSVfile>  Output from diffMeth2.py'
    print '  [up|down]  Report only results in this direction (def. up)'
    sys.exit(-1)
  try:
    fIn = open(args[0], 'rU')
  except IOError:
    print 'Error! Cannot open', args[0]
    sys.exit(-1)
  fOut = open(args[1], 'w')
  up = 1
  if args[2] == 'down':
    up = 0

  # write lines where 0->1 is > 0.1 (abs),
  #   and 0->2 is less than half 0->1
  count = 0
  fOut.write(fIn.readline())
  for line in fIn:
    spl = line.split('\t')
    if spl[-3] == 'NA' or spl[-2] == 'NA': continue
    val1 = float(spl[-3])
    val2 = float(spl[-2])
    if up:
      if val1 >= 0.1 and val2 < val1 / 2:
        fOut.write(line)
        count += 1
    else:
      if val1 <= -0.1 and val2 > val1 / 2:
        fOut.write(line)
        count += 1
  fIn.close()
  fOut.close()
  print 'Regions found:', count

if __name__ == '__main__':
  main()


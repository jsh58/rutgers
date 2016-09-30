#!/usr/bin/python

# JMG 9/27/16
# Reformat expression data to generate a heatmap.

import sys

def main():
  args = sys.argv[1:]
  if len(args) < 2:
    sys.stderr.write('Usage: python %s  <gene_exp.diff>  ' % sys.argv[0]
      + '<output>  [<base_sample>]  [excluded_samples_csv]\n')
    sys.exit(-1)
  fIn = open(args[0], 'rU')
  fOut = open(args[1], 'w')

  # load optional sample name to compare against (e.g. 'Control')
  sample = ''
  if len(args) > 2:
    sample = args[2]

  # load samples to be excluded from analysis
  excl = []
  if len(args) > 3:
    excl = args[3].split(',')

  # check header for issues
  header = fIn.readline().rstrip()
  spl = header.split('\t')
  if spl[4] != 'sample_1' or spl[5] != 'sample_2' \
      or spl[6] != 'status' or spl[9] != 'log2(fold_change)' \
      or spl[13] != 'significant':
    sys.stderr.write('Error! Input not properly formatted\n')
    sys.exit(-1)

  # parse input file
  val = {}     # for saving log2-fc values
  signif = {}  # is gene significant?
  maxVal = 0.0 # max fold-change (for converting 'inf')
  for line in fIn:
    spl = line.rstrip().split('\t')

    # check for excluded sample
    if excl and (spl[4] in excl or spl[5] in excl):
      continue
    #if sample and spl[4] != sample and spl[5] != sample:
    #  continue

    gene = spl[0]
    if gene not in val:
      val[gene] = []
      signif[gene] = 0
    if spl[13] == 'yes':
      signif[gene] = 1
      if spl[9][-3:] != 'inf' and abs(float(spl[9])) > maxVal:
        maxVal = abs(float(spl[9]))

    # save log2-fc values
    if spl[6] != 'OK':
      val[gene].append((spl[4], spl[5], '0'))  # put 0 if test not done
    else:
      val[gene].append((spl[4], spl[5], spl[9]))

  fIn.close()

  # print outputs for signif genes
  if not sample:
    for gene in sorted(val):
      if signif[gene]:
        for v in val[gene]:
          fOut.write('%s\t%s\t%s\t%s\n' % (gene, v[0], v[1], v[2]))
    fOut.close()
    sys.exit(0)

  # if given baseline sample name, report log2-fc
  #   differences with that sample only
  # first, print header
  order = []
  header = 'Symbol'
  for v in val[gene]:
    if v[0] == sample:
      header += '\t' + v[1]
      order.append(v[1])
    elif v[1] == sample:
      header += '\t' + v[0]
      order.append(v[0])
  if not order:
    sys.stderr.write('Error! Baseline sample %s not found\n' % sample)
    sys.exit(-1)
  fOut.write(header + '\n')

  # print values
  for gene in sorted(val):
    if not signif[gene]:
      continue

    # save result values to 'res' list
    res = [''] * len(order)
    for v in val[gene]:
      if v[0] == sample:
        try:
          idx = order.index(v[1])
        except ValueError:
          continue

        # convert 'inf' to maxVal
        if v[2][-3:] == 'inf':
          res[idx] = v[2].replace('inf', str(maxVal))
        else:
          res[idx] = v[2]

      elif v[1] == sample:
        # samples in reverse order
        try:
          idx = order.index(v[0])
        except ValueError:
          continue

        # convert 'inf' to maxVal
        if v[2] == '-inf':
          res[idx] = str(maxVal)
        elif v[2] == 'inf':
          res[idx] = '-' + str(maxVal)
        elif v[2]:
          res[idx] = str(-float(v[2]))  # negate b/c order is backwards
        else:
          res[idx] = ''

    fOut.write('%s\t' % gene + '\t'.join(res) + '\n')

  fOut.close()

if __name__ == '__main__':
  main()

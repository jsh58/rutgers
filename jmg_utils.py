#!/usr/bin/python

# JMG 6/28/16
# Wrappers for various functions.

#import sys
#import gzip

def getInt(arg):
  '''
  Convert given argument to int.
  '''
  try:
    val = int(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to int\n' % arg)
    sys.exit(-1)
  return val

def getFloat(arg, minVal = None, maxVal = None):
  '''
  Convert given argument to float. Ensure it is within range.
  '''
  try:
    val = float(arg)
  except ValueError:
    sys.stderr.write('Error! Cannot convert %s to float\n' % arg)
    sys.exit(-1)
  if (minVal != None and val < minVal) or \
      (maxVal != None and val > maxVal):
    sys.stderr.write('Error! Value %f is outside of range [%f,%f]\n' \
      % (val, minVal, maxVal))
    sys.exit(-1)
  return val

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

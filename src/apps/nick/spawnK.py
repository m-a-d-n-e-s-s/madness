#!/usr/bin/python
import sys, math
#
if len(sys.argv) != 4:
  sys.exit("Requires arguments: min max dx");
def floatRange(a, b, inc):
  """
  Returns a list containing an arithmetic progression of floats.
  This is simply a float version of the built-in range(a, b, step)
  function.  The result is [ , ) as always in Python.
  """
  try: x = [float(a)]
  except: return False
  for i in range(1, int(math.ceil((b - a ) / inc))):
    x. append(a + i * inc)
  return x

xMIN = float(sys.argv[1])
xMAX = float(sys.argv[2])
dx = float(sys.argv[3])
for y in floatRange(xMIN, xMAX, dx):
  for z in floatRange(xMIN, xMAX, dx):
    print "0 %3.1f %3.1f" % (y, z)

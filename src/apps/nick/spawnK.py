#!/usr/bin/python
import sys, math
#cart = False corresponds to spherical coordinates
cart = False
zaxis = False
eps = 0.001
if(cart):
  if(len(sys.argv) != 7):
    sys.exit("Requires arguments: xMIN xMAX dx zMIN zMAX dz");
if(not cart):
  if(len(sys.argv) != 5):
    sys.exit("Requires arguments: kMIN dk kMAX   nTH=PI/nTH");
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

if(cart):
  xMIN = float(sys.argv[1])
  xMAX = float(sys.argv[2])
  dx = float(sys.argv[3])
  zMIN = float(sys.argv[4])
  zMAX = float(sys.argv[5])
  dz = float(sys.argv[6])
  for x in floatRange(xMIN, xMAX+eps, dx):
    for z in floatRange(zMIN, zMAX+eps, dz):
      #if(math.sqrt(x*x + z*z)<min(xMAX,zMAX)):
      print "0 %3.1f %3.1f" % (x, z)
else:
  kMIN = float(sys.argv[1])
  dk = float(sys.argv[2])
  kMAX = float(sys.argv[3])
  PI = 3.14159265358979
  dTH = PI/float(sys.argv[4])
  for k in floatRange(kMIN, kMAX+eps, dk):
    if( math.fabs(k)<eps ):
      break
    for th in floatRange(0, PI+eps, dTH):
      print "0 %13.12f %13.12f" % (k*math.sin(th), k*math.cos(th))

# z-axis resolution (this messes up numerical integration)
# Add resolution to the z-oaxis
if(zaxis):
#  print "k%dk  "
  for k in floatRange(kMIN, kMAX+eps, dk/5):
#    print math.fabs((k + eps/100)%dk), "   "
    if( math.fabs((k + eps/100)%dk) > eps ):
      print "0 0 %5.4f" % k


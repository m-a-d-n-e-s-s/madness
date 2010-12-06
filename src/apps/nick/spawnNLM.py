#!/usr/bin/python
import sys
if(len(sys.argv) != 2):
    sys.exit("Requires nMAX as an argument");
for n in range(0,int(sys.argv[1])+1):
    for l in range (0, n):
        print n, l, "0"

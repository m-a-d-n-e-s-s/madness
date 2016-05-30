#!/usr/bin/env python3
"""Generate vectorized mtxm implementations."""

import logging
import argparse
import sys
from itertools import permutations, product

from codegen.driver import tester_gen
from codegen.bgqdriver import tester_gen_bgq
from codegen.mtxm import *

# TODO: call bgq version of tester gen if needed

def tester(m, versions):
    """Generate .cc files and Makefile to test correctness and
    performance of versions"""
    n = 1
    with open("Makefile", 'w') as mf:
        if m.have_bgp:
            print("\tHPM=/soft/apps/UPC/lib/libhpm.a", file=mf)
            print("\tCXX=tmpixlcxx_r", file=mf)
            print("\tCXXFLAGS=-g -O3 -qarch=450d -qtune=450 -qthreaded", file=mf)
        elif m.have_bgq:
            print("\tCXX=/bgsys/drivers/ppcfloor/comm/xl/bin/mpixlcxx_r", file=mf)
            print("\tCXXFLAGS=-I/bgsys/drivers/ppcfloor/ -O5 -qhot=level=1 -qsimd=auto -qtune=qp -qarch=qp", file=mf)
        else:
            print("\tCXX=g++", file=mf)
            march = "ssse3"
            extra = ''
            if type(m) == MTXMAVX:
                march = "avx"
            if type(m) == MTXMAVX2:
                march = "avx2"
                extra = "-march=core-avx2"
            print("\tCXXFLAGS=-O3 {} -m{} -lrt -lm".format(extra, march), file=mf)
        print("all:", file=mf)

        while versions:
            bunch_size = 1
            if m.have_bgp:
                bunch_size = 64
            elif m.have_bgq:
                bunch_size = 128
            cur, versions = versions[:bunch_size], versions[bunch_size:]
            with open("tune_mtxm_{}.cc".format(n), 'w') as f:
                if m.have_bgq:
                    tester_gen_bgq(f, cur, m.gen, m.complex_a, m.complex_b)
                else:
                    tester_gen(f, cur, m.gen, m.complex_a, m.complex_b, m.have_bgp)
            print("\t$(CXX) -c $(CXXFLAGS) tune_mtxm_{0}.cc -o tune_mtxm_{0}.o".format(n), file=mf)
            print("\t$(CXX) $(CXXFLAGS) tune_mtxm_{0}.o $(HPM) -o tune_mtxm_{0}.x".format(n), file=mf)
            if m.have_bgq:
                print("\tqsub -A MADNESS_MPQC_esp -t 10 -n {0} --mode c1 --cwd `pwd` -O tune_mtxm_{1} ./tune_mtxm_{1}.x".format(len(cur), n), file=mf)
            n += 1

def main():
    # Parse Arguments
    parser = argparse.ArgumentParser(description='What does this program do?')
    parser.add_argument('-l', '--log', dest='loglevel', default='ERROR',
            help='log level to output (default: ERROR)')
    parser.add_argument('--log-file', dest='logfile',
            help='destination file to write logging output')
    parser.add_argument('-a', '--complex-a', action='store_true', default=False,
            help='A (left) matrix is complex', dest='cxa')
    parser.add_argument('-b', '--complex-b', action='store_true', default=False,
            help='B (right) matrix is complex', dest='cxb')
    parser.add_argument('-m', '--arch', default='sse', choices=['sse', 'avx', 'avx2', 'bgp', 'bgq'],
            help='Target architecture')
    parser.add_argument('-n', '--name', default='mtxmq',
            help='Name of function to generate')
    parser.add_argument('-t', '--test', default=False, action='store_true',
            help='Generate performance and correctness testing code')
    args = parser.parse_args()

    try:
        loglevel = getattr(logging, args.loglevel.upper())
    except:
        parser.error("Invalid LOGLEVEL")

    # Configure logging
    logging.basicConfig(level=loglevel,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            filename=args.logfile)
    logger = logging.getLogger(__name__)

    if args.arch == 'sse':
        M = MTXMSSE
        bests = {
                (False, False): ("jik", {'i':2, 'j':14, 'k':1}, args.name),
                (True, False) : ("jik", {'i':2, 'j':12, 'k':1}, args.name),
                (False, True) : ("jik", {'i':2, 'j':14, 'k':1}, args.name),
                (True, True)  : ("ijk", {'i':2, 'j':14, 'k':1}, args.name)}
        versions = [(lo, {'i':i, 'j':j, 'k':1}, '{}_i{}j{}k1'.format(lo, i,j)) for i in range(2,21,2) for j in range(2,21,2) for lo in ["ijk", "jik"] if i*j <= 60]
    elif args.arch == 'avx':
        M = MTXMAVX
        bests = {
                (False, False): ("jik", {'i':2, 'j':20, 'k':1}, args.name), # jik_i3j12k1 is a bit better, except on 400 20 20 much slower
                (True, False) : ("ijk", {'i':4, 'j':12, 'k':1}, args.name), # ijk_i4j8k1 is better for 400 20 20 but worse on others
                (False, True) : ("ijk", {'i':3, 'j':12, 'k':1}, args.name), # jik_i3j12k1 best on avg across tested sizes, ijk_i4j8k1 best on 400 20 20
                (True, True)  : ("ijk", {'i':2, 'j':16, 'k':1}, args.name)} # jik is better for squares, ijk better for 400 20 20 and similar
        versions = [(lo, {'i':i, 'j':j, 'k':1}, '{}_i{}j{}k1'.format(lo, i,j)) for i in range(1,10) for j in range(4,21,4) for lo in ["ijk", "jik"] if i*j <= 60]
    elif args.arch == 'avx2':
        M = MTXMAVX2
        bests = {
                (False, False): ("jik", {'i':2, 'j':20, 'k':1}, args.name),
                (True, False) : ("jik", {'i':3, 'j':16, 'k':1}, args.name),
                (False, True) : ("jik", {'i':3, 'j':16, 'k':1}, args.name),
                (True, True)  : ("ijk", {'i':2, 'j':20, 'k':1}, args.name)}
        versions = [(lo, {'i':i, 'j':j, 'k':1}, '{}_i{}j{}k1'.format(lo, i,j)) for i in range(1,10) for j in range(4,21,4) for lo in ["ijk", "jik"] if i*j <= 60]
    elif args.arch == 'bgp':
        M = MTXMBGP
        bests = {
                (False, False): ("ijk", {'i':8, 'j':4, 'k':1}, args.name),
                (True, False) : ("ijk", {'i':6, 'j':8, 'k':1}, args.name),
                (False, True) : ("jik", {'i':2, 'j':20, 'k':1}, args.name), # ijk_i9j4k1 is better for (m*m,m)T*(m*m) 400 20 20 by almost 30%, but generally quite a bit worse
                (True, True)  : ("ijk", {'i':5, 'j':8, 'k':1}, args.name)}
        versions = [(lo, {'i':i, 'j':j, 'k':1}, '{}_i{}j{}k1'.format(lo, i,j)) for i in range(1,10) for j in range(4,21,4) for lo in ["ijk", "jik"] if i*j <= 60]
    elif args.arch == 'bgq':
        M = MTXMBGQ
        bests = {
                (False, False): ("jik", {'i':4, 'j':20, 'k':1}, args.name),
                (True, False) : ("ijk", {'i':6, 'j':8, 'k':1}, args.name),
                (False, True) : ("jik", {'i':4, 'j':20, 'k':1}, args.name),
                (True, True)  : ("jik", {'i':3, 'j':28, 'k':1}, args.name)}
        versions = [(lo, {'i':i, 'j':j, 'k':1}, '{}_i{}j{}k1'.format(lo, i,j)) for i in range(1,14) for j in range(4,29,4) for lo in ["ijk", "jik"] if i*j <= 90]

    m = M(args.cxa, args.cxb)
    best = bests[(args.cxa, args.cxb)]

    if args.test:
        tester(m, versions)
    else:
        m.gen(sys.stdout, *best)


if __name__ == '__main__':
    main()


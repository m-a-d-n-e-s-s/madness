#!/usr/bin/env python3

import sys
import subprocess
sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from madjsoncompare import madjsoncompare

if __name__ == "__main__":

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory: @SRCDIR@")

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile="@SRCDIR@/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments=' --geometry=he --wf=cis'
    dft_arguments=' --dft="maxiter=1; econv=1.e-4; dconv=1.e-3; prefix='+prefix+'"'
    other_arguments=' --response="thresh=1.e-3; maxiter=10; guess_maxiter=0; econv=1; dconv=1; guess_excitation_operators=dipole+; excitations=0"'
    cmd='rm '+outputfile+'; ./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)
    print("finished with run")
    print(p.stdout)
    exitcode=p.returncode
    print("exitcode ",exitcode)

    sys.exit(p.returncode)
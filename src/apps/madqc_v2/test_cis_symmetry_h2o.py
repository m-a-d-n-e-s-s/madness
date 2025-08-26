#!/usr/bin/env python3

import sys
import subprocess
sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import madjsoncompare, cleanup, skip_on_small_machines

if __name__ == "__main__":

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory: @SRCDIR@")

    if (skip_on_small_machines()):
        print("Skipping this verylong test on small machines")
        sys.exit(77)

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile="@SRCDIR@/"+prefix+"_b1.calc_info.ref.json"

    # run test
    global_arguments=' --molecule=h2o --wf=cis --prefix='+prefix
    dft_arguments=' --dft="k=8; localize=canon;"
    other_arguments=' --tdhf="freeze=1; thresh=1.e-3; econv=1.e-3; dconv=1.e-2"'
    cleanup(prefix)  # Clean up previous output files
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)

    print("finished with run")
    print(p.stdout)

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",1,"excitations",0,"irrep"],1.0)
    cmp.compare(["tasks",1,"excitations",0,"omega"],1.e-4)
    cmp.compare(["tasks",1,"excitations",0,"oscillator_strength_length"],1.e-3)
    cmp.compare(["tasks",1,"excitations",0,"oscillator_strength_velocity"],1.e-3)

    referencefile="@SRCDIR@/"+prefix+"_a2.calc_info.ref.json"
    dft_arguments=' --dft="k=8; localize=canon; prefix='+prefix+'; no_compute=1"'
    other_arguments=' --tdhf="irrep=a2; freeze=1; thresh=1.e-3; econv=1.e-3; dconv=1.e-2; restart=no_restart"'
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)
    print("finished with run")
    print(p.stdout)


    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",1,"excitations",0,"irrep"],1.0)
    cmp.compare(["tasks",1,"excitations",0,"omega"],1.e-4)
    cmp.compare(["tasks",1,"excitations",0,"oscillator_strength_length"],1.e-3)
    cmp.compare(["tasks",1,"excitations",0,"oscillator_strength_velocity"],1.e-3)
    print("final success: ",cmp.success)


    sys.exit(cmp.exitcode())

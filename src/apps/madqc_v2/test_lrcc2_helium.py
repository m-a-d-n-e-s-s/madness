#!/usr/bin/env python3

import sys
import subprocess
import argparse

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import madjsoncompare, cleanup, skip_on_small_machines

if __name__ == "__main__":

    # skip test on small machines
    # sys.exit(77) is used to indicate that the test was skipped, cf AddScriptedTests.cmake
    if (skip_on_small_machines()):
        print("Skipping this verylong test on small machines")
        sys.exit(77)

    # get command line arguments
    parser=argparse.ArgumentParser(description='command line arguments for this test case')
    # default value will be set by cmake
    parser.add_argument("--reference_directory",action="store",default="@SRCDIR@",help="the directory with the reference file in json format")
    args=parser.parse_args()

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory:",args.reference_directory)

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile=args.reference_directory+"/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments=' --geometry=he --wf=cc2'
    dft_arguments=' --dft="maxiter=10; econv=1.e-5; dconv=1.e-3; prefix='+prefix+'; k=5"'
    other_arguments=' --cc2="freeze 0; calc_type=lrcc2; iter_max=2"'
    cleanup(prefix)  # Clean up previous output files
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)

    p=subprocess.run(cmd,shell=True,stdout=None, stderr=subprocess.PIPE , universal_newlines=True)
    print("finished with run")
    exitcode=p.returncode
    print("exitcode ",exitcode)


    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",0,"model"],1.e-4)
    cmp.compare(["tasks",0,"properties","energy"],1.e-4)
    cmp.compare(["tasks",1,"model"],1.e-4)
    cmp.compare(["tasks",1,"nfreeze"],1.e-2)
    cmp.compare(["tasks",1,"cc2_correlation_energy"],1.e-2)
    cmp.compare(["tasks",1,"excitations",0,"omega"],1.e-2)  # lrcc2 es
    cmp.compare(["tasks",2,"mp2_correlation_energy"],1.e-2)
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode() + exitcode)
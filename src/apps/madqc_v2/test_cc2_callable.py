#!/usr/bin/env python3

import sys
import subprocess
import argparse

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from madjsoncompare import madjsoncompare

if __name__ == "__main__":

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
    dft_arguments=' --dft="maxiter=1; econv=1.e-4; dconv=1.e-3; prefix='+prefix+'; k=5"'
    other_arguments=' --cc2="freeze 1"'
    cmd='rm '+outputfile+' reference.00000; ./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)

    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)
    print("finished with run")
    print(p.stdout)
    exitcode=p.returncode
    print("exitcode ",exitcode)


    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",0,"model"],1.e-4)
    cmp.compare(["tasks",0,"properties","energy"],1.e-4)
    cmp.compare(["tasks",1,"model"],1.e-4)
    cmp.compare(["tasks",1,"nfreeze"],1.e-2)
    cmp.compare(["tasks",1,"mp2_correlation_energy"],1.e-2)
    print("final success: ",cmp.success)

    # sys.exit(cmp.exitcode())
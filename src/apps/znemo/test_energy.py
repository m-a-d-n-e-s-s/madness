#!/usr/bin/env python3

import sys
import subprocess
import argparse

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import madjsoncompare

if __name__ == "__main__":

    # get command line arguments
    parser=argparse.ArgumentParser(description='command line arguments for this test case')
    # default value will be set by cmake
    parser.add_argument("--reference_directory",action="store",default="@SRCDIR@",help="the directory with the reference file in json format")
    args=parser.parse_args()

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory:",args.reference_directory)

    prefix='madtest1'
    outputfile=prefix+'.calc_info.json'
    referencefile=args.reference_directory+"/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments=' --geometry=he'
    dft_arguments=' --dft="maxiter=1; econv=1.e-4; dconv=1.e-3; prefix='+prefix+'"'
    other_arguments=' --complex="physical_B=-1.0" --geometry="no_orient=true; source_name=he"'
    cmd='rm '+outputfile+' reference.00000 restartaodata; ./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)
    print("finished with run")

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["driver"],1.e-4)
    cmp.compare(["model"],1.e-4)
    cmp.compare(["return_energy"],1.e-2)
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode())
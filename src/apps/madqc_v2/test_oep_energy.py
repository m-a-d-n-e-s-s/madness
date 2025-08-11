#!/usr/bin/env python3

import sys
import subprocess
import argparse

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import madjsoncompare, cleanup

if __name__ == "__main__":

    # get command line arguments
    parser=argparse.ArgumentParser(description='command line arguments for this test case')
    # default value will be set by cmake
    parser.add_argument("--reference_directory",action="store",default="@SRCDIR@",help="the directory with the reference file in json format")
    args=parser.parse_args()

    # some user output
    print("Testing @BINARY@/@TESTCASE")
    print(" reference files found in directory:",args.reference_directory)

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile=args.reference_directory+"/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments=' --geometry=be --wf=oep'
    dft_arguments=' --dft="maxiter=3; econv=1.e-4; dconv=1.e-3; k=7; prefix='+prefix+'"'
    other_arguments=' --oep="model=oaep; oep_maxiter=3"'

    # cleanup previous output files
    cleanup(prefix)
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)

    print("finished with run")
    print(p.stdout)
    exitcode=p.returncode
    print("program ended successfully: ",exitcode==0)

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",0,"scf_total_energy"],1.e-4)
    cmp.compare(["tasks",0,"properties","energy"],1.e-4)
    cmp.compare(["tasks",0,"scf_eigenvalues_a","vals",0],1.e-4)
    cmp.compare(["tasks",1,"scf_total_energy"],1.e-4)
    cmp.compare(["tasks",1,"properties","energy"],1.e-4)
    cmp.compare(["tasks",1,"scf_eigenvalues_a","vals",0],1.e-4)
    cmp.compare(["tasks",1,"model"],1.e-4)
    cmp.compare(["tasks",1,"scf_total_energy"],1.e-4)
    cmp.compare(["tasks",1,"drho"],1.e-4)
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode() + exitcode)
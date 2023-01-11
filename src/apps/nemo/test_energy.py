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
    print("Testing nemo/test_localization")
    print(" reference files found in directory:",args.reference_directory)

    prefix='madtest1'
    outputfile=prefix+'.scf_info.json'
    referencefile=args.reference_directory+"/"+prefix+".scf_info.ref.json"

    # run test
    cmd='rm '+outputfile+'; nemo --input=input --geometry=he --dft="maxiter=1; econv=1.e-4; dconv=1.e-3; prefix='+prefix+'"'
    print("executing \n ",cmd)
    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    print("finished with run")

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare([0,"scf_k"],1.e-4)
    cmp.compare([0,"scf_threshold"],1.e-4)
    cmp.compare([0,"scf_energy"],1.e-2)
    cmp.compare([0,"scf_energy"],1.e-4)
    cmp.compare([0,"scf_dipole_moment","dims"],1.e-4)
    cmp.compare([0,"scf_dipole_moment","vals",0],1.e-4)
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode())
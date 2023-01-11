#!/usr/bin/env python3

import sys
import subprocess
import argparse

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from madjsoncompare import madjsoncompare



def localizer_run(localizer):
    prefix='madtest'+localizer
    outputfile=prefix+'.scf_info.json'

    refdir=args.reference_directory
    referencefile=refdir+"/"+prefix+".scf_info.ref.json"


    cmd='rm '+outputfile+'; nemo --input=input --geometry=h2o --dft="maxiter=10; econv=3.e-5; k=8; localize='+localizer+'; dconv=1.e-3; prefix='+prefix+'"'
    print("executing \n ",cmd)
    output=subprocess.run(cmd,shell=True,capture_output=True, text=True).stdout
    print("finished with run1")

    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare([0,"scf_k"],1.e-4)
    cmp.compare([0,"scf_threshold"],1.e-4)
    cmp.compare([0,"scf_energy"],1.e-2)
    cmp.compare([0,"scf_energy"],1.e-4)
    cmp.compare([0,"scf_dipole_moment","dims"],1.e-4)
    cmp.compare([0,"scf_dipole_moment","vals",0],1.e-4)
    return cmp.success

if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='command line arguments for this test case')
    # default value will be set by cmake
    parser.add_argument("--reference_directory",action="store",default="@SRCDIR@",help="the directory with the reference file in json format")
    args=parser.parse_args()


    print("Testing nemo/test_localization")
    print(" reference files found in directory:",args.reference_directory)
    success=localizer_run('canon')
    success=localizer_run('boys') and success
    success=localizer_run('new') and success

    print("final success: ",success)
    ierr=0
    if (not success):
        ierr=1
    sys.exit(ierr)
#!/usr/bin/env python3

import sys
import subprocess
sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import madjsoncompare, cleanup

if __name__ == "__main__":

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory: @SRCDIR@")

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile="@SRCDIR@/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments='--molecule=he --wf=nemo --prefix='+prefix
    dft_arguments=' --dft="maxiter=10; econv=1.e-5; dconv=1.e-3;"'
    other_arguments=''
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    cleanup(prefix)
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)
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
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode() + exitcode)

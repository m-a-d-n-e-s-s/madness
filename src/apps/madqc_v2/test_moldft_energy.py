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
    global_arguments=' --molecule="source_name=he; eprec=1.e-6" --prefix='+prefix
    dft_arguments=' --dft="maxiter=10; econv=1.e-4; dconv=1.e-3;"'
    other_arguments=''
    cleanup(prefix)  # Clean up previous output files
    cmd='./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)

    print("finished with run")
    print(p.stdout)
    exitcode=p.returncode
    print("exitcode ",exitcode)

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",0,"scf_total_energy"],1.e-4)
    print("final success: ",cmp.success)

    sys.exit(p.returncode + exitcode)
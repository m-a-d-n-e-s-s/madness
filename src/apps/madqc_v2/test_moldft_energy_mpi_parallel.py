#!/usr/bin/env python3

import sys
import subprocess
sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from madjsoncompare import madjsoncompare

def cleanup(prefix):
    """Remove output files and directories created during the test."""
    cmd = f'rm -r {prefix}.calc_info.json {prefix}'
    print("Cleaning up with command:", cmd)
    subprocess.run(cmd, shell=True)

if __name__ == "__main__":

    # some user output
    print("Testing @BINARY@/@TESTCASE@")
    print(" reference files found in directory: @SRCDIR@")

    prefix='mad_@BINARY@_@TESTCASE@'
    outputfile=prefix+'.calc_info.json'
    referencefile="@SRCDIR@/"+prefix+".calc_info.ref.json"

    # run test
    global_arguments=' --geometry=he'
    dft_arguments=' --dft="maxiter=1; econv=1.e-4; dconv=1.e-3; prefix='+prefix+'"'
    other_arguments=''
    cleanup(prefix)  # Clean up previous output files
    cmd='MAD_NUM_THREADS=1 mpirun -np 2 ./@BINARY@ '+global_arguments + dft_arguments  + other_arguments
    print("executing \n ",cmd)
#    p=subprocess.run(cmd,shell=True,capture_output=True, text=True)
    p=subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE , universal_newlines=True)

    print("finished with run")
    print(p.stdout)
    exitcode=p.returncode
    print("exitcode ",exitcode)

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["tasks",0,"properties","energy"],1.e-1)
    cmp.compare(["tasks",0,"metadata","mpi_size"],1.e-4)
    print("final success: ",cmp.success)

    sys.exit(p.returncode + exitcode)
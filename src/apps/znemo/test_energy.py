#!/usr/bin/env python3

import sys
import subprocess
import argparse
import os

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

    # cleanup leftovers from previous runs
    for stale in (outputfile, "reference.00000", "restartaodata"):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    # run test
    cmd = [
        "./@BINARY@",
        "--molecule=he",
        "--no_orient=true",
        f"--dft=maxiter=1; econv=1.e-4; dconv=1.e-3; prefix={prefix}",
        "--complex=physical_B=-1.0",
    ]
    env = os.environ.copy()
    try:
        env["MAD_NUM_THREADS"] = str(min(64, max(1, int(env.get("MAD_NUM_THREADS", "8")))))
    except ValueError:
        env["MAD_NUM_THREADS"] = "8"
    print("using MAD_NUM_THREADS =", env["MAD_NUM_THREADS"])
    print("executing \n ", " ".join(cmd))
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8', errors='replace',
        env=env,
    )
    print("finished with run")
    print(p.stdout)
    if p.stderr:
        print(p.stderr)
    if p.returncode != 0:
        print("znemo exited with non-zero status:", p.returncode)
        sys.exit(p.returncode)
    if not os.path.exists(outputfile):
        fallback = "mad.calc_info.json"
        if os.path.exists(fallback):
            print("Expected output file not found:", outputfile)
            print("Using fallback output file:", fallback)
            outputfile = fallback
        else:
            print("Expected output file not found:", outputfile)
            sys.exit(1)

    # compare results
    cmp=madjsoncompare(outputfile,referencefile)
    cmp.compare(["driver"],1.e-4)
    cmp.compare(["model"],1.e-4)
    cmp.compare(["return_energy"],1.e-2)
    print("final success: ",cmp.success)

    sys.exit(cmp.exitcode())

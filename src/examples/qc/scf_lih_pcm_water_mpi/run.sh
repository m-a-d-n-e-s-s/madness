#!/bin/sh
# Same launch conventions as scf_he_hf_mpi/run.sh: MAD_NUM_THREADS is the worker
# count PER RANK (plus one comm thread per rank), and --bind-to none keeps
# OpenMPI from pinning every thread of a rank to one core.
exec env MAD_NUM_THREADS=${MPI_WORKERS:-3} \
    ${MPIEXEC:-mpiexec} ${MPIEXEC_BIND---bind-to none} -np ${NP:-2} \
    ${MADQC:-madqc} --wf=scf scf_lih_pcm_water_mpi.in

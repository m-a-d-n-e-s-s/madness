#!/bin/sh
# Thread budget: MAD_NUM_THREADS is workers PER RANK, and with nranks > 1 each
# rank additionally spawns a MADNESS comm thread. 2 ranks x (3 workers + 1 comm)
# fills an 8-hwthread machine exactly. Starving a rank (MAD_NUM_THREADS=1) makes
# the one worker service both compute and communication, which is very slow.
#
# --bind-to none is load-bearing: OpenMPI binds each rank to a single core by
# default, so those threads land on one hwthread. The symptom is not slowness but
# a stall -- the SCF wedges in SCF::loadbal -> redistribute -> fence, burning CPU
# and making no progress. Other MPI vendors spell this differently or not at all;
# override MPIEXEC_BIND (set it empty to drop the flag).
# Deliberately NOT defaulted from an inherited MAD_NUM_THREADS: that variable
# conventionally means "threads for this job", whereas here it is per rank, so
# inheriting it would launch NP x that many workers and oversubscribe. Set
# MPI_WORKERS to change the per-rank count.
exec env MAD_NUM_THREADS=${MPI_WORKERS:-3} \
    ${MPIEXEC:-mpiexec} ${MPIEXEC_BIND---bind-to none} -np ${NP:-2} \
    ${MADQC:-madqc} --wf=scf scf_he_hf_mpi.in

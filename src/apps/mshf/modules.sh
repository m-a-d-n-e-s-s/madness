module load gcc
module load mvapich2
module load mkl

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/projects/hpcsoft/toss2/wolf/mvapich2/2.2_gcc-6.4.0/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/projects/hpcsoft/toss2/common/intel-clusterstudio/2017.1.024/compilers_and_libraries_2017/linux/mkl/lib/intel64
export POOL_NTHREAD=16
export OMP_NUM_THREADS=1
export MAD_WAIT_TIMEOUT=600
export MAD_BIND="1 0 2"
export MV2_ENABLE_AFFINITY=0 


###module load intel-mpi
###export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/projects/hpcsoft/toss3/snow/mvapich2/2.2_gcc-6.4.0/lib
###export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/projects/hpcsoft/toss3/common/intel-clusterstudio/2017.1.024/impi/2017.1.132/intel64/lib
###export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/projects/hpcsoft/toss3/common/intel-clusterstudio/2017.1.024/compilers_and_libraries_2017/linux/mkl/lib/intel64

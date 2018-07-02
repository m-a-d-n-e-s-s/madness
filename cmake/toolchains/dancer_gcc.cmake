# Set the system name so CMake uses the appropriate platform settings.
# NOTE: The platforms settings for BlueGeneP are the same for BlueGeneQ 
# set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# # Set environment paths
# set(IBM_DIR    "$ENV{IBM_MAIN_DIR}")
# set(XLF_DIR    "${IBM_DIR}/xlf/bg/14.1")
# set(XLSMP_DIR  "${IBM_DIR}/xlsmp/bg/3.1")
# set(ESSL_DIR   "/soft/libraries/essl/current/essl/5.1")
# set(LAPACK_DIR "/soft/libraries/alcf/current/xl/LAPACK")

#set (PARSEC_DIR "/home/aguermou/devel_build_gcc/install")

# # V1R2M0
# #set(MPI_DIR   "/bgsys/drivers/ppcfloor/comm/gcc")
# #set(PAMI_DIR  "/bgsys/drivers/ppcfloor/comm/sys")
# # V1R2M1
# #set(GCC_DIR    "/bgsys/drivers/toolchain/V1R2M2_base_4.7.2/gnu-linux-4.7.2")
# # V1R2M2
#set(GCC_DIR    "/bgsys/drivers/toolchain/V1R2M2_base_4.7.2-efix14/gnu-linux-4.7.2-efix014")
# set(MPI_DIR    "/bgsys/drivers/ppcfloor/comm")
# set(PAMI_DIR   "/bgsys/drivers/ppcfloor/comm")
# set(SPI_DIR    "/bgsys/drivers/ppcfloor/spi")

# Set compilers
set(CMAKE_C_COMPILER       "/opt/gcc-5.1/bin/gcc")
set(CMAKE_CXX_COMPILER     "/opt/gcc-5.1/bin/g++")
set(CMAKE_Fortran_COMPILER "/opt/gcc-5.1/bin/gfortran")
set(MPI_C_COMPILER         "mpicc")
set(MPI_CXX_COMPILER       "mpicxx")

# Set compile flags
set(CMAKE_C_FLAGS_INIT             "-std=c99" CACHE STRING "Inital C compile flags")
set(CMAKE_C_FLAGS_DEBUG            "-g -Wall" CACHE STRING "Inital C debug compile flags")
set(CMAKE_C_FLAGS_MINSIZEREL       "-Os -DNDEBUG" CACHE STRING "Inital C minimum size release compile flags")
set(CMAKE_C_FLAGS_RELEASE          "-O3 -DNDEBUG" CACHE STRING "Inital C release compile flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g -Wall" CACHE STRING "Inital C release with debug info compile flags")
set(CMAKE_CXX_FLAGS_INIT           " -m64 -mcx16 -g -Wall" CACHE STRING "Inital C++ compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall" CACHE STRING "Inital C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" CACHE STRING "Inital C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" CACHE STRING "Inital C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall" CACHE STRING "Inital C++ release with debug info compile flags")


# Set library

#set(XL_LIBRARIES ${XLSMP_DIR}/bglib64/libxlsmp.a)
#set(XLF_LIBRARIES ${XLF_DIR}/bglib64/libxlf90_r.a;${XLF_DIR}/bglib64/libxlfmath.a;${XLF_DIR}/bglib64/libxlopt.a;${XLF_DIR}/bglib64/libxl.a;-ldl;-lm)
set(BLAS_LIBRARIES /opt/intel/mkl/lib/intel64/libmkl_sequential.so;/opt/intel/mkl/lib/intel64/libmkl_core.so;/opt/intel/mkl/lib/intel64/libmkl_gf_lp64.so;-lm)
set(LAPACK_LIBRARIES ${BLAS_LIBRARIES})
set(LAPACK_COMPILE_DEFINITIONS HAVE_INTEL_MKL=1 CACHE STRING "LAPACK preprocessor definitions")
set(FORTRAN_INTEGER_SIZE "4" CACHE STRING "Set Fortran integer size in bytes")
set(HAVE_SPINLOCKS OFF CACHE BOOL "Enable if pthread lib supports spinlocks in pmrrr")

set(PAPI_INCLUDE_DIR "/opt/papi-5.4.3/include")
set(PAPI_LIBRARY "/opt/papi-5.4.3/lib/libpapi.so")


##############################################################

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
#set(CMAKE_FIND_ROOT_PATH
#    /bgsys/drivers/ppcfloor/
#    ${MPI_DIR}
#    ${PAMI_DIR}
#    ${SPI_DIR}
#    ${GCC_DIR}
#    ${CLANG_DIR}
#    ${IBM_DIR}
#    ${XLF_DIR}
#    ${XLSMP_DIR}
#    ${ESSL_DIR})

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
#set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
#set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

##############################################################

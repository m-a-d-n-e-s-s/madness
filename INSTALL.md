## Introduction

Linux and MacOS are supported with x86, Arm64, and IBM Power processors. GPUs are not yet utilized.

MADNESS uses CMake to configure the build. Assuming that necessary prerequisites (below) are installed on your system in default locations and the source has been downloaded into the directory `/path/to/madness/source`, you can make a directory (outside the source tree) to build in and configure the build as follows
```
mkdir build
cd build
cmake /path/to/madness/source
make applications
```
The default make target builds only the numerical library and underlying runtime.  To build applications (e.g., `moldft`, `nemo`) specify either `applications` (for all of them) or the name of the desired application on the make command.  The target `everything` does what you expect. You can run executables and use libraries from the build directory, but to install into the default location (`/usr`) use `make install` (after first building the applications).

If required libraries are not in default locations or if you wish to override defaults, you may have to set CMake cache variables as described below.  For instance, to build a debug version, without MPI, and with installation in `/home/me/madinstall`
```
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_MPI=OFF -DCMAKE_INSTALL_PREFIX=/home/me/madinstall /path/to/madness/source
make applications
make install
```


Boolean values for cache variables (specified to CMake using the `-DVARIABLE_NAME` notation in the example above) are 
considered true if the constant is 1, ON, YES, TRUE, Y, or a non-zero number; or false if the constant is 0, OFF, NO, FALSE, N, or IGNORE.

## Prerequisites

Fast BLAS and linear algebra libraries are essential. These must be sequential (single thread) implementations since MADNESS uses tasks/threads for parallelism and invokes the BLAS within a single-threaded task.  On X86, we recommend the free [Intel MKL library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.8bsxug), which is auto detected on all platforms if the environment variable `MKLROOT` is set.  On MacOS, the Apple [Accelerate](https://developer.apple.com/documentation/accelerate) framework (installed as part of [Xcode](https://developer.apple.com/xcode/)) is also autodetected.  AMD ACML has not been tested in a while but can be enabled with the CMake variables below.  Other libraries need to have their link (and possibly also header) paths and flags provided in the CMake variables for the compiler and linker.  On ARM, we recommend the ARM performance library available with their optimized LLVM compiler and enabled with the `-armpl`` compiler+linker flag. The OpenBLAS libary for ARM (not on x86) had a major peformance problem when we last tested it in circa 2020 due to a mutex around a memory block shared by all threads (this issue was reported and might be fixed by now).

An MPI library is required by default, and should be autodected primarily by looking for the `mpicxx` command to compile C++ code with MPI.  If you wish to overide this, use the appropriate cmake or environment variables (below).  You can also disable use of MPI using `-DENABLE_MPI=OFF` --- in this case you can still use all cores within a shared memory computer.

MADNESS will autodetect the [Intel TBB library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html#gs.8byhgg), which is available for free from Intel or via standard Linux package managers (even on ARM).  TBB provides a fast task pool.  If this is not detected, MADNESS will employ its own task pool.  [PaRSEC](https://icl.utk.edu/parsec/) can also be used (see variables below) but is not recommended unless you are using MADNESS with the [Template Task Graph](https://github.com/TESSEorg/ttg).

### Minimal Ubuntu packages 

If you are using a different distribution, just change the name of the command to install packages.

From a vanilla desktop install of Ubuntu (tested with 23.0.4) 
```
   sudo apt install build-essential cmake gcc g++ git intel-mkl mpich 
```
If you want to build the documentation then also
```
   sudo apt install doxygen graphviz texlive
```

### Minimal MacOS packages

To be written


The below summarizes some of the most useful CMake variables.

## Compiler variables

The following CMake cache variables are used to specify the compilers, compiler
flags, and linker flags.

* CMAKE_C_COMPILER --- C compiler [default=search]
* CMAKE_CXX_COMPILER --- C++ compiler [default=search]
* CMAKE_ASM_COMPILER --- Assembly compiler [default=search]
* MPI_CXX_COMPILER --- MPI C++ compiler wrapper [default=search]
* MPI_C_COMPILER --- MPI C compiler wrapper [default=search]
      
You can specify compile flags with the following variables. These cache variables
are optional, and it is typically not necessary to specify these variables.

* CMAKE_C_FLAGS --- Compile flags passed to the C compiler
* CMAKE_CXX_FLAGS --- Compile flags passed to the C++ compiler
* CMAKE_EXE_LINKER_FLAGS --- Linker flags to be used to create executables.
* CMAKE_STATIC_LINKER_FLAGS --- Linker flags to be used to create static 
      libraries.
* CMAKE_SHARED_LINKER_FLAGS --- Linker flags to be used to create shared 
      libraries.

## Build options

The following CMake cache variables turn features on and off.

* CMAKE_BUILD_TYPE --- Debug or Release
* ENABLE_GENTENSOR --- Enable generic tensors; only useful if need
                       compressed 6-d tensors, e.g. in MP2 [default=OFF]
* ENABLE_TASK_PROFILER - Enable task profiler that collects per-task start and 
      stop times. [default=OFF]
* ENABLE_WORLD_PROFILE --- Enables world profiling [default=OFF]
* ENABLE_MEM_STATS --- Gather memory statistics (expensive) (default=OFF)
* ENABLE_TENSOR_BOUNDS_CHECKING --- Enable checking of bounds in tensors ... 
      slow but useful for debugging [default=OFF]
* ENABLE_TENSOR_INSTANCE_COUNT --- Enable counting of allocated tensors for 
      memory leak detection [default=OFF]
* ENABLE_SPINLOCKS --- Enables use of spinlocks instead of mutexes (faster 
      unless over subscribing processors) [default=ON]
* ENABLE_NEVER_SPIN --- Disables use of spinlocks (notably for use inside
      virtual machines [default=OFF]
* ENABLE_BSEND_ACKS --- Use MPI Send instead of MPI Bsend for huge message 
      acknowledgements [default=ON]
* BUILD_TESTING --- Enables unit tests targets [default=ON]
* FORTRAN_INTEGER_SIZE --- Set the fortran integer size (4 or 8 bytes) used for 
      BLAS and LAPACK function calls [default=4]
* ASSERTION_TYPE --- Define MADNESS assertion behavior
      (abort|assert|disable|throw) [default=throw]
* MPI_THREAD --- Thread level for MPI (multiple|serialized)" [default=multiple]
* BUILD_SHARED_LIBS --- Enable shared libraries. This option is only available
      if the platform supports shared libraries; if that's true and MADNESS_ASSUMES_ASLR_DISABLED is ON (see below) the default is ON,
      otherwise the default is OFF.
* MADNESS_BUILD_MADWORLD_ONLY --- whether to build the MADNESS runtime only; if `ON`, discovery of BLAS/LAPACK
      and building of numerical components and applications will be disabled [default=`OFF`]
* MADNESS_BUILD_LIBRARIES_ONLY --- whether to build the MADNESS libraries only; if `ON`,
      building of numerical components and applications will be disabled and
      the value of `MADNESS_BUILD_MADWORLD_ONLY` ignored [default=`OFF`]

## External libraries

The following CMake cache variables enable the use of external libraries with
MADNESS. If the WITH_* variable is set to "ON" by default, failure to find the
external library is not an error. If you explicitly set a WITH_* variable to 
"ON" when the default is set to "OFF," an error will occur if the library is
not found.

* CMAKE_PREFIX_PATH - A semicolon seperated list of paths that are used when 
      searching for external libraries and dependencies. You may use this CMake
      cache variable to specify the prefix for any of the dependencies, or you
      may specify path for individual components below.

* ENABLE_MPI --- Enable use of MPI, should specify MPI_CXX_COMPILER or MPI_C_COMPILER
                 explicitly or have them in PATH; if not found will use StubMPI and
                 be limited to 1 process [default=ON]

In the following section, each optional library privides four variables that
the user can use to enable cmake to find the correct dependencies: 
  * ENABLE_<LIB> --- Enable the library feathres (ON|OFF)
  * <LIB>_ROOT_DIR --- Prefix path used to search for the external library.
  * <LIB>_INCLUDE_DIR --- The external library include path. By default the
        include path is ${<LIB>_ROOT_DIR}/include, if <LIB>_ROOT_DIR is
        specified in the configure command.
  * <LIB>_LIBRARY --- The external library path. By default the
        include path is ${<LIB>_ROOT_DIR}/(lib64|lib), if <LIB>_ROOT_DIR is
        specified in the configure command.
If the <LIB>_ROOT_DIR, <LIB>_INCLUDE_DIR, and <LIB>_LIBRARY will be used to
search for specific dependencies. If the external library is not found in these
given paths, or if the paths are not given, CMake will search the paths in 
CMAKE_PREFIX_PATH as well as other system paths.

### Library of Exchange-Correlation DFT functionals (LIBXC):

* ENABLE_LIBXC --- Enables use of the libxc library of density functionals.
      [default=ON]
* LIBXC_ROOT_DIR --- The install prefix for LIBXC.
* LIBXC_INCLUDE_DIR --- The path to the LIBXC include directory.
* LIBXC_LIBRARY --- The path to the LIBXC library.

E.g.,
```
cmake -DENABLE_LIBXC=ON -DLIBXC_LIBRARIES=/home/rjh/install/lib/libxc.a -DLIBXC_INCLUDE_DIRS=/home/rjh/install/include ../madness
```

### Intel Threading Building Blocks (TBB):

* Define MADNESS_TASK_BACKEND=TBB --- this should auto detect TBB if it is installed in a standard location.  If it is not, you'll have to also define
  * TBB_ROOT_DIR --- The install prefix for TBB.  If TBB_ROOT_DIR is not given, it will be set to the value of the TBBROOT environment variable if it is set.  And if `cmake` still does not detect things you can set
  * TBB_INCLUDE_DIR --- The path to the TBB include directory
  * TBB_LIBRARY --- The path to the TBB library directory. By default, the library
      search path is ${TBB_ROOT_DIR}/(lib/intel64/gcc4.4|lib) on Linux and
      ${TBB_ROOT_DIR}/(lib/libc++|lib) on OS X, if TBB_ROOT_DIR is specified in
      the configure command.
  * Optionally, MADNESS_EXPLOIT_TBB_PRIORITY --- If ON, MADNESS will try to use Intel TBB task priorities [default=OFF]


### Intel Math Kernel Library (MKL):

* ENABLE_MKL --- Search for Intel MKL for BLAS and LAPACK support [default=ON]
* MKL_ROOT_DIR --- The install prefix for MKL.
* MKL_LIBRARY --- The path to the MKL library directory.

If MKL_ROOT_DIR is not given, it will be set to the value of the MKLROOT environment variable if it is set.

### AMD Core Math Library (ACML):

*This is out of date since ACML is now AOCL.*

* ENABLE_ACML --- Search for AMD math library for BLAS and LAPACK support
      [default=ON]
* ACML_ROOT_DIR --- The install prefix for ACML.
* ACML_LIBRARY --- The path to the ACML library directory.

### Google Performance Tools (Gperftools):

* ENABLE_GPERFTOOLS --- Enable use of gperftools, including tcmalloc.
      [default=OFF]
* ENABLE_TCMALLOC_MINIMAL --- Enable use of the minimal tcmalloc library only.
      [default=OFF]
* GPERFTOOLS_ROOT_DIR --- The install prefix for gperftools.
* GPERFTOOLS_INCLUDE_DIR --- The path to the gperftools include directory.
* GPERFTOOLS_LIBRARY --- The path to the gperftools library directory.

If GPERFTOOLS_ROOT_DIR is not given, it will be set to the value of the GPERFTOOLS_DIR environment variable if it is set.

### Libunwind:

* ENABLE_LIBUNWIND --- Force detection of gperftools [default=OFF, i.e. Libunwind will be searched for when needed]
* LIBUNWIND_DIR --- The install prefix for Libunwind.

If LIBUNWIND_DIR is not given, it will be set to the value of the LIBUNWIND_DIR environment variable if it is set.

### IntegratorXX for numerical integration via DFT grids
* ENABLE_INTEGRATORXX --- Enables use of IntegratorXX
* INTEGRATORXX_ROOT_DIR --- The install prefix for IntegratorXX
* INTEGRATORXX_INCLUDE_DIR --- The path to the IntegratorXX include directory (should be added automatically when the correct PCM_ROOT_DIR is given)

IntegratorXX is a library for numerical integration via DFT grids. It is used in the MP3 code of Madness for 
generating low-rank representations of 6D functions. 
If IntegratorXX absent, a Gaussian-distributed random grid will be used, leading to slightly varying results in 
different MP3 runs.


### Polarizable Conitinuum Solver (PCM):

* ENABLE_PCM --- Enables use of PCM
* PCM_ROOT_DIR --- The install prefix for PCM 
* PCM_INCLUDE_DIR --- The path to the PCM include directory (should be added automatically when the correct PCM_ROOT_DIR is given)
* PCM_LIBRARY --- The path to the PCM library (should be added automatically when the correct PCM_ROOT_DIR is given)
set either PCM_ROOT_DIR or manually set PCM_INCLUDE_DIR and PCM_LIBRARY
See also
madness/CMakeLists.txt
madness/external/pcm.cmake
madness/modules/FindPCM.cmake
madness/src/apps/chem/CMakeLists.txt

### Performance Application Programming Interface (PAPI):

* ENABLE_PAPI --- Enables use of PAPI [default=OFF]
* PAPI_ROOT_DIR --- The install prefix for PAPI.
* PAPI_INCLUDE_DIR --- The path to the PAPI include directory.
* PAPI_LIBRARY --- The path to the PAPI library directory.

### Elemental parallel linear algebra library:

*This has not been tested in some time.*
      
Elemental provides optional distributed-memory linear algebra for some MADNESS application codes.
MADNESS source includes (modified) Elemental v0.84, which has been validated to work with
the few MADNESS apps that can use Elemental. You can instruct MADNESS to download and compile
a more recent version of Elemental, if desired, but the apps will not use Elemental then.
Such bundling is currently necessary if your code uses the MADworld runtime AND Elemental because
madness::initialize will call El::initialize() .

* ENABLE_ELEMENTAL --- Enable Elemental [default=OFF].
* ELEMENTAL_TAG --- specifies the version of Elemental to be downloaded, compiled, and installed alongside
                    MADNESS (numerical codes of MADNESS will not use Elemental).
                    If not set, will use the included Elemental source.

### Parallel Runtime Scheduling and Execution Controller (PaRSEC):

Recommended only for TTG development.

* ENABLE_PARSEC --- Enables use of PaRSEC as the task scheduler [default=OFF]. The use of Intel TBB should be disabled
                    to use PaRSEC.

If ENABLE_PARSEC is set but PaRSEC is not found, it will be built from source.
    
## MADNESS Runtime and the Address Space Layout Randomization (ASLR)

If you can run with one process and are getting almost immediate segmentation violations with more processes, this might be the issue.

ASLR (Linux relevant documentation [here](https://linux-audit.com/linux-aslr-and-kernelrandomize_va_space-setting/))
is a standard technique for increasing platform security implemented by the OS kernel and/or
the dynamic linker. By randomizing both where the shared libraries are loaded as well as (when enabled) the absolute
position of the executable in memory (such executables are known as position-independent executables). Until recently
MADNESS could only be used with MPI on platforms with ASLR if built with static libraries (for MADNESS code; system and other libraries could still be shared).  However, static libraries make it hard to integrate with Python and other frameworks that demand shared libraries. 

If properly configured and built, MADNESS can now be used on ASLR platforms using either static (the default and easiest) or shared libraries. Use the following variables to control the ASLR-related aspects of MADNESS runtime.

* MADNESS_ASSUMES_ASLR_DISABLED --- MADNESS runtime will assume that the Address Space Layout Randomization (ASLR) is off.
      By default MADNESS_ASSUMES_ASLR_DISABLED is set to OFF (i.e. MADNESS will assume that ASLR is enabled);
      this will cause all libraries by default to be static (BUILD_SHARED_LIBS=OFF)
      and compiled as position-independent code (CMAKE_POSITION_INDEPENDENT_CODE=ON).
      This will also enable a runtime check for ASLR.
* CMAKE_POSITION_INDEPENDENT_CODE --- This standard CMake variable controls whether targets are compiled by default
      as position-independent code or not. If BUILD_SHARED_LIBS=OFF need to set this to ON if want to use the MADNESS
      libraries to build shared libraries or position-independent executables.

To make things more concrete, consider the following 2 scenarios:
* Platform with ASLR disabled --- set MADNESS_ASSUMES_ASLR_DISABLED=ON to set defaults correctly and enable the ASLR check.
      BUILD_SHARED_LIBS can be set to ON (to produce shared libraries, e.g. to save space) or to OFF to
      produce static libraries. If the static libraries will be linked into shared libraries set
      CMAKE_POSITION_INDEPENDENT_CODE=ON, otherwise CMAKE_POSITION_INDEPENDENT_CODE will be set to OFF for maximum efficiency
      of function calls.
* Platform with ASLR enabled --- this is the default. Setting BUILD_SHARED_LIBS=ON in this scenario will produce
      executables that can only be safely used with 1 MPI rank, thus BUILD_SHARED_LIBS will be defaulted to OFF (i.e.
      MADNESS libraries will be build as static libraries). CMAKE_POSITION_INDEPENDENT_CODE is by default set to ON,
      thus MADNESS libraries can be linked into position-independent executables safely. MADNESS libraries can also be
      linked into a shared library and used with more than 1 MPI rank, provided that *ALL* code using MADNESS is part of the *SAME* shared library.
      E.g. to link MADNESS into a Python module compile MADNESS and all libraries using MADNESS as shared libraries
      (with CMAKE_POSITION_INDEPENDENT_CODE=ON) and link them all together into a single module.

## Warning about fast memory allocators

Summary: Only use fast memory allocators if you are using just 1 MPI process or have configured without MPI.

Depending on the calculation and the number of threads being used, MADNESS can receive about a 10% or even more speedup from fast memory allocators such as tcmalloc, jemalloc, tbbmalloc, etc.  However, these **do not work with MPI over InfiniBand** and probably most other transport layers.  It can appear to work, and then fail with either wrong numbers or MPI errors. The reason is that IB requires that memory be pinned and hence MPI introduces its own memory allocator(s) to manage this.  By overriding the allocator, you will break the guarantee that memory is pinned.  

## Toolchain files

**Use of these files is now deprecated --- configuration should usually work without this.**  However, they can be useful if all else fails or on "bleeding-edge" supercomputers with non-standard software environments.

MADNESS provides toolchain files for select systems. They contain the platform specific settings
neccessary to build on the given platform. The toolchain files are included with
the MADNESS source in the cmake/toolchains directory.

* CMAKE_TOOLCHAIN_FILE --- Specifies the path (including the file name) to the
      toolchain file.

For example, to specify the toolchain file for Mira:

    $ cmake -D CMAKE_TOOLCHAIN_FILE=/path/to/madness/source/cmake/toolchains/mira-gcc-essl.cmake \
        /path/to/madness/source
    

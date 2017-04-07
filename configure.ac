AC_PREREQ(2.59)
# was 2.61 but jaguar is not up to date so backported which
# involved replacing AC_PROG_SED and UINT*_T type checks.
AC_INIT([MADNESS], [0.10.1], [https://github.com/m-a-d-n-e-s-s/madness/issues], , [https://github.com/m-a-d-n-e-s-s/madness])
AC_CONFIG_SRCDIR([configure.ac])

# This is where autoconf automatically generated files go
AC_CONFIG_AUX_DIR([config])

# This is where local macros will be stored (by us)
AC_CONFIG_MACRO_DIR([config])

# Default options for automake ??
AM_INIT_AUTOMAKE([-Wall -Werror nostdinc foreign])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# We want #defines output into a header file
AC_CONFIG_HEADERS([src/madness/config.h])

# Detect specific systems so can provide their whacky defaults
acx_default_enable_all_static=no
ACX_CRAYXT
ACX_CRAYXE
ACX_MAC
ACX_IBMBGP
ACX_IBMBGQ
#AC_AIX
ACX_WITH_STUBMPI

# Set default thread processor-affinity
BIND=${BIND-"-1 -1 -1"}
AC_DEFINE_UNQUOTED([MAD_BIND_DEFAULT], ["$BIND"], [The default binding for threads])


# This sets host, host_cpu, host_vendor, host_os
AC_CANONICAL_HOST
AC_MSG_NOTICE([HOST INFORMATION:   $host    $host_cpu    $host_vendor    $host_os])
AC_DEFINE_UNQUOTED([HOST_CPU], ["$host_cpu"], [Defines the host cpu (x86, x86_64, ...).])
AC_DEFINE_UNQUOTED([HOST_SYSTEM], ["$host_os"], [Defines the host os (linux-gnu, ...).])

# Pick out CPUs for which we have special source files (assembler)
use_x86_64_asm=no
use_x86_32_asm=no

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if !(defined(__x86_64__) || defined(_M_X64))
#error Not x86_64
#endif
]])], [use_x86_64_asm=yes])

if test $use_x86_64_asm = no; then
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if !(defined(__i386) || defined(_M_IX86))
#error Not x86
#endif
]])], [use_x86_32_asm=yes])
fi

# Define AM conditional for use of aforementioned assembler
AM_CONDITIONAL([USE_X86_64_ASM], [test $use_x86_64_asm = yes])
AM_CONDITIONAL([USE_X86_32_ASM], [test $use_x86_32_asm = yes])

# All source will be built with the MPI C++ and C compilers
ACX_MPI

initial_CFLAGS="$CFLAGS"
initial_CXXFLAGS="$CXXFLAGS"

AC_PROG_CC([$MPICC])
AC_PROG_CXX([$MPICXX])
AC_PROG_F77([$MPIF77]) # f90 is a superset of f77
AC_LANG_CPLUSPLUS
AC_MSG_NOTICE([Compiler choices after 'AC LANG' ... $CXX and $CC])

CFLAGS="$initial_CFLAGS"
CXXFLAGS="$initial_CXXFLAGS"

# Check that the MPICC Compiler can link correctly.
AC_LANG_PUSH([C])
AC_MSG_CHECKING([that the $MPICC linker works])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[]], [[]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR(["Unable to link with $MPICC])])
AC_LANG_POP

# Check that the MPICXX Compiler can link correctly.
AC_MSG_CHECKING([that the $MPICXX linker works])
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[]], [[]])],
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR(["Unable to link with $MPICXX])])

# Initialize Libtools (for shared libraries)
LT_INIT
LT_LANG([C++])
AC_SUBST(enable_shared)
AC_SUBST(enable_static)
AC_ARG_ENABLE([all-static],
    [AC_HELP_STRING([--enable-all-static@<:@=yes|no@:>@],
      [Link with all static libraries --- likely to fail due to missing system static libraries. Do you really instead want --disable-shared which will avoid building the shared libaries and link against our static libraries?.])],
    [
      case $enableval in
      yes)
        acx_enable_all_static="yes"
      ;;
      no)
        acx_enable_all_static="no"
      ;;
      *)
        AC_MSG_ERROR([Invalid argument '$enableval' given for --enable-all-static])
      esac
    ],
    [acx_enable_all_static="$acx_default_enable_all_static"]
  )

AM_CONDITIONAL([LINK_ALL_STATIC], [test $acx_enable_all_static = yes])

ACX_DETECT_CXX
# Check for C++11 support, and
# add the -std=c++11 or -std=c++0x flag to CXXFLAGS if necessary.
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

# Set debugging
ACX_ENABLE_DEBUGGING

# Set optimization
ACX_ENABLE_OPTIMIZATION

# Set warning flags by compiler vendor
ACX_ENABLE_WARN

# Selelect best compilation flags
ACX_ENABLE_OPTIMAL

# Get optional external libraries inplace so that building will partially check them
#ACX_WITH_LIBUNWIND ... no longer needed for google perf?
ACX_WITH_GOOGLE_PERF
ACX_WITH_GOOGLE_TEST
ACX_WITH_LIBXC
ACX_WITH_PCM
ACX_WITH_TBB
ACX_ENABLE_GENTENSOR

# Set task profiler option
ACX_ENABLE_TASK_PROFILER

# Check for basic build functionality
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AM_PROG_AS
#AC_PROG_SED
AC_CHECK_PROG([HAVE_SED],[sed],[YES],[NO])
if test $HAVE_SED = YES; then
   AC_PATH_PROG(SED,sed)
fi

# We need these to build the documentation
AC_CHECK_PROG([HAVE_DOXYGEN],[doxygen],[YES],[NO])
if test $HAVE_DOXYGEN = YES; then
   AC_PATH_PROG(DOXYGEN,doxygen)
fi
AC_CHECK_PROG([HAVE_DOT],[dot],[YES],[NO])
if test $HAVE_DOT = YES; then
   AC_PATH_PROG(DOT,dot)
fi
AC_CHECK_PROG([HAVE_PDFLATEX],[pdflatex],[YES],[NO])
if test $HAVE_PDFLATEX = YES; then
   AC_PATH_PROG(PDFLATEX,pdflatex)
fi
AC_CHECK_PROG([HAVE_LATEX],[latex],[YES],[NO])
if test $HAVE_LATEX = YES; then
   AC_PATH_PROG(LATEX,latex)
fi
AM_CONDITIONAL([CAN_BUILD_DOC], [test $HAVE_DOXYGEN = YES])

AC_HEADER_TIME
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h sys/param.h sys/time.h unistd.h])
AC_CHECK_HEADERS([bits/atomicity.h ext/atomicity.h])
AC_HEADER_STDBOOL
# Check c++11 headers
# MUST include CXXFLAGS into CPPFLAGS to choose proper C++ standard switch
ac_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $CXXFLAGS"
AC_CHECK_HEADERS([array type_traits tuple initializer_list])
CPPFLAGS="$ac_save_CPPFLAGS"

# Types
AC_TYPE_PID_T
AC_C_RESTRICT
ACX_CHECK_TLS
AC_TYPE_SIZE_T
#AC_TYPE_UINT64_T
#AC_TYPE_INT64_T
AC_CHECK_TYPES([uint64_t, int64_t, uint32_t, int32_t, uint16_t, uint8_t, long long])
AC_C_VOLATILE
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])
AC_CHECK_SIZEOF([void *])
AC_SUBST([AUTOCONF_SIZEOF_VOID_P],[$ac_cv_sizeof_void_p])

# Library functions
AC_CHECK_FUNCS([sleep random execv perror gettimeofday memmove memset pow sqrt strchr strdup getenv])
AC_FUNC_FORK
# EFV: when using Intel compiler on OS X malloc is found to be non-GNU-compatible, which breaks gcc's cstdlib
#AC_FUNC_MALLOC
AC_FUNC_ERROR_AT_LINE
ACX_POSIX_MEMALIGN

# Check for Elemental
ACX_WITH_ELEMENTAL

# Set the default fortran integer type
AC_ARG_WITH([fortran-int32],
            [AC_HELP_STRING([--with-fortran-int32], [Set FORTRAN integer size to 32 bit. By default, assume 32-bit integers.])],
            [
              case $withval in
              yes)
                acx_fortran_int_size=4
              ;;
              no)
                acx_fortran_int_size=8
              ;;
              *)
                AC_MSG_ERROR([Unknown value provided for --with-fortran-int32])
              esac
            ],
            [acx_fortran_int_size=4])

case $acx_fortran_int_size in
4)
  AC_DEFINE([MADNESS_FORTRAN_DEFAULT_INTEGER_SIZE],[4],[Default Fortran integer size.])
  AC_SUBST([MADNESS_FORTRAN_DEFAULT_INTEGER_SIZE],[4])
  AC_MSG_NOTICE([Default Fortran integer: integer*4])
;;
8)
  if test x$ac_cv_type_int64_t = xno; then
    AC_MSG_ERROR([The compiler does not support int64_t, which is required to link with fortran libraries that use 8 byte integers. Rerun configure with '--with-fortran-int32' option.])
  fi
  AC_DEFINE([MADNESS_FORTRAN_DEFAULT_INTEGER_SIZE],[8],[Default Fortran integer size.])
  AC_SUBST([MADNESS_FORTRAN_DEFAULT_INTEGER_SIZE],[8])
  AC_MSG_NOTICE([Default Fortran integer: integer*8])
esac

ACX_WITH_MKL

# Checks for Fortran linking conventions and BLAS+LAPACK at the same time
ACX_FORTRAN_SYMBOLS

# Compiler quirks
ACX_UNQUALIFIED_STATIC_DECL
ACX_STD_ABS

# Check for xterm
AC_CHECK_PROG(XTERM,xterm,xterm,no)
if test $XTERM = 'xterm'; then
  AC_DEFINE(HAVE_XTERM,[1],[If set indicates xterm is in path])
fi

# Check for gdb
AC_CHECK_PROG(GDB,gdb,gdb,no)
if test $GDB = 'gdb'; then
  AC_DEFINE(HAVE_GDB,[1],[If set indicates gdb is in path])
fi

# How to handle MADNESS exceptions.
AC_ARG_ENABLE(madex,
              [AC_HELP_STRING([--enable-madex=arg],
                [Specifies behavior of MADNESS assertions: throw (default), abort, assert, disable])],
              [madex=$enableval],
              [madex=throw])
case $madex in
  abort)
    AC_MSG_NOTICE([MADNESS assertions will abort])
    AC_DEFINE([MADNESS_ASSERTIONS_ABORT], [1], [Set if MADNESS assertions abort])
  ;;
  assert)
    AC_MSG_NOTICE([MADNESS assertions will assert])
    AC_DEFINE([MADNESS_ASSERTIONS_ASSERT], [1], [Set if MADNESS assertions assert])
  ;;
  disable)
    AC_MSG_NOTICE([MADNESS assertions are disabled])
    AC_DEFINE([MADNESS_ASSERTIONS_DISABLE], [1], [Set if MADNESS assertions disabled])
  ;;
  *)
    AC_MSG_NOTICE([MADNESS assertions will throw])
    AC_DEFINE([MADNESS_ASSERTIONS_THROW], [1], [Set if MADNESS assertions throw])
  ;;
esac

# Determine if we should use the instrumented new/delete
AC_ARG_ENABLE([memstats],
              [AC_HELP_STRING([--enable-memstats], [If yes, gather memory statistics (expensive)])],
              [memstats=$enableval],
              [memstats=no])
if test "$memstats" = "yes"; then
  AC_MSG_NOTICE([MADNESS will gather memory statistics (expensive)])
  AC_DEFINE([WORLD_GATHER_MEM_STATS], [1], [Set if MADNESS gathers memory statistics])
fi

AC_ARG_ENABLE([profile], 
              [AC_HELP_STRING([--enable-profile], [Enables profiling])], 
              [AC_MSG_NOTICE([Enabling profiling]); AC_DEFINE(WORLD_PROFILE_ENABLE, [1], [Define if should enable profiling])], 
              [])

AC_ARG_ENABLE([tensor-bounds-checking], 
              [AC_HELP_STRING([--enable-tensor-bounds-checking],
                [Enable checking of bounds in tensors ... slow but useful for debugging])], 
              [AC_MSG_NOTICE([Enabling tensor bounds checking]); AC_DEFINE(TENSOR_BOUNDS_CHECKING, [1], [Define if should enable bounds checking in tensors])], 
              [])

AC_ARG_ENABLE([tensor-instance-count], 
              [AC_HELP_STRING([--enable-tensor-instance-count],
                [Enable counting of allocated tensors for memory leak detection])], 
              [AC_MSG_NOTICE([Enabling tensor instance counting]); AC_DEFINE(TENSOR_INSTANCE_COUNT, [1], [Define if should enable instance counting in tensors])], 
              [])

AC_ARG_ENABLE([spinlocks], 
              [AC_HELP_STRING([--enable-spinlocks],
                [Enables use of spinlocks instead of mutexs (faster unless over subscribing processors)])], 
              [AC_MSG_NOTICE([Enabling use of spinlocks]); AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])], 
              [])

AC_ARG_ENABLE([never-spin], 
              [AC_HELP_STRING([--enable-never-spin],
                [Disables use of spinlocks (notably for use inside virtual machines)])], 
              [AC_MSG_NOTICE([Disabling use of spinlocks]); AC_DEFINE(NEVER_SPIN, [1], [Define if should use never use spinlocks])], 
              [])

AC_ARG_WITH([papi], 
            [AC_HELP_STRING([--with-papi], [Enables use of PAPI])], 
            [AC_MSG_NOTICE([Enabling use of PAPI]); AC_DEFINE(HAVE_PAPI,[1], [Define if have PAPI])], 
            [])

AC_ARG_WITH([parsec],
            [AC_HELP_STRING([--with-parsec], [Build with support for PaRSEC])])
AS_IF([test "x$with_parsec" != "xno"],
      [AS_IF([test "x$with_parsec" != "xyes"],
             [PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:$with_parsec/lib/pkgconfig"])
       PKG_CHECK_MODULES([PaRSEC], [parsec],
                         [AC_DEFINE([HAVE_PARSEC], [1], [Define if have PaRSEC])
                          # Remove the trailing -I from the PaRSEC pkgconfig file
                          PaRSEC_CFLAGS=`echo ${PaRSEC_CFLAGS} | sed 's#-I$##'`
                          CFLAGS="$CFLAGS ${PaRSEC_CFLAGS}"
                          # Don't propagate the -std=c99 C flag to C++
                          PARSEC_CXXFLAGS=`echo ${PaRSEC_CFLAGS} | sed 's#-std=c99##'`
                          CXXFLAGS="$CXXFLAGS ${PARSEC_CXXFLAGS}"
                          LIBADD="$LIBADD ${PaRSEC_LIBS}"],
                         [with_parsec="no"])
       ])
AM_CONDITIONAL([USE_PARSEC], [test "x$with_parsec" != "xno"])

# Capture configuration info for user by compiled code
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_USER],["$USER"],[User that configured the code])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_DATE],["`date`"],[Date of configuration])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_CXX], ["$CXX"],[Configured C++ compiler])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_CXXFLAGS], ["$CXXFLAGS"],[Configured C++ compiler flags])
AC_DEFINE_UNQUOTED([MADNESS_CONFIGURATION_HOST], ["$HOSTNAME"],[Configured on this machine])

# Pass additional configuration information to be used in source 
# ALSO: undefine the standard definitions that every autoheader-using package will define and thus clash
AC_DEFINE_UNQUOTED([MADNESS_PACKAGE_NAME],["$PACKAGE_NAME"],[The full name of this package])
AC_DEFINE_UNQUOTED([MADNESS_PACKAGE_VERSION],["$PACKAGE_VERSION"],[The version of this package])
AC_DEFINE_UNQUOTED([MADNESS_PACKAGE_URL],["$PACKAGE_URL"],[The URL of the home page of this package])
AC_DEFINE_UNQUOTED([MADNESS_PACKAGE_BUGREPORT],["$PACKAGE_BUGREPORT"],[The address where to send bug reports for this package])

# Not yet enabled the user to configure this stuff but probably should
AC_ARG_ENABLE([bsend-ack],
            [AC_HELP_STRING([--disable-bsend-ack], [Use MPI Send instead of MPI Bsend for huge message acknowledgements.])],
            [
              if test "x$enableval" != xno; then
                AC_DEFINE([MADNESS_USE_BSEND_ACKS],[1],[Define when MADNESS should use Bsend for huge message acknowledgements.])
              fi
            ],
            [
             AC_DEFINE([MADNESS_USE_BSEND_ACKS],[1],[Define when MADNESS should use Bsend for huge message acknowledgements.])
            ])

# Set the MPI thread level used by MADNESS
# the default is not MPI_THREAD_MULTIPLE
MADNESS_MPI_THREAD_LEVEL=MPI_THREAD_MULTIPLE
AC_ARG_WITH([mpi-thread],
            [AC_HELP_STRING([--with-mpi-thread], [Thread level for MPI (multiple,serialized)])],
            [
              case $withval in
                multiple)
                  MADNESS_MPI_THREAD_LEVEL=MPI_THREAD_MULTIPLE
                  ;;
                serialized)
                  MADNESS_MPI_THREAD_LEVEL=MPI_THREAD_SERIALIZED
                  ;;
                *)
                  AC_MSG_ERROR([Invalid value for --with-mpi-thread ($withval)])
                  ;;
              esac
            ])
AC_DEFINE_UNQUOTED([MADNESS_MPI_THREAD_LEVEL], [$MADNESS_MPI_THREAD_LEVEL], [Thread-safety level requested from MPI by MADNESS])
AC_MSG_NOTICE([MPI thread-safety requested: $MADNESS_MPI_THREAD_LEVEL])
# to use Elemental MUST use MPI_THREAD_MULTIPLE
if test X$MADNESS_HAS_ELEMENTAL = X1 && test X$MADNESS_MPI_THREAD_LEVEL != XMPI_THREAD_MULTIPLE; then
  AC_MSG_ERROR([Use of Elemental precludes --with-mpi-thread=serialized, must use --with-mpi-thread=multiple])
fi

# Display options
AC_MSG_NOTICE([CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([CFLAGS   = $CFLAGS])
AC_MSG_NOTICE([CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([LDFLAGS  = $LDFLAGS])
AC_MSG_NOTICE([LIBS     = $LIBS])
AC_MSG_NOTICE([build shared libraries: $enable_shared])
AC_MSG_NOTICE([build static libraries: $enable_static])
AC_MSG_NOTICE([link all static libraries: $acx_enable_all_static])

# Set the library suffix used when generating doc/devsamp/Makefile.sample
if test $enable_shared = yes; then

if test $ON_A_MAC = yes; then
# OS X library suffix
AC_SUBST(LIB_SUFFIX,.dylib)
else
# Linux library suffix
AC_SUBST(LIB_SUFFIX,.so)
fi

else

# Set the static library file suffix (same on Linux and OS X)
AC_SUBST(LIB_SUFFIX,.a)

fi


# These are the files that must be generated/preprocessed ... as you add files
# please maintain the directory tree file order and indenting.
AC_CONFIG_FILES([
  Makefile
    config/Makefile
    config/MADNESS.pc
    config/madness-config.cmake
    config/madness-config-version.cmake
    config/madness-project.cmake
    src/Makefile
      src/madness/Makefile
        src/madness/external/Makefile
          src/madness/external/gtest/Makefile
          src/madness/external/muParser/Makefile
          src/madness/external/tinyxml/Makefile
        src/madness/world/Makefile
        src/madness/misc/Makefile
        src/madness/tensor/Makefile
        src/madness/mra/Makefile
      src/apps/Makefile
        src/apps/chem/Makefile
        src/apps/tdse/Makefile
        src/apps/moldft/Makefile
      src/examples/Makefile
        src/examples/compiler/Makefile-prog
    doc/Makefile
    doc/doxygen.cfg
    doc/devsamp/Makefile.sample
])

# The following make files are not generated by default. If you want to build
# these applications, move the appropriate file into the AC_CONFIG_FILES list
# above.

#        src/apps/hf/Makefile
#        src/apps/exciting/Makefile
#        src/apps/polar/Makefile
#        src/apps/nick/Makefile
#        src/apps/interior_bc/Makefile

AC_OUTPUT

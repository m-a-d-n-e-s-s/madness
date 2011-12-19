# This function is used to add performace or system specific compiler flags to
# to CFLAGS and CXXFLAGS environment variables. Users are expected to specify 
# their own optimization flags for unknown compilers, unknown systems, or
# special use cases by adding appropriate values to CFLAGS and CXXFLAGS. It
# should not be used for warning (e.g. -Wall), optimization (e.g. -O3), or debug
# (e.g. -g) flags. Those handled in acx_enable_warn.m4, acx_enable_optimization.m4, 
# and acx_enable_debugging.m4 respectively.
AC_DEFUN([ACX_ENABLE_OPTIMAL], [
  acx_enable_optimal=""
  acx_enable_optimal_save_cxxflags="$CXXFLAGS"
  acx_enable_optimal_flags=""
  acx_enable_optimal_compiler="$CXXVENDOR"
  
  # Allow the user to enable or disable optimal flag
  AC_ARG_ENABLE([optimal],
    [AC_HELP_STRING([--enable-optimal@<:@=yes|no|GNU|clang|Pathscale|Portland|Intel|IBM@:>@],
      [Auto detect optimal CXXFLAGS for compiler and known systems.@<:@default=yes@:>@])],
    [
      case $enableval in
      yes)
        acx_enable_optimal="yes"
      ;;
      no)
        acx_enable_optimal="no"
      ;;
      *)
        acx_enable_optimal="yes"
        acx_enable_optimal_compiler="$enableval"
      esac
    ],
    [acx_enable_optimal="yes"]
  )

  # Set the flags for the specific compilers and systems
  if test $acx_enable_optimal != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    case $acx_enable_optimal_compiler in
      GNU)
        # Delete trailing -stuff from X.X.X-stuff then parse
        CXXVERSION=[`$CXX -dumpversion | sed -e 's/-.*//'`]
        CXXMAJOR=[`echo $CXXVERSION | sed -e 's/\.[.0-9a-zA-Z\-_]*//'`]
        CXXMINOR=[`echo $CXXVERSION | sed -e 's/[0-9]*\.//' -e 's/\.[0-9]*//'`]
        CXXMICRO=[`echo $CXXVERSION | sed -e 's/[0-9]*\.[0-9]*\.//'`]
        AC_MSG_NOTICE([Setting compiler flags for GNU C++ major=$CXXMAJOR minor=$CXXMINOR micro=$CXXMICRO])

        # Flags for all GCC variants
        acx_enable_optimal_flags="$acx_enable_optimal_flags -ffast-math"
        if test $enable_cpp0x = "yes"; then
          acx_enable_optimal_flags="$acx_enable_optimal_flags -std=c++0x"
        fi

        # Add GCC system specific flags
        if test "x$HAVE_CRAYXT" = xyes; then
          ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [-march=barcelona],
            [acx_enable_optimal_flags="$acx_enable_optimal_flags -march=barcelona"])
        elif  test "x$HAVE_IBMBGP" = xyes; then
          acx_enable_optimal_flags=""
        else 
          ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [-march=native],
            [acx_enable_optimal_flags="$acx_enable_optimal_flags -march=native"])
        fi

        # Add flags for Intel x86 architectures. 
        case $host_cpu in
          ??86*)
            acx_enable_optimal_flags="$acx_enable_optimal_flags -mfpmath=sse -msse -mpc64"
           ;;
        esac
      ;;
      
      clang)
        acx_enable_optimal_flags="$acx_enable_optimal_flags"
      ;;

      Pathscale)
        acx_enable_optimal_flags="$acx_enable_optimal_flags"
        if test "x$HAVE_CRAYXT" = xyes; then
          ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [-march=barcelona],
            [acx_enable_optimal_flags="$acx_enable_optimal_flags -march=barcelona"])
        fi
      ;;

      Portland)
        acx_enable_optimal_flags="$acx_enable_optimal_flags -fastsse -Mflushz -Mcache_align -Drestrict=__restrict" 
        AC_MSG_NOTICE([Appending -pgf90libs to LIBS so can link against Fortran BLAS/linalg])
        LIBS="$LIBS -pgf90libs"
        if test "x$HAVE_CRAYXT" = xyes; then
          ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [-tp barcelona-64],
            [acx_enable_optimal_flags="$acx_enable_optimal_flags -tp barcelona-64"])
        fi
      ;;

      Intel)
        acx_enable_optimal_flags="$acx_enable_optimal_flags -ip -no-prec-div -mkl -ansi"
        if test $enable_cpp0x = "yes"; then
          acx_enable_optimal_flags="$acx_enable_optimal_flags -std=c++0x"
        fi
#-use-intel-optimized-headers -fp-model fast=2 -inline-level=2
      ;;
      
      IBM)
        acx_enable_optimal_flags="$acx_enable_optimal_flags"
        if test "x$HAVE_IBMBGP" = xyes; then
          ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [ -qtune=450 -qarch=450d -qlanglvl=extended],
            [acx_enable_optimal_flags="$acx_enable_optimal_flags -qtune=450 -qarch=450d -qlanglvl=extended "])
        fi
      ;;
      
      *)
        AC_MSG_WARN([Optimal flags not set for $acx_enable_optimal_compile compiler])
      ;;
    esac

    # Test the flags and add them to flag variables if successful.
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$acx_enable_optimal_flags],
      [CXXFLAGS="$CXXFLAGS $acx_enable_optimal_flags"],
      [AC_MSG_WARN([$CXX does not accept $acx_enable_optimal_flags, no optimal flags will be used.])])
  fi

])

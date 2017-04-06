AC_DEFUN([ACX_WITH_MKL], [

  AC_ARG_VAR([MKLROOT], [Intel MKL install directory.])

  AC_ARG_WITH([mkl],
    [AS_HELP_STRING([--with-mkl@<:@=yes|no|check|DIR@:>@],
      [Enables use of Intel MKL @<:@default=check@:>@])],
    [
      case $withval in
      yes)
        acx_with_mkl="yes"
      ;;
      no)
        acx_with_mkl="no"
      ;;
      *)
        acx_with_mkl="$withval"
      ;;
      esac
    ],
    [
      acx_with_mkl="check"
    ]
  )
  
  if test $acx_with_mkl != "no"; then
    AC_MSG_CHECKING([Intel MKL])

    case $acx_with_mkl in
    yes|check)
      if test -d "$MKLROOT"; then
        acx_mkl_root="$MKLROOT"
      elif test -d "/opt/intel/mkl"; then
        acx_mkl_root="/opt/intel/mkl"
      else
        acx_mkl_root="no"
      fi
    ;;
    *)
      if test -d "$acx_with_mkl"; then
        acx_mkl_root="$acx_with_mkl"
      else
        acx_mkl_root="no"
      fi
    ;;
    esac
    
    # Show results
    AC_MSG_RESULT([$acx_mkl_root])
    
    if test $acx_mkl_root != "no"; then
      # Set the MKL library name based on the integer type 
      if test $acx_fortran_int_size == "4"; then
        acx_mkl_int_type="lp64"
      else
        acx_mkl_int_type="ilp64"
      fi
      
      AC_DEFINE([HAVE_INTEL_MKL], [1], [Define if using Intel MKL math library])                         

      # Add MKL libraries to the LIBS variable
      if test $ON_A_MAC = "no"; then
        LIBS="$LIBS -L$acx_mkl_root/lib/intel64 -Wl,--start-group -lmkl_intel_$acx_mkl_int_type -lmkl_core -lmkl_sequential -Wl,--end-group -lpthread -lm -ldl"
      else
        LIBS="$LIBS -L$acx_mkl_root/lib -lmkl_intel_$acx_mkl_int_type -lmkl_core -lmkl_sequential -lpthread -lm"
      fi
    else
      # Return an error if not checking
      if test $acx_with_mkl != "check"; then
        AC_MSG_ERROR([Intel MKL install directory not found or invalid.])
      fi
    fi
  fi

])

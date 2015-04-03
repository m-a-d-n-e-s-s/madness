AC_DEFUN([ACX_FORTRAN_SYMBOLS], [
# Dubiously checks for Fortran linking conventions and BLAS+LAPACK at the same time
# mostly to avoid the need for having a fortran compiler installed

# Check for no underscore first since IBM BG ESSL seems to define dgemm with/without underscore
# but dsyev only without underscore ... but AMD ACML also defines both but with different
# interfaces (fortran and c) ... ugh.  Hardwire linking for bgp and restore to original order.

       AC_MSG_NOTICE([Checking Fortran-C linking conventions (dgemm -> ?)])
       fsym=no

       if test $host = "powerpc-bgp-linux-gnu"; then
          fsym="lc"
          AC_MSG_NOTICE([Fortran linking convention is $fsym])
          AC_DEFINE([FORTRAN_LINKAGE_LC],[1],[Fortran-C linking convention lower case (no underscore)])
       fi
       if test $host = "powerpc64-bgq-linux-gnu"; then
          fsym="lc"
          AC_MSG_NOTICE([Fortran linking convention is $fsym])
          AC_DEFINE([FORTRAN_LINKAGE_LC],[1],[Fortran-C linking convention lower case (no underscore)])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm_],[fsym="lcu"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm],[fsym="lc"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([dgemm__],[fsym="lcuu"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([DGEMM],[fsym="uc"])
       fi
       if test $fsym = no; then
           AC_CHECK_FUNC([DGEMM_],[fsym="ucu"])
       fi

# Well there is nothing in the existing path that gives us what we are
# looking for so try looking for some standard examples so that common
# Linux, Apple and configurations work automatically.  We save the
# BLAS library name instead of adding it directly to LIBS since it
# will need to be appened after any LAPACK library that is yet to
# be found.

# OS X
    if test $fsym$ON_A_MAC = noyes; then
        LDFLAGS="$LDFLAGS -framework Accelerate"
        fsym="lcu"
        AC_MSG_NOTICE([Using Accelerate framework for BLAS support])
    fi

# Linux
    BLASLIB=""
    if test $fsym = no; then
        AC_LANG_SAVE
        AC_LANG([C++])
        for blaslib in openblas blas; do
            AC_CHECK_LIB([$blaslib], 
                         [dgemm_], 
                         [fsym="lcu"; BLASLIB="-l$blaslib"; AC_MSG_NOTICE([Found dgemm_ in $blaslib]); break], 
                         [AC_MSG_NOTICE([Unable to find dgemm_ in $blaslib])],
                         [-lpthread])
        done
        AC_LANG_RESTORE
    fi
 
# others ... insert here or extend above for loop if correct symbol is dgemm_
    if test $fsym = no; then
        AC_MSG_ERROR([Could not find dgemm with any known linking conventions])
    fi

    AC_MSG_NOTICE([Fortran linking convention is $fsym]) 

# Now verify that we have at least one of the required lapack routines and again attempt to search for candidate libraries if nothing is found

       if test $fsym = lc; then
           AC_DEFINE([FORTRAN_LINKAGE_LC],[1],[Fortran-C linking convention lower case (no underscore)])
           lapacksym=dsyev
       fi
       if test $fsym = lcu; then
           AC_DEFINE([FORTRAN_LINKAGE_LCU],[1],[Fortran-C linking convention lower case with single underscore])
           lapacksym=dsyev_
       fi
       if test $fsym = lcuu; then
           AC_DEFINE([FORTRAN_LINKAGE_LCUU],[1],[Fortran-C linking convention lower case with double underscore])
           lapacksym=dsyev__
       fi
       if test $fsym = uc; then
           AC_DEFINE([FORTRAN_LINKAGE_UC],[1],[Fortran-C linking convention upper case])
           lapacksym=DSYEV
       fi
       if test $fsym = ucu; then
           AC_DEFINE([FORTRAN_LINKAGE_UCU],[1],[Fortran-C linking convention upper case with single underscore])
           lapacksym=DSYEV_
       fi

       LAPACKLIB=""       
       foundlapack=no
       AC_CHECK_FUNC([$lapacksym],[foundlapack=yes],AC_MSG_NOTICE([Could not find dsyev with selected linking convention in default library path]))

       LAPACKLIB=""       
       if test $foundlapack = no; then
           AC_LANG_SAVE
           AC_LANG([C++])
           for lapacklib in lapack; do
               AC_CHECK_LIB([$lapacklib], 
                            [$lapacksym], 
                            [foundlapack=yes; LAPACKLIB="-l$lapacklib"; AC_MSG_NOTICE([Found $lapacksym in $lapacklib]); break], 
                            [AC_MSG_NOTICE([Unable to find $lapacksym in $lapackib])],
                            [$BLASLIB -lpthread])
           done
           AC_LANG_RESTORE
       fi

       if test $foundlapack = no; then
           AC_MSG_NOTICE([Could not find $lapacksym in any known library ... specify LAPACK library via LIBS])
       fi

       LIBS="$LIBS $LAPACKLIB $BLASLIB"
])






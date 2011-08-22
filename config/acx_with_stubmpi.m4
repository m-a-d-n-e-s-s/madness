AC_DEFUN([ACX_WITH_STUBMPI],
[
  acx_with_stubmpi=no
  AC_ARG_WITH([stubmpi],
    [AS_HELP_STRING([--with-stubmpi], [Build without MPI ... i.e., stubbing it out.])],
    [
      case $withval in
      yes)
        acx_with_stubmpi=yes
      ;;
      no)
        acx_with_stubmpi=no
      ;;
      *)
        acx_with_stubmpi=yes
      ;;
      esac
    ])
    if test $acx_with_stubmpi = yes; then
        AC_DEFINE(STUBOUTMPI,[1],[If defined header disable MPI by including stubmpi.h])
        MPICC="$CC"
	MPICXX="$CXX"
        if test "x$CC" = x; then
           MPICC=gcc;
        fi
        if test "x$CXX" = x; then
           MPICXX=g++;
        fi
        echo "Stubbing out MPI with MPICXX=$MPICXX MPICC=$MPICC"
    fi
])

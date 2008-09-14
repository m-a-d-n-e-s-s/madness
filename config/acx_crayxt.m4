AC_DEFUN([ACX_CRAYXT], [
        # If on a Cray XT 
        #   - defines HAVE_CRAYXT=1 in headers 
        #   - defines HAVE_CRAYXT=yes in the script
        #   - sets MPICXX=CC and MPICC=cc if the user has not already set them
        AC_CHECK_FILE([/proc/cray_xt],[HAVE_CRAYXT=yes AC_DEFINE([HAVE_CRAYXT],[1],[Defined if we are running on a Cray-XT])],[])
        if test "x$HAVE_CRAYXT" = xyes; then
                AC_DEFINE(AMD_QUADCORE_TUNE,[1],"Target for tuning mtxmq kernels")
                if test "x$MPICC" = x; then
                        echo "Choosing MPICC=cc for Cray-XT"
                        MPICC=cc;
                fi
                if test "x$MPICXX" = x; then
                        echo "Choosing MPICXX=CC for Cray-XT"
                        MPICXX=CC;
                fi
        fi
])





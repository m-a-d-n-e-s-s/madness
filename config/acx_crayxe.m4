AC_DEFUN([ACX_CRAYXE], [
        # If on a Cray XE 
        #   - defines HAVE_CRAYXE=1 in headers 
        #   - defines HAVE_CRAYXE=yes in the script
        #   - sets MPICXX=CC and MPICC=cc if the user has not already set them
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        AC_TRY_COMPILE(,[
          #ifdef __CRAYXE
          int ok;
          #else
          choke me
          #endif
        ],[HAVE_CRAYXE=yes AC_DEFINE([HAVE_CRAYXE],[1],[Defined if we are running on an Cray XE])],[echo This is not a Cray XE.])
        if test "x$HAVE_CRAYXE" = xyes; then
                AC_DEFINE(AMD_QUADCORE_TUNE,[1],"Target for tuning mtxmq kernels")
                if test "x$MPICC" = x; then
                        echo "Choosing MPICC=cc for Cray XE"
                        MPICC=cc;
                fi
                if test "x$MPICXX" = x; then
                        echo "Choosing MPICXX=CC for Cray XE"
                        MPICXX=CC;
                fi
                echo "int main(){return 0;}" > dddd.cc
                CC dddd.cc -lacml
                if test $? = 0; then
                        echo "AMD ACML library detected"
                        LIBS="$LIBS -lacml"
                        AC_DEFINE(HAVE_ACML,[1],[Define if AMD math library available - ACML])
                fi
                BIND="1 0 2"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])





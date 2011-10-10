AC_DEFUN([ACX_IBMBGQ],[
        # If on an IBMBGQ
        #   - defines HAVE_IBMBGQ=1 in headers 
        #   - defines HAVE_IBMBGQ=yes in the script
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        #   - sets the host architecture to powerpc-bgq-linux-gnu
        #
        AC_TRY_COMPILE(,[
          #ifdef __bgq__
          int ok;
          #else
          choke me
          #endif
        ],[HAVE_IBMBGQ=yes AC_DEFINE([HAVE_IBMBGQ],[1],[Defined if we are running on an IBM Blue Gene/Q])],[echo This is not an IBM Blue Gene/Q.])
        if test "x$HAVE_IBMBGQ" = xyes; then
                echo "IBM Blue Gene/Q system detected"
                host="powerpc64-bgq-linux"
                host_triplet="powerpc64-bgq-linux"

                BIND="-1 -1 -1"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])


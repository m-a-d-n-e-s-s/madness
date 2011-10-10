AC_DEFUN([ACX_IBMBGP],[
        # If on an IBMBGP
        #   - defines HAVE_IBMBGP=1 in headers 
        #   - defines HAVE_IBMBGP=yes in the script
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        #   - sets the host architecture to powerpc-bgp-linux-gnu
        #
        #AC_CHECK_FILE([/bgsys], [HAVE_IBMBGP=yes AC_DEFINE([HAVE_IBMBGP],[1],[Defined if we are running on an IBM Blue Gene/P])],[])
        AC_TRY_COMPILE(,[
          #ifdef __bgp__
          int ok;
          #else
          choke me
          #endif
        ],[HAVE_IBMBGP=yes AC_DEFINE([HAVE_IBMBGP],[1],[Defined if we are running on an IBM Blue Gene/P])],[echo This is not an IBM Blue Gene/P.])
        if test "x$HAVE_IBMBGP" = xyes; then
                echo "IBM Blue Gene/P system detected"
                host="powerpc-bgp-linux"
                host_triplet="powerpc-bgp-linux"

                BIND="-1 -1 -1"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])


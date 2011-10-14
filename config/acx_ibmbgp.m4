AC_DEFUN([ACX_IBMBGP],[
        # If on an IBMBGP
        #   - defines HAVE_IBMBGP=1 in headers 
        #   - defines HAVE_IBMBGP=yes in the script
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        #   - sets the host architecture to powerpc-bgp-linux-gnu
        #
        #AC_CHECK_FILE([/bgsys], [HAVE_IBMBGP=yes AC_DEFINE([HAVE_IBMBGP],[1],[Defined if we are running on an IBM Blue Gene/P])],[])
        echo "int main(){ #ifdef __bgp__ \ return 0; #else \ choke me #endif }" > __bgp__.cc
        mpicxx __bgp__.cc
        if test $? = 0; then
                echo "IBM Blue Gene/P detected"
                HAVE_IBMBGP=yes
                AC_DEFINE(HAVE_IBMBGP,[1],[Defined if we are running on an IBM Blue Gene/P])
        fi
        if test "x$HAVE_IBMBGP" = xyes; then
                host="powerpc-bgp-linux"
                host_triplet="powerpc-bgp-linux"

                BIND="-1 -1 -1"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])


AC_DEFUN([ACX_IBMBGQ],[
        # If on an IBMBGQ
        #   - defines HAVE_IBMBGQ=1 in headers 
        #   - defines HAVE_IBMBGQ=yes in the script
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        #   - sets the host architecture to powerpc-bgq-linux-gnu
        #
        #AC_CHECK_FILE([/bgsys], [HAVE_IBMBGQ=yes AC_DEFINE([HAVE_IBMBGQ],[1],[Defined if we are running on an IBM Blue Gene/Q])],[])
        echo "int main(){ #ifdef __bgq__ \ return 0; #else \ choke me #endif }" > __bgq__.cc
        mpicxx __bgq__.cc
        if test $? = 0; then
                echo "IBM Blue Gene/Q detected"
                HAVE_IBMBGQ=yes
                AC_DEFINE(HAVE_IBMBGQ,[1],[Defined if we are running on an IBM Blue Gene/Q])
        fi
        if test "x$HAVE_IBMBGQ" = xyes; then
                host="powerpc-bgq-linux"
                host_triplet="powerpc-bgq-linux"

                BIND="-1 -1 -1"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])


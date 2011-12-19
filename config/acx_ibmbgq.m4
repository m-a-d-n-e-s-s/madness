AC_DEFUN([ACX_IBMBGQ],[
        # If on an IBMBGQ
        #   - defines HAVE_IBMBGQ=1 in headers 
        #   - defines HAVE_IBMBGQ=yes in the script
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        #   - sets the host architecture to powerpc-bgq-linux-gnu
        #
        #AC_CHECK_FILE([/bgsys], [HAVE_IBMBGQ=yes AC_DEFINE([HAVE_IBMBGQ],[1],[Defined if we are running on an IBM Blue Gene/Q])],[])
        echo "int main()"       > __bgq__.cc
        echo "{"               >> __bgq__.cc
        echo "#ifdef __bgq__"  >> __bgq__.cc
        echo "return 0;"       >> __bgq__.cc
        echo "#else"           >> __bgq__.cc
        echo "choke me"        >> __bgq__.cc
        echo "#endif"          >> __bgq__.cc
        echo "}"               >> __bgq__.cc
        mpicxx __bgq__.cc >& /dev/null
        if test $? = 0; then
                echo "IBM Blue Gene/Q detected"
                HAVE_IBMBGQ=yes
                AC_DEFINE(HAVE_IBMBGQ,[1],[Defined if we are running on an IBM Blue Gene/Q])
        fi
        /bin/rm __bgq__.cc
        if test "x$HAVE_IBMBGQ" = xyes; then
                host="powerpc64-bgq-linux"
                host_triplet="powerpc64-bgq-linux"

                BIND="-1 -1 -1"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])


AC_DEFUN([ACX_CRAYXT], [
        # If on a Cray XT 
        #   - defines HAVE_CRAYXT=1 in headers 
        #   - defines HAVE_CRAYXT=yes in the script
        #   - sets MPICXX=CC and MPICC=cc if the user has not already set them
        #   - sets thread binding to "1 0 2"
        #   - enables spinlocks
        echo "int main()"       > __crayxt.cc
        echo "{"               >> __crayxt.cc
        echo "#ifdef __CRAYXT" >> __crayxt.cc
        echo "return 0;"       >> __crayxt.cc
        echo "#else"           >> __crayxt.cc
        echo "choke me"        >> __crayxt.cc
        echo "#endif"          >> __crayxt.cc
        echo "}"               >> __crayxt.cc
        CC __crayxt.cc >& /dev/null
        if test $? = 0; then
                echo "Cray XT detected"
                HAVE_CRAYXT=yes
                AC_DEFINE(HAVE_CRAYXT,[1],[Defined if we are running on an Cray XT])
        fi
        /bin/rm __crayxt.cc
        if test "x$HAVE_CRAYXT" = xyes; then
                AC_DEFINE(AMD_QUADCORE_TUNE,[1],"Target for tuning mtxmq kernels")
                if test "x$MPICC" = x; then
                        echo "Choosing MPICC=cc for Cray XT"
                        MPICC=cc;
                fi
                if test "x$MPICXX" = x; then
                        echo "Choosing MPICXX=CC for Cray XT"
                        MPICXX=CC;
                fi
                echo "int main(){return 0;}" > __acml.cc
                CC __acml.cc -lacml >& /dev/null
                if test $? = 0; then
                        echo "AMD ACML library detected"
                        LIBS="$LIBS -lacml"
                        AC_DEFINE(HAVE_ACML,[1],[Define if AMD math library available - ACML])
                fi
                /bin/rm __acml.cc
                BIND="1 0 2"
                AC_DEFINE(USE_SPINLOCKS, [1], [Define if should use spinlocks])
        fi
])





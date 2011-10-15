AC_DEFUN([ACX_CRAYXE], [
        # If on a Cray XE 
        #   - defines HAVE_CRAYXE=1 in headers 
        #   - defines HAVE_CRAYXE=yes in the script
        #   - sets MPICXX=CC and MPICC=cc if the user has not already set them
        #   - sets thread binding to "1 0 2" TODO: this has to be wrong on AMD Magny Cours
        #   - enables spinlocks
        echo "int main()"       > __crayxe.cc
        echo "{"               >> __crayxe.cc
        echo "#ifdef __CRAYXE" >> __crayxe.cc
        echo "return 0;"       >> __crayxe.cc
        echo "#else"           >> __crayxe.cc
        echo "choke me"        >> __crayxe.cc
        echo "#endif"          >> __crayxe.cc
        echo "}"               >> __crayxe.cc
        CC __crayxe.cc >& /dev/null
        if test $? = 0; then
                echo "Cray XE detected"
                HAVE_CRAYXE=yes
                AC_DEFINE(HAVE_CRAYXE,[1],[Defined if we are running on an Cray XE])
        fi
        /bin/rm __crayxe.cc
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





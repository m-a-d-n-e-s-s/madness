AC_DEFUN([AC_DETECT_CXX], [
    # Sets environment variable CXXVENDOR to one of
    # [GNU,Intel,Portland,Pathscale,unknown]
    CXXVENDOR=unknown
    if test $CXXVENDOR = unknown; then
        $CXX --version 2>&1 | egrep -q "GCC|GNU|gcc|gnu|g\+\+"
        if test $? = 0; then
           echo GNU C++ compiler detected by MADNESS
           CXXVENDOR=GNU
        fi
    fi
    if test $CXXVENDOR = unknown; then
        $CXX --version 2>&1 | grep -q "Intel"
        if test $? = 0; then
           echo Intel C++ compiler detected by MADNESS
           CXXVENDOR=Intel
        fi
    fi
    if test $CXXVENDOR = unknown; then
        $CXX --version 2>&1 | grep -q "Portland"
        if test $? = 0; then
           echo 'Portland C++ compiler detected by MADNESS'
           #echo 'Working around Portland std::exit bug'
           #echo '#ifdef __cplusplus' >> confdefs.h
           #echo '#include <cstdlib>' >> confdefs.h
           #echo 'namespace std {void exit(int){::exit(int);}}' >> confdefs.h
           CXXVENDOR=Portland
        fi
    fi
    if test $CXXVENDOR = unknown; then
        $CXX -v 2>&1 | grep -q "Pathscale"
        if test $? = 0; then
           echo Pathscale C++ compiler detected by MADNESS
           CXXVENDOR=Pathscale
        fi
    fi
    if test $CXXVENDOR = unknown; then
        echo Unknown vendor for C++ compiler
    fi
])





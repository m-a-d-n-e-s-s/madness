AC_DEFUN([ACX_ENABLE_CPP0XTR1],
[
        # default is to use compiler support for either/both C++0x or TR1 headers
        # this is indicated by $enable_cpp0x=yes and $enable_cpptr1=yes
        #
        # other values mean disabled
        AC_ARG_ENABLE([cpp0x],  AS_HELP_STRING([--disable-cpp0x],  [Disable use of compiler support for C++0x (default is enabled or yes)]))
        AS_IF([test "x$enable_cpp0x" != "xno"], [enable_cpp0x="yes"])

        AC_ARG_ENABLE([cpptr1], AS_HELP_STRING([--disable-cpptr1], [Disable use of compiler provided TR1 headers (default is enabled or yes)]))
        AS_IF([test "x$enable_cpptr1" != "xno"], [enable_cpptr1="yes"])

        echo "enable_cpp0x=$enable_cpp0x  enable_cpptr1=$enable_cpptr1" 
])

        

AC_DEFUN([ACX_WITH_BOOST],
[
  AC_ARG_WITH([boost],
            [AS_HELP_STRING([--with-boost@<:@=DIR@:>@], [DIR where the Boost library is installed.])],
            [
              case $with_boost in
                yes)
                  AC_MSG_ERROR([You must specify the install path for Boost.])
                ;;
                no)
                ;;
                *)
                  CPPFLAGS="-I$with_boost/include/boost/tr1/tr1 -I$with_boost/include $CPPFLAGS"
                  AC_DEFINE([MADNESS_USE_BOOST_TR1], [1], 
                    [Madness will use Boost.TR1 where the comiler lacks support for TR1.])
                    
                  # Check for the pressence of the Boost TR1 header files.
                  AC_CHECK_HEADER([boost/tr1/tr1/memory], [], [AC_MSG_ERROR([Unable to find the Boost TR1 memory header file.])])
                ;;
              esac
            ])
])

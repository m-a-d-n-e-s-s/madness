AC_DEFUN([ACX_WITH_BOOST],
[
  acx_with_boost=""
  AC_ARG_WITH([boost],
    [AS_HELP_STRING([--with-boost@<:@=Install DIR@:>@], [Build with Boost TR1 library.])],
    [
      case $withval in
      yes)
        AC_MSG_ERROR([You must specify the install path for Boost.])
      ;;
      no)
        acx_with_boost="no"
      ;;
      *)
        acx_with_boost=$withval
        CPPFLAGS="-I$withval/include/boost/tr1/tr1 -I$withval/include -I$withval/boost/tr1/tr1 -I$withval $CPPFLAGS"
        AC_DEFINE([MADNESS_HAS_BOOST_TR1], [1], 
          [Madness will use Boost.TR1 where the compiler lacks support for TR1.])
                    
        # Check for the pressence of the Boost TR1 header files.
        AC_CHECK_HEADER([boost/tr1/tr1/memory], [],
          [AC_MSG_ERROR([Unable to find the Boost TR1 memory header file.])])
        AC_CHECK_HEADER([boost/make_shared.hpp],
          [AC_DEFINE([MADNESS_HAS_BOOST_MAKE_SHARED], [1], [Madness has Boost make_shared and allocate_shared available.])],
          [AC_MSG_ERROR([Unable to find the Boost make_shared / allocate header file.])])
        AC_CHECK_HEADER([boost/tr1/tr1/type_traits], [],
          [AC_MSG_ERROR([Unable to find the Boost TR1 type_traits header file.])])
        AC_CHECK_HEADER([boost/tr1/tr1/array], [],
          [AC_MSG_ERROR([Unable to find the Boost TR1 array header file.])])
    ;;
  esac
  ], [acx_with_boost="no"])
])

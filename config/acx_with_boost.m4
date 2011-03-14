AC_DEFUN([ACX_WITH_BOOST],
[
  acx_with_boost=no
  AC_ARG_WITH([boost],
    [AS_HELP_STRING([--with-boost@<:@=Install DIR@:>@], [Build with Boost TR1 library.])],
    [
      case $withval in
      yes)
        acx_with_boost=yes
      ;;
      no)
        acx_with_boost=no
      ;;
      *)
        acx_with_boost=yes
        CPPFLAGS="-I$withval/include $CPPFLAGS"
      ;;
      esac
    ],
    [acx_with_boost=no]
  )
  
  if test "$acx_with_boost" != no; then
    # Check for the pressence of the Boost TR1 header files.
    AC_CHECK_HEADER([boost/tr1/memory.hpp], [],
      [AC_MSG_ERROR([Unable to find the Boost TR1 memory header file.])])
    AC_CHECK_HEADER([boost/make_shared.hpp],
      [AC_DEFINE([MADNESS_HAS_BOOST_MAKE_SHARED], [1], [Madness has Boost make_shared and allocate_shared available.])],
      [AC_MSG_ERROR([Unable to find the Boost make_shared / allocate header file.])])
    AC_CHECK_HEADER([boost/tr1/type_traits.hpp], [],
      [AC_MSG_ERROR([Unable to find the Boost TR1 type_traits header file.])])
    AC_CHECK_HEADER([boost/tr1/array.hpp], [],
      [AC_MSG_ERROR([Unable to find the Boost TR1 array header file.])])

    AC_DEFINE([MADNESS_HAS_BOOST_TR1], [1], 
      [Madness will use Boost.TR1 where the compiler lacks support for TR1.])
  fi
])

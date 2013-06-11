AC_DEFUN([ACX_WITH_ELEMENTAL],
[
  acx_with_elemental=no
  AC_ARG_WITH([elemental],
    [AS_HELP_STRING([--with-elemental@<:@=DIR@:>@], [Build with Elemental headers.])],
    [
      case $withval in
      yes)
        acx_with_elemental=yes
      ;;
      no)
        acx_with_elemental=no
      ;;
      *)
        acx_with_elemental=yes
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        LIBS="$LIBS -L$withval/lib -lelemental -lpmrrr -lelem-dummy-lib $LIBS"
      ;;
      esac
    ]
  )
  
  if test "$acx_with_elemental" != no; then
    # Check for the pressence of Elemental header files.
    AC_CHECK_HEADER([elemental.hpp], [],
      [AC_MSG_ERROR([Unable to find the elemental.hpp  header file.])])
    AC_DEFINE([MADNESS_HAS_ELEMENTAL], [1], 
      [Madness will use Elemental for parallel linear algebra operations])
  fi
])

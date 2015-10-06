#
# Use LIBXSMM (https://github.com/hfp/libxsmm) for mTxm.
# LIBXSMM supports only double (and single, which is not relevant,
# so this does not change anything for the complex case.
#
AC_DEFUN([ACX_WITH_LIBXSMM],
[
  acx_with_libxsmm=no
  AC_ARG_WITH([libxsmm],
    [AS_HELP_STRING([--with-libxsmm@<:@=DIR@:>@], [Build with LIBXSMM.])],
    [
      case $withval in
      yes)
        acx_with_libxsmm=yes
      ;;
      no)
        acx_with_libxsmm=no
      ;;
      *)
        acx_with_libxsmm=yes
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        LIBXSMM_LIBS="-L$withval/lib -lxsmm "
	LIBS="$LIBXSMM_LIBS $LIBS"
      ;;
      esac
    ]
  )
  
  if test "$acx_with_libxsmm" != no; then
    # Check for the presence of libxsmm header files.
    AC_CHECK_HEADER([libxsmm.h], [],
      [AC_MSG_ERROR([Unable to find the libxsmm.h header file.])])
    AC_DEFINE([MADNESS_USE_LIBXSMM], [1], 
      [Madness will use LIBXSMM for real double case of mTxm.])
    MADNESS_HAS_LIBXSMM=1
  fi
])

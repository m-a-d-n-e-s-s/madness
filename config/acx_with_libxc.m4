AC_DEFUN([ACX_WITH_LIBXC], [
  acx_with_libxc=""
  AC_ARG_WITH([libxc],
    [AS_HELP_STRING([--with-libxc@<:@=Install DIR@:>@],
      [Enables use of the libxc library of density functionals])],
    [
      case $withval in
      yes)
        acx_with_libxc="yes"
      ;;
      no)
        acx_with_libxc="no"
      ;;
      *)
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        LIBS="$LIBS -L$withval/lib"
        acx_with_libxc="$withval"
      esac
    ],
    [acx_with_libxc="no"]
  )
  if test $acx_with_libxc != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    AC_CHECK_LIB([xc], [xc_func_end], [LIBS="$LIBS -lxc"], [AC_MSG_ERROR(["Unable to link with libxc])])
    AC_DEFINE([MADNESS_HAS_LIBXC], [1], [Define if using libxc])
    AC_LANG_RESTORE
  fi
])

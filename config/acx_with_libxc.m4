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
    [acx_with_libxc="yes"]
  )

  if test $acx_with_libxc != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    AC_CHECK_HEADERS([xc.h xc_funcs.h], [], 
                     [acx_with_libxc=no
                      AC_MSG_NOTICE([Unable to include with xc.h or xc_func.h])])
    AC_CHECK_LIB([xc], [xc_func_end], [], 
                 [acx_with_libxc=no
                  AC_MSG_NOTICE([Unable to link with libxc])])


    AC_LANG_RESTORE
  fi

  if test $acx_with_libxc != "no"; then
    AC_DEFINE([MADNESS_HAS_LIBXC], [1], [Define if using libxc])                         
  fi

  AM_CONDITIONAL([MADNESS_HAS_LIBXC], [test $acx_with_libxc != "no"])
])

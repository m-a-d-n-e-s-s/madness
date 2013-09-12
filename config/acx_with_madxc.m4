AC_DEFUN([ACX_WITH_MADXC], [
  acx_with_madxc=""
  AC_ARG_WITH([madxc],
    [AS_HELP_STRING([--with-madxc@<:@=Install DIR@:>@],
      [Enables use of the madxc library of density functionals])],
    [
      case $withval in
      yes)
        acx_with_madxc="yes"
      ;;
      no)
        acx_with_madxc="no"
      ;;
      *)
        CPPFLAGS="-I$withval $CPPFLAGS"
        LIBS="$LIBS -L$withval"
        acx_with_madxc="$withval"
      esac
    ],
    [acx_with_madxc="no"]
  )
  if test $acx_with_madxc != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    AC_CHECK_HEADERS([libMADxc.h], [], [AC_MSG_ERROR(["Unable to include with  libMADxc.h])])
    AC_CHECK_LIB([MADxc], [rks_c_vwn5_ ], [], [AC_MSG_ERROR(["Unable to link with libMADxc])])
    AC_DEFINE([MADNESS_HAS_MADXC], [1], [Define if using libMADxc])
    AC_LANG_RESTORE
  fi

  AM_CONDITIONAL([MADNESS_HAS_MADXC], [test $acx_with_madxc != "no"])

])


AC_DEFUN([ACX_WITH_PCM], [
  acx_with_pcm=""
  AC_ARG_WITH([pcm],
    [AS_HELP_STRING([--with-pcm@<:@=Install DIR@:>@],
      [Enables use of the polarizable contiuum model library PCM])],
    [
      case $withval in
      yes)
        acx_with_pcm="yes"
      ;;
      no)
        acx_with_pcm="no"
      ;;
      *)
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        LIBS="$LIBS -L$withval/lib"
        acx_with_pcm="$withval"
      esac
    ],
    [acx_with_pcm="yes"]
  )

  if test $acx_with_pcm != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    AC_CHECK_HEADERS([PCMSolver/pcmsolver.h PCMSolver/PCMInput.h], [], 
                     [acx_with_pcm=no
                      AC_MSG_NOTICE([Unable to include with pcmsolver.h or PCMInput.h])])
    AC_CHECK_LIB([pcm],[pcmsolver_is_compatible_library], [], 
                 [acx_with_pcm=no
                  AC_MSG_NOTICE([Unable to link with pcm])])


    AC_LANG_RESTORE
  fi

  if test $acx_with_pcm != "no"; then
    AC_DEFINE([MADNESS_HAS_PCM], [1], [Define if using pcm])                         
  fi

  AM_CONDITIONAL([MADNESS_HAS_PCM], [test $acx_with_pcm != "no"])
])

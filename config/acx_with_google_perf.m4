AC_DEFUN([ACX_WITH_GOOGLE_PERF], [
  acx_with_google_perf=""
  AC_ARG_WITH([google-perf],
    [AS_HELP_STRING([--without-google-perf],
      [Disables use of Google fast malloc.profiler/heapchecker; --with-google-perf=directory overrides install path])],
    [
      case $withval in
      yes)
        acx_with_google_perf="yes"
      ;;
      no)
        acx_with_google_perf="no"
      ;;
      *)
        LIBS="$LIBS -L$withval/lib"
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        acx_with_google_perf="$withval"
      esac
    ],
    [acx_with_google_perf="yes"]
  )

  if test $acx_with_google_perf != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])

    AC_CHECK_LIB([tcmalloc_and_profiler], [ProfilerStart], 
                 [LIBS="-lprofiler -ltcmalloc $LIBS" 
                  have_profiler=yes
                  AC_MSG_NOTICE([Google profiler and fast malloc found])
                  AC_DEFINE([MADNESS_HAS_GOOGLE_PERF], [1], [Define if using Google PerformanceTools])],
                 [have_profiler=no
                  AC_MSG_NOTICE(["Unable to link with libprofiler+libtcmalloc])])

    if test $have_profiler = "no" ; then
      AC_CHECK_LIB([tcmalloc_minimal], [malloc], 
                   [LIBS="-ltcmalloc_minimal $LIBS"
                    AC_MSG_NOTICE([Google minimal fast malloc found])
                    AC_DEFINE([MADNESS_HAS_GOOGLE_PERF_MINIMAL], [1], [Define if using Google PerformanceTools without libunwind])],
                   [AC_MSG_NOTICE(["Unable to link with libtcmalloc_minimal])])
    fi

    AC_LANG_RESTORE
  fi
])

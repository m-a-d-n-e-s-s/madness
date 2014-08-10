AC_DEFUN([ACX_WITH_GOOGLE_TEST], [
  acx_with_google_test=""
  AC_ARG_WITH([google-test],
    [AS_HELP_STRING([--with-google-test@<:@=yes|no@:>@],
      [Enables use of Google unit test @<:@default=no@:>@])],
    [
      case $withval in
      yes)
        acx_with_google_test="yes"
      ;;
      *)
        acx_with_google_test="no"
      esac
    ],
    [
      acx_with_google_test="no"
    ]
  )

  AC_ARG_VAR([GTEST_CPPFLAGS], [C-like preprocessor flags for Google Test.])
  AC_ARG_VAR([GTEST_CXXFLAGS], [C++ compile flags for Google Test.])
  AC_ARG_VAR([GTEST_LDFLAGS], [Linker path and option flags for Google Test.])
  AC_ARG_VAR([GTEST_LIBS], [Library linking flags for Google Test.])

  if test $acx_with_google_test != "no"; then
    if test $acx_with_boost != "no"; then
      GTEST_CPPFLAGS="$GTEST_CPPFLAGS -DGTEST_HAS_TR1_TUPLE=0 -DGTEST_USE_OWN_TR1_TUPLE=1"
    fi

    # Set preprocessor and build variables
    AC_DEFINE([MADNESS_HAS_GOOGLE_TEST], [1], [Define if should use Google unit testing])
    AC_SUBST([GTEST_CPPFLAGS])
    AC_SUBST([GTEST_CXXFLAGS])
    AC_SUBST([GTEST_LDFLAGS])
    AC_SUBST([GTEST_LIBS])

  fi

  AM_CONDITIONAL([MADNESS_HAS_GOOGLE_TEST], [test $acx_with_google_test != no])
])

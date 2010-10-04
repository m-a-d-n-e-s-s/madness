AC_DEFUN([ACX_WITH_GOOGLE_TEST], [
  acx_with_google_test=""
  AC_ARG_WITH([google-test],
    [AS_HELP_STRING([--with-google-test@<:@=Install DIR@:>@],
      [Enables use of Google unit test])],
    [
      case $withval in
      yes)
        acx_with_google_test="yes"
      ;;
      no)
        acx_with_google_test="no"
      ;;
      *)
        CPPFLAGS="$CPPFLAGS -I$withval/include"
        LDFLAGS="$LDFLAGS -L$withval/lib"
        acx_with_google_test="$withval"
      esac
    ],
    [acx_with_google_test="no"]
  )
  if test $acx_with_google_test != "no"; then
    AC_LANG_SAVE
    AC_LANG([C++])
    if test $acx_with_boost != "no"; then
      AC_DEFINE([GTEST_HAS_TR1_TUPLE], [0], [Define this as 0 when using Boost TR1 and Google Test])
      AC_DEFINE([GTEST_USE_OWN_TR1_TUPLE], [1], [Define this as 1 when using Boost TR1 and Google Test])
    fi
    AC_CHECK_HEADER([gtest/gtest.h], [], [AC_MSG_ERROR([Unable to compile with Google Test.])])
    AC_MSG_CHECKING([for ::testing::InitGoogleTest in -lgtest])
    LIBS="$LIBS -lgtest"
    AC_LINK_IFELSE(
      [AC_LANG_PROGRAM(
        [[namespace testing { void InitGoogleTest(int*,char**); }]],
        [[int i = 0; char** c = 0; ::testing::InitGoogleTest(&i,c);]])
      ],
      [AC_MSG_RESULT([yes])],
      [AC_MSG_RESULT([no])
       AC_MSG_ERROR(["Unable to link with libgtest])])
    AC_DEFINE([MADNESS_HAS_GOOGLE_TEST], [1], [Define if should use Google unit testing])
    AC_LANG_RESTORE
  fi
])
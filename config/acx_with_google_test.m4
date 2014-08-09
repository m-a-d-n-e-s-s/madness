AC_DEFUN([ACX_WITH_GOOGLE_TEST], [
  acx_with_google_test=""
  AC_ARG_WITH([google-test],
    [AS_HELP_STRING([--with-google-test@<:@=Install DIR|yes|no@:>@],
      [yes/no (en)/(dis)ables use of Google unit test; or specify different gtest installation dir to override use of internal version @<:@default=no@:>@])],
    [
      case $withval in
      yes)
        acx_with_google_test="yes"
        GTEST_CPPFLAGS="-I`pwd`/src/lib/gtest/include"
        GTEST_LDFLAGS="-L`pwd`/src/lib/gtest"
        GTEST_LIBS="-lMADgtest"
      ;;
      no)
        acx_with_google_test="no"
      ;;
      *)
        acx_with_google_test="yes"
        GTEST_CPPFLAGS="-I$withval/include"
        GTEST_LDFLAGS="-L$withval/lib"
        GTEST_LIBS="-lgtest"
        acx_with_google_test="$withval"
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
    # If not using our default internal google test version then check user
    # provied version works
    if test $acx_with_google_test != "yes"; then
      AC_LANG_SAVE
      AC_LANG([C++])
    
      # Save and set compile flags for checks
      acx_with_google_test_cppflags_save=$CPPFLAGS
      acx_with_google_test_cxxflags_save=$CXXFLAGS
      acx_with_google_test_ldflags_save=$LDFLAGS
      acx_with_google_test_libs_save=$LIBS
      CPPFLAGS="$GTEST_CPPFLAGS $CPPFLAGS"
      CXXFLAGS="$GTEST_CXXFLAGS $CXXFLAGS"
      LDFLAGS="$GTEST_LDFLAGS $LDFLAGS"
      LIBS="$GTEST_LIBS $LIBS"
    
      AC_CHECK_HEADER([gtest/gtest.h], [], [AC_MSG_ERROR([Unable to compile with user provided Google Test.])])
    
      AC_MSG_CHECKING([for ::testing::InitGoogleTest in -lgtest])
      AC_LINK_IFELSE(
        [AC_LANG_PROGRAM(
          [[namespace testing { void InitGoogleTest(int*,char**); }]],
          [[int i = 0; char** c = 0; ::testing::InitGoogleTest(&i,c);]])
        ],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_ERROR(["Unable to link with user provided libgtest])])

      # Restore the original flags
      CPPFLAGS=$acx_with_google_test_cppflags_save
      CXXFLAGS=$acx_with_google_test_cxxflags_save
      LDFLAGS=$acx_with_google_test_ldflags_save
      LIBS=$acx_with_google_test_libs_save

      AC_LANG_RESTORE
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

AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_CACHE_CHECK([for std::shared_ptr], [acx_cv_std_shard_ptr], [
    AC_LANG_SAVE
    AC_LANG([C++])
    AC_TRY_COMPILE(
      [#include <memory>],
      [std::shared_ptr<int> p;],
      [acx_cv_std_shard_ptr=yes],
      [acx_cv_std_shard_ptr=no]
    )
    AC_LANG_RESTORE
  ])
  if test "$acx_cv_std_shard_ptr" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_SHARED_PTR],[1],[define if std::shared_ptr is available.])
  fi
  
  AC_CACHE_CHECK([for std::tr1::shared_ptr], [acx_cv_std_tr1_shard_ptr], [
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_TRY_COMPILE(
      [#include <memory>],
      [std::tr1::shared_ptr<int> p;],
      [acx_cv_std_tr1_shard_ptr=yes],
      [acx_cv_std_tr1_shard_ptr=no]
    )
    AC_LANG_RESTORE
  ])
  if test "$acx_cv_std_tr1_shard_ptr" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
  fi
  if test "$acx_cv_std_shard_ptr$acx_cv_std_tr1_shard_ptr" = nono; then
    AC_MSG_ERROR([Unable to find std::shared_ptr or std::tr1::shared_ptr. Reconfigure Madness with --with-boost option.])
  fi
])

AC_DEFUN([ACX_CHECK_TYPE_TRAITS], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  AC_CACHE_CHECK([for std type traits], [acx_cv_std_type_traits], [
    AC_TRY_LINK(
      [#include <type_traits>],
      [bool test = std::is_same<int, double>::value;],
      [acx_cv_std_type_traits=yes],
      [acx_cv_std_type_traits=no],
    )
  ])
  if test "$acx_cv_std_type_traits" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
  fi
  
  AC_CACHE_CHECK([for std::tr1 type traits], [acx_cv_std_tr1_type_traits], [
    AC_TRY_LINK(
      [#include <type_traits>],
      [bool test = std::tr1::is_same<int, double>::value;],
      [acx_cv_std_tr1_type_traits=yes],
      [acx_cv_std_tr1_type_traits=no]
    )
  ])
  if test "$acx_cv_std_tr1_type_traits" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1::shared_ptr is available.])
  fi
  if test "$acx_cv_std_type_traits$acx_cv_std_tr1_type_traits" = nono; then
    AC_MSG_ERROR([Unable to find std or std::tr1 type traits. Reconfigure Madness with --with-boost option.])
  fi
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_TR1],
[
  ACX_CHECK_SHARED_PTR
  ACX_CHECK_TYPE_TRAITS
])

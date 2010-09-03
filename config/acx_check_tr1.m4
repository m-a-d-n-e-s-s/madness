AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_CACHE_CHECK([for std::shared_ptr], [acx_cv_std_shard_ptr], [
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    AC_TRY_COMPILE(
      [#include <memory>],
      [std::shared_ptr<int> p;],
      [acx_cv_std_shard_ptr=yes],
      [acx_cv_std_shard_ptr=no]
    )
    AC_LANG_RESTORE
  ])
  if test "$acx_cv_std_shard_ptr" = yes; then
    AC_DEFINE([MADNESS_HAVE_STD_SHARED_PTR],[1],[define if std::shared_ptr is available.])
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
    AC_DEFINE([MADNESS_HAVE_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
  fi
])

AC_DEFUN([ACX_CHECK_TR1],
[
  ACX_CHECK_SHARED_PTR
])

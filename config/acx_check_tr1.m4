AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_CACHE_CHECK([for std::shared_ptr], [acx_cv_std_shard_ptr], [
    AC_TRY_COMPILE(
      [#include <memory>],
      [std::shared_ptr<int> p;],
      [acx_cv_std_shard_ptr=yes],
      [acx_cv_std_shard_ptr=no]
    )
  ])
  if test "$acx_cv_std_shard_ptr" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_SHARED_PTR],[1],[define if std::shared_ptr is available.])
  fi
  
  # Check for shared_ptr in std::tr1 namespace
  AC_CACHE_CHECK([for std::tr1::shared_ptr], [acx_cv_std_tr1_shard_ptr], [
    AC_TRY_COMPILE(
      [#include <memory>],
      [std::tr1::shared_ptr<int> p;],
      [acx_cv_std_tr1_shard_ptr=yes],
      [acx_cv_std_tr1_shard_ptr=no]
    )
  ])
  if test "$acx_cv_std_tr1_shard_ptr" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
  fi
  
  # Generate a configure error if we cannot find shared_ptr
  if test "$acx_cv_std_shard_ptr$acx_cv_std_tr1_shard_ptr" = nono; then
    AC_MSG_ERROR([std::shared_ptr and std::tr1::shared_ptr are not supported. Reconfigure Madness to use Boost with the --with-boost option.])
  fi
  
  #Check for std::make_shared and std::allocate_shared
  AC_CACHE_CHECK([for std::make_shared and std::allocate_shared], [acx_cv_std_make_shared], [
    AC_TRY_COMPILE(
      [#include <memory>],
      [using ::std::make_shared; using ::std::allocate_shared;],
      [acx_cv_std_make_shared=yes],
      [acx_cv_std_make_shared=no]
    )
  ])
  if test "$acx_cv_std_make_shared" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_MAKE_SHARED],[1],
      [define if std::make_shared and std::allocate_shared are available.])
  fi
 
  # make_shared and allocate_shared are not a part of TR1 so we do not check for it here.
  
  # Generate a configure error if we cannot find suitable make_shared and allocate_shared support.
  if test "$ac_cv_header_boost_make_shared_hpp$acx_cv_std_make_shared" = nono; then
    AC_MSG_ERROR([std::make_shared and boost::make_shared are not available. Reconfigure Madness to use Boost with the --with-boost option.])
  fi
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_TYPE_TRAITS], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  AC_CACHE_CHECK([for std type traits], [acx_cv_std_type_traits], [
    AC_TRY_COMPILE(
      [#include <type_traits>],
      [
      typedef std::is_same<int, double> sameT;
      typedef std::is_integral<int> integralT;
      typedef std::is_floating_point<int> is_floating_pointT;
      typedef std::is_fundamental<int> fundamentalT;
      typedef std::is_pointer<int> pointerT;
      typedef std::is_array<int> arrayT;
      typedef std::is_reference<int> referenceT;
      typedef std::is_const<int> constT;
      typedef std::add_const<int> add_constT;
      typedef std::remove_const<int> remove_constT;
      typedef std::remove_reference<int> remove_referenceT;
      ],
      [acx_cv_std_type_traits=yes],
      [acx_cv_std_type_traits=no],
    )
  ])
  if test "$acx_cv_std_type_traits" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
  fi
  
  AC_CACHE_CHECK([for std::tr1 type traits], [acx_cv_std_tr1_type_traits], [
    AC_TRY_COMPILE(
      [#include <type_traits>],
      [
      typedef std::tr1::is_same<int, double> sameT;
      typedef std::tr1::is_integral<int> integralT;
      typedef std::tr1::is_floating_point<int> is_floating_pointT;
      typedef std::tr1::is_fundamental<int> fundamentalT;
      typedef std::tr1::is_pointer<int> pointerT;
      typedef std::tr1::is_array<int> arrayT;
      typedef std::tr1::is_reference<int> referenceT;
      typedef std::tr1::is_const<int> constT;
      typedef std::tr1::add_const<int> add_constT;
      typedef std::tr1::remove_const<int> remove_constT;
      typedef std::tr1::remove_reference<int> remove_referenceT;
      ],
      [acx_cv_std_tr1_type_traits=yes],
      [acx_cv_std_tr1_type_traits=no]
    )
  ])
  if test "$acx_cv_std_tr1_type_traits" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
  fi
  if test "$acx_cv_std_type_traits$acx_cv_std_tr1_type_traits" = nono; then
    AC_MSG_ERROR([std and std::tr1 type traits are not supported. Reconfigure Madness to use Boost with the --with-boost option.])
  fi
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_ARRAY], [
  AC_LANG_SAVE
  AC_LANG([C++])

  AC_CACHE_CHECK([for std::array], [acx_cv_std_array], [
    AC_TRY_COMPILE(
      [#include <array>],
      [std::array<int, 3> a;],
      [acx_cv_std_array=yes],
      [acx_cv_std_array=no]
    )
  ])
  if test "$acx_cv_std_array" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_ARRAY],[1],[define if std::array is available.])
  fi
  
  AC_CACHE_CHECK([for std::tr1::array], [acx_cv_std_tr1_array], [
    AC_TRY_COMPILE(
      [#include <array>],
      [std::tr1::array<int, 3> a;],
      [acx_cv_std_tr1_array=yes],
      [acx_cv_std_tr1_array=no]
    )
  ])
  if test "$acx_cv_std_tr1_array" = yes; then
    AC_DEFINE([MADNESS_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])
  fi
  if test "$acx_cv_std_array$acx_cv_std_tr1_array" = nono; then
    AC_MSG_ERROR([std::array and std::tr1::array are not supported. Reconfigure Madness to use Boost with the --with-boost option.])
  fi
  
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_TR1],
[
  ACX_CHECK_SHARED_PTR
  ACX_CHECK_TYPE_TRAITS
  ACX_CHECK_ARRAY
])

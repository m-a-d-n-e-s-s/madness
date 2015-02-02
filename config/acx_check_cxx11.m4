AC_DEFUN([ACX_CHECK_VARIADIC_TEMPLATE], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for variadic template support])

  
  # Check for std::shared_ptr in <memory>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[template <typename ... T> void func(const T&... t) { }]],
        [[func(1, 2, 3);]]
      )
    ],
    [acx_variadic_template=yes],
    [acx_variadic_template=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_variadic_template])
  
  if test "$acx_variadic_template" = no; then
      AC_MSG_ERROR([$CXX does not support variadic templates])
  fi
])

AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for std::shared_ptr support])

  
  # Check for std::shared_ptr in <memory>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <memory>]],
        [[std::shared_ptr<int> p = std::make_shared<int>(1);]]
      )
    ],
    [acx_shared_ptr=yes],
    [acx_shared_ptr=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_shared_ptr])
  
  if test "$acx_shared_ptr" = no; then
      AC_MSG_ERROR([$CXX does not support std::shared_ptr])
  fi
  
])

AC_DEFUN([ACX_CHECK_TYPE_TRAITS], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for type traits support])

  
  # Check for type traits in <type_traits> and std namespace
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <type_traits>]],
        [[typedef std::is_same<int, double> sameT;]]
      )
    ],
    [acx_type_traits=yes],
    [acx_type_traits=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_type_traits])
  
  if test "$acx_type_traits" = no; then
      AC_MSG_ERROR([$CXX does not support type traits])
  fi
  
])

AC_DEFUN([ACX_CHECK_ARRAY], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for std::array support])

  
  # Check for type traits in <type_traits> and std namespace
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <array>]],
        [[std::array<int,10> a;]]
      )
    ],
    [acx_array=yes],
    [acx_array=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_array])
  
  if test "$acx_array" = no; then
      AC_MSG_ERROR([$CXX does not support std::array])
  fi
  
])

AC_DEFUN([ACX_CHECK_HASH], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for std::hash support])

  
  # Check for type traits in <type_traits> and std namespace
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <functional>]],
        [[std::hash<int> h; h(1);]]
      )
    ],
    [acx_hash=yes],
    [acx_hash=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_hash])
  
  if test "$acx_hash" = no; then
      AC_MSG_ERROR([$CXX does not support std::hash])
  fi
  
])

AC_DEFUN([ACX_CHECK_TUPLE], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for std::tuple support])

  
  # Check for type traits in <type_traits> and std namespace
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <tuple>]],
        [[std::tuple<int, double, unsigned long> t(1, 2.0, 3ul);]]
      )
    ],
    [acx_tuple=yes],
    [acx_tuple=no]
  )
  
  AC_LANG_RESTORE
  
  # post make_shared results
  AC_MSG_RESULT([$acx_tuple])
  
  if test "$acx_tuple" = no; then
      AC_MSG_ERROR([$CXX does not support std::tuple])
  fi
  
])

AC_DEFUN([ACX_CHECK_CXX11],
[
  ACX_CHECK_VARIADIC_TEMPLATE
  ACX_CHECK_SHARED_PTR
  ACX_CHECK_TYPE_TRAITS
  ACX_CHECK_ARRAY
  ACX_CHECK_HASH
  ACX_CHECK_TUPLE
])

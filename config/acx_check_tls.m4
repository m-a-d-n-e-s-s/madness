AC_DEFUN([ACX_CHECK_TLS],[
  # Check for thread local storage support.
  # thread_local, __thread_local, __declspec(thread)
  
  AC_LANG_PUSH([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([for C++ thread_local keyword])
  acx_check_tls=no
  
  # Check for the key word thread_local (C++11)
  AC_LINK_IFELSE(  
    [
      AC_LANG_PROGRAM(
        [[thread_local int i = 0;]],
        [[i = 1;]]
      )
    ],
    [acx_check_tls="thread_local"]
  )
  
  # Check for the key word __thread
  if test "$acx_check_tls" = no; then
    AC_LINK_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[__thread int i = 0;]],
          [[i = 1;]]
        )
      ],
      [
        acx_check_tls="__thread"
        AC_DEFINE([thread_local],[__thread],[Define the thread_local key word.])
      ]
    )
  fi
  
   # Check for the key word __declspec(thread)
#  if test "$acx_check_tls" = no; then
#    AC_LINK_IFELSE(
#      [
#        AC_LANG_PROGRAM(
#          [[__declspec(thread) int i = 0;]],
#          [[i = 1;]]
#        )
#      ],
#      [
#        acx_check_tls="__declspec(thread)"
#        AC_DEFINE([thread_local],[__declspec(thread)],[Define the thread_local key word.])
#      ]
#    )
#  fi
  
  if test "$acx_check_tls" = no; then
    AC_DEFINE([thread_local],[],[Define the thread_local key word.])
  fi

  
  AC_MSG_RESULT([$acx_check_tls])
  
  AC_LANG_POP
])


AC_DEFUN([ACX_ENABLE_DEBUGGING], [
  ac_enable_debugging="no"
  ac_enable_debugging_flags=""
  AC_ARG_ENABLE([debugging],
    [AC_HELP_STRING([--enable-debugging@<:@=yes|no|OPTION@:>@],
      [Enable debugging C and C++ compilers @<:@default=no@:>@]) ],
    [
      case $enableval in
        yes)
          ac_enable_debugging="yes"
          ac_enable_debugging_flags="-g"
        ;;
        no)
        ;;
        *)
          ac_enable_debugging="yes"
          ac_enable_debugging_flags="-g$enableval"
        ;;
      esac
    ])

  if test $ac_enable_debugging != no; then
    ACX_CHECK_COMPILER_FLAG([C], [CFLAGS], [$ac_enable_debugging_flags],
      [CFLAGS="$CFLAGS $ac_enable_debugging_flags"],
      [AC_MSG_WARN([$CC does not accept $ac_enable_debugging_flags, no debugging flags will be used.])])
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$ac_enable_debugging_flags],
      [CXXFLAGS="$CXXFLAGS $ac_enable_debugging_flags"],
      [AC_MSG_WARN([$CXX does not accept $ac_enable_debugging_flags, no debugging flags will be used.])])
  fi
])
# This function is used to add debug flags (e.g. -g) to CFLAGS and CXXFLAGS 
# environment variables. Users are expected to specify their own debug flags for
# special use cases by adding appropriate values to CFLAGS and CXXFLAGS.
AC_DEFUN([ACX_ENABLE_DEBUGGING], [
  acx_enable_debugging="no"
  acx_enable_debugging_flags=""

  # Allow the user to enable or disable debugging flag
  AC_ARG_ENABLE([debugging],
    [AC_HELP_STRING([--enable-debugging@<:@=yes|no|LEVEL@:>@],
      [Enable debugging C and C++ compilers. You can also specify debug level (e.g. 3). @<:@default=no@:>@]) ],
    [
      case $enableval in
        yes)
          acx_enable_debugging="yes"
          acx_enable_debugging_flags="-g"
        ;;
        no)
        ;;
        *)
          acx_enable_debugging="yes"
          acx_enable_debugging_flags="-g$enableval"
        ;;
      esac
    ])

  # Test the flags and add them to flag variables if successful.
  if test $acx_enable_debugging != no; then
    ACX_CHECK_COMPILER_FLAG([C], [CFLAGS], [$acx_enable_debugging_flags],
      [CFLAGS="$CFLAGS $acx_enable_debugging_flags"],
      [AC_MSG_WARN([$CC does not accept $acx_enable_debugging_flags, no debugging flags will be used.])])
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$acx_enable_debugging_flags],
      [CXXFLAGS="$CXXFLAGS $acx_enable_debugging_flags"],
      [AC_MSG_WARN([$CXX does not accept $acx_enable_debugging_flags, no debugging flags will be used.])])
  fi
])
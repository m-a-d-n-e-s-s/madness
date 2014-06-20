# This function is used to add the gentensor flags and CXXFLAGS. These flags are 
# necessary for the correlated quantum chemistry (e.g. mp2), but gentensor
# can't handle complex tensors. In this case disable gentensor (default behavior)
AC_DEFUN([ACX_ENABLE_GENTENSOR], [
  acx_enable_gentensor="no"
  acx_enable_gentensor_flags=""

  # Allow the user to enable or disable gentensor flag
  AC_ARG_ENABLE([gentensor],
    [AC_HELP_STRING([--enable-gentensor@<:@=yes|no],
      [Enable gentensor C and C++ compilers.]) ],
    [
      case $enableval in
        yes)
          acx_enable_gentensor="yes"
          acx_enable_gentensor_flags="-DUSE_GENTENSOR"
        ;;
        no)
        ;;
        *)
        ;;
      esac
    ])

  # Test the flags and add them to flag variables if successful.
  if test $acx_enable_gentensor != no; then
    ACX_CHECK_COMPILER_FLAG([C], [CFLAGS], [$acx_enable_gentensor_flags],
      [CFLAGS="$CFLAGS $acx_enable_gentensor_flags"],
      [AC_MSG_WARN([$CC does not accept $acx_enable_gentensor_flags, no gentensor flags will be used.])])
    ACX_CHECK_COMPILER_FLAG([C++], [CXXFLAGS], [$acx_enable_gentensor_flags],
      [CXXFLAGS="$CXXFLAGS $acx_enable_gentensor_flags"],
      [AC_MSG_WARN([$CXX does not accept $acx_enable_gentensor_flags, no gentensor flags will be used.])])
  fi
])

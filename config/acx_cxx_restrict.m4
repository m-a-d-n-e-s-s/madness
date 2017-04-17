dnl @synopsis ACX_CXX_RESTRICT
dnl
dnl This macro determines whether the C++ compiler supports the
dnl "restrict" keyword introduced in C99, or an equivalent. Does
dnl nothing if the compiler accepts the keyword. Otherwise, if the
dnl compiler supports an equivalent (like GCC's __restrict__) defines
dnl "restrict" to be that. Otherwise, defines "restrict" to be empty.
dnl
dnl @category Obsolete
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2005-05-31
dnl @license GPLWithACException

AC_DEFUN([ACX_CXX_RESTRICT],
[AC_CACHE_CHECK([for C++ restrict keyword], acx_cv_cxx_restrict,
[acx_cv_cxx_restrict=unsupported
 AC_LANG_SAVE
 AC_LANG_PUSH([C++])
 # Try the official restrict keyword, then gcc's __restrict__, then
 # SGI's __restrict.  __restrict has slightly different semantics than
 # restrict (it's a bit stronger, in that __restrict pointers can't
 # overlap even with non __restrict pointers), but I think it should be
 # okay under the circumstances where restrict is normally used.
 for acx_kw in restrict __restrict__ __restrict; do
   AC_TRY_COMPILE([], [float * $acx_kw x;], [acx_cv_cxx_restrict=$acx_kw; break])
 done
 AC_LANG_RESTORE
])
 if test "$acx_cv_cxx_restrict" != "restrict"; then
   acx_kw="$acx_cv_cxx_restrict"
   if test "$acx_kw" = unsupported; then acx_kw=""; fi
   AC_DEFINE_UNQUOTED(MADNESS_RESTRICT, $acx_kw, [Define to equivalent of C99 restrict keyword, or to nothing if this is not supported.  Do not define if restrict is supported directly.])
  fi
  AC_MSG_RESULT([$acx_cv_cxx_restrict])
  AC_LANG_POP
])

AC_DEFUN([ACX_POSIX_MEMALIGN], [
    AC_MSG_CHECKING([for posix_memalign])
    AC_CHECK_FUNC([posix_memalign],[gotpm=1], [gotpm=0])
    if test $gotpm==1; then
        AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [Set if have posix_memalign])
    else
        AC_LANG_PUSH([C++])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
extern "C"  int posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
int main() {
    void *m;
    posix_memalign(&m, 16, 1024);
    return 0;
}
        ]])],
        [ AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [Set if have posix_memalign])
          AC_DEFINE(MISSING_POSIX_MEMALIGN_PROTO, [1], [Set if the posix_memalign prototype is missing]) 
          AC_MSG_RESULT([yes]) ],
         [AC_MSG_RESULT([no])])
        AC_LANG_POP([C++])
    fi
])

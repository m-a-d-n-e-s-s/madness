AC_DEFUN([ACX_POSIX_MEMALIGN], [
    AC_MSG_CHECKING([for existence of posix_memalign])
    AC_CHECK_FUNC([posix_memalign],[gotpm=1], [gotpm=0])
    if test $gotpm = 1; then
        AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [Set if have posix_memalign])
    else
        AC_MSG_WARN([[   posix_memalign NOT FOUND ... enabling override of new/delete ... THIS WILL BE SLOW ]])
        AC_DEFINE([WORLD_GATHER_MEM_STATS], [1], [Set if MADNESS gathers memory statistics])
    fi

    if test $gotpm = 1; then
        AC_MSG_CHECKING([if missing declaration of posix_memalign in stdlib.h])
        AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([[
#include <stddef.h>
#include <stdlib.h>
extern "C"  int posix_memalign(void **memptr, size_t alignment, size_t size);
int main() {
    void *m;
    posix_memalign(&m, 16, 1024);
    return 0;
}
        ]],
        [ AC_DEFINE(MISSING_POSIX_MEMALIGN_PROTO, [1], [Set if the posix_memalign prototype is missing]) 
          AC_MSG_RESULT([yes]) ],
         [AC_MSG_RESULT([no])])
        AC_LANG_POP([C++])
    fi
])

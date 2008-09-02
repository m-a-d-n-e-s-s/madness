AC_DEFUN([ACX_GNU_HASHMAP], [
    AC_LANG_PUSH([C++])
    AC_MSG_CHECKING([for GNU hashmap and associated namespace])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[        
#include <ext/hash_map>
using __gnu_cxx::hash_map; ]])],

        [ AC_DEFINE(HAVE_GNU_HASHMAP,[1],[Enable if have GNU hashmap]) 
          AC_DEFINE(GNU_HASHMAP_NAMESPACE, __gnu_cxx, [GNU hashmap namespace]) 
          AC_MSG_RESULT([yes])],

        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[        
#include <ext/hash_map>
using _STLP_STD::hash_map; ]])],
            [ AC_DEFINE(HAVE_GNU_HASHMAP,[1],[Enable if have GNU hashmap]) 
              AC_DEFINE(GNU_HASHMAP_NAMESPACE, _STLP_STD, [GNU hashmap namespace]) 
              AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])]))
    AC_LANG_POP([C++])
])





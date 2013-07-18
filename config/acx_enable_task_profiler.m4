AC_DEFUN([ACX_ENABLE_TASK_PROFILER], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  acx_enable_task_profiler=""
  
  # Allow the user to enable or disable warnings
  AC_ARG_ENABLE([task-profiler],
    [AC_HELP_STRING([--enable-task-profiler],
      [Enable task profiler that collects per-task start and stop times.])],
    [
      case $enableval in
      yes)
        acx_enable_task_profiler="yes"
      ;;
      no)
        acx_enable_task_profiler="no"
      ;;
      *)
        acx_enable_task_profiler="yes"
      esac
    ],
    [acx_enable_task_profiler="no"]
  )
  
  # Automatically specify the warning flags for known compilers. 
  if test $acx_enable_task_profiler != "no"; then
    AC_CHECK_HEADER([execinfo.h], [], [AC_MSG_ERROR([execinfo.h is required by the task profiler.])])
    AC_CHECK_HEADER([cxxabi.h], [], [AC_MSG_ERROR([cxxabi.h is required by the task profiler.])])
    AC_DEFINE([MADNESS_TASK_PROFILING],[1],[Define to enable task profiler.])
    if test $ON_A_MAC = "no"; then
      case $acx_enable_optimal_compiler in
        GNU)
          LDFLAGS="$LDFLAGS -rdynamic"
        ;;
        
        clang)
          LDFLAGS="$LDFLAGS -rdynamic"
        ;;
    
        Pathscale)
        ;;
    
        Portland)
        ;;
    
        Intel)
        ;;
          
        IBM)
        ;;
      
        *)
        ;;
      esac
    fi
  fi
  
  AC_LANG_RESTORE
])

if(ENABLE_TASK_PROFILER OR ENABLE_GPERFTOOLS)

  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # libunwind is not supported on OS X or Windows
    set(LIBUNWIND_FOUND FALSE)
  else()
    find_package(Libunwind)
  endif()
    
  # Set the output variables
  if(LIBUNWIND_FOUND)
    set(MADNESS_HAS_LIBUNWIND 1)
  endif()
      
endif()
if(WITH_TASK_PROFILER OR WITH_GPERFTOOLS)

  if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    if(WITH_GPERFTOOLS)
      # Libunwind 0.99 is required when using gperftools on Linux 
      find_package(Libunwind 0.99 REQUIRED)
    else()
      # Libunwind is optional if we are only using the task profiler
      find_package(Libunwind)
    endif()
  else()
    # libunwind is not supported on OS X or Windows
    set(LIBUNWIND_FOUND FALSE)
  endif()
    
  # Set the output variables
  if(LIBUNWIND_FOUND)
    set(MADNESS_HAS_LIBUNWIND 1)
  endif()
      
endif()
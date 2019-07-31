if(NOT LIBUNWIND_FOUND AND (ENABLE_LIBUNWIND OR ENABLE_TASK_PROFILER))

  find_package(Libunwind)
    
  # Set the output variables
  if(LIBUNWIND_FOUND)
    set(MADNESS_HAS_LIBUNWIND 1)
  endif()
      
endif()
if(ENABLE_GPERFTOOLS OR ENABLE_TCMALLOC_MINIMAL)
  
  if(ENABLE_GPERFTOOLS)
    find_package(Gperftools COMPONENTS tcmalloc OPTIONAL_COMPONENTS profiler)
  else()
    find_package(Gperftools REQUIRED COMPONENTS tcmalloc_minimal)
  endif()

  # Set the config.h variables
  if(GPERFTOOLS_FOUND AND ENABLE_TCMALLOC_MINIMAL)
    set(MADNESS_HAS_GOOGLE_PERF_MINIMAL 1)
  endif()
  if(LIBUNWIND_FOUND)
    set(MADNESS_HAS_LIBUNWIND 1)
  endif()
      
endif()
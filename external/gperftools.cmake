if(WITH_GPERFTOOLS)
  find_package(Gperftools)

  if(GPERFTOOLS_FOUND)
    if(GPERFTOOLS_tcmalloc_minimal_FOUND AND NOT GPERFTOOLS_tcmalloc_FOUND)
      set(MADNESS_HAS_GOOGLE_PERF_MINIMAL 1)
    endif()
  endif()
      
endif()
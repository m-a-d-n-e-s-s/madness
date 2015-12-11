if(ENABLE_TBB)
  find_package(TBB)
  
  if(TBB_FOUND)
    set(HAVE_INTEL_TBB 1)
  endif()
  
endif()
if(ENABLE_TBB)
  find_package(TBB 4.3.5)
  
  if(TBB_FOUND)
    set(HAVE_INTEL_TBB 1)
    message(STATUS "Found TBB: ${TBB_LIBRARIES}")
  endif()
  
endif()
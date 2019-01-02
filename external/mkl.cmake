if(ENABLE_MKL)
  find_package(MKL)
  
  if(MKL_FOUND)
    set(HAVE_INTEL_MKL 1)
  endif()
  
endif()
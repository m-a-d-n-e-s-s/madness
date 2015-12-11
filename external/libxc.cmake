if(ENABLE_LIBXC)

  find_package(Libxc)
    
  # Set the output variables
  if(LIBXC_FOUND)
    set(MADNESS_HAS_LIBXC 1)
  endif()
      
endif()
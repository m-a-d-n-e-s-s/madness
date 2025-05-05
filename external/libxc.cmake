if(ENABLE_LIBXC)

  find_package(Libxc CONFIG)
    
  # Set the output variables
  if(TARGET Libxc::xc)
    set(MADNESS_HAS_LIBXC 1)
  endif()
      
endif()

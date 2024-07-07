if(ENABLE_INTEGRATORXX)

  find_package(IntegratorXX)

  # Set the output variables
  if(INTEGRATORXX_FOUND)
    set(MADNESS_HAS_INTEGRATORXX 1)
  endif()
      
endif()

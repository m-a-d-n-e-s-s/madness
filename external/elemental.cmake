if(ENABLE_ELEMENTAL)
  find_package(Elemental 0.84 EXACT REQUIRED 
      COMPONENTS pmrrr OPTIONAL_COMPONENTS lapack-addons)
      
  set(MADNESS_HAS_ELEMENTAL 1)
  if(ELEMENTAL_VERSION VERSION_LESS 0.85)
    set(HAVE_ELEMENTAL_H 1)
  else()
    set(HAVE_EL_H 1)
  endif()

endif()
if(ENABLE_BOOST)

  find_package(Boost 1.34  REQUIRED COMPONENTS math_tr1)

  # Set the output variables
  if(BOOST_FOUND)
	  set(MADNESS_HAS_BOOST 1)
  endif()

endif()


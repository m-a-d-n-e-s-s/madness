if(ENABLE_BOOST)

    cmake_minimum_required(VERSION 3.15.0)  # for Boost::headers
    find_package(Boost 1.35 REQUIRED)

  # Set the output variables
  if(TARGET Boost::headers)
	  set(MADNESS_HAS_BOOST 1)
  endif()

endif()


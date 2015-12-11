if(ENABLE_PAPI)
  
  find_package(Papi REQUIRED)
      
  # Set the output variables
  set(HAVE_PAPI 1)

endif()

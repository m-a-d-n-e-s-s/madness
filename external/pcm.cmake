if(ENABLE_PCM)

  find_package(PCM COMPONENTS pcm)
    
  # Set the output variables
  if(PCM_FOUND)
    set(MADNESS_HAS_PCM 1)
  endif()
      
endif()

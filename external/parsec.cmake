if(ENABLE_PARSEC)
  find_package(PARSEC REQUIRED)

  if(PARSEC_FOUND)
    set(HAVE_PARSEC 1)
  endif()
endif()
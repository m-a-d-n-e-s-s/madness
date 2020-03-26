if(ENABLE_PARSEC)
  find_package(PARSEC REQUIRED)

  if(TARGET PaRSEC::parsec)
    set(HAVE_PARSEC 1)
  else()
    message(FATAL_ERROR "PaRSEC found but PaRSEC::parsec target not provided; using legacy PaRSEC build?")
  endif()
endif()
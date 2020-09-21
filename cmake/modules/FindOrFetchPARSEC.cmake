find_package(PaRSEC CONFIG QUIET COMPONENTS parsec HINTS ${PaRSEC_ROOT_DIR})

if (NOT TARGET PaRSEC::parsec)

  # configure PaRSEC
  set(SUPPORT_FORTRAN OFF CACHE BOOL "Disable Fortran support in PaRSEC")
  set(CMAKE_CROSSCOMPILING OFF)
  set(CMAKE_SYSTEM_PROCESSOR ${CMAKE_HOST_SYSTEM_PROCESSOR})
  set(PARSEC_WITH_DEVEL_HEADERS ON CACHE BOOL "Install PaRSEC headers")

  FetchContent_Declare(
          PARSEC
          GIT_REPOSITORY      https://bitbucket.org/herault/parsec.git
          GIT_TAG             ${MADNESS_TRACKED_PARSEC_TAG}
  )
  FetchContent_MakeAvailable(PARSEC)
  FetchContent_GetProperties(PARSEC
          SOURCE_DIR PARSEC_SOURCE_DIR
          BINARY_DIR PARSEC_BINARY_DIR
          )
  set_property(DIRECTORY ${PARSEC_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL TRUE)

endif(NOT TARGET PaRSEC::parsec)

# postcond check
if (NOT TARGET PaRSEC::parsec)
  message(FATAL_ERROR "FindOrFetchPARSEC could not make PaRSEC::parsec target available")
endif(NOT TARGET PaRSEC::parsec)

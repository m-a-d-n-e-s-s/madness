find_package(PaRSEC CONFIG QUIET COMPONENTS parsec HINTS ${PaRSEC_ROOT_DIR})

if (TARGET PaRSEC::parsec)

  message(STATUS "Found PaRSEC CONFIG at ${PaRSEC_CONFIG}")

else (TARGET PaRSEC::parsec)

  # configure PaRSEC
  set(SUPPORT_FORTRAN OFF CACHE BOOL "Disable Fortran support in PaRSEC")
  set(CMAKE_CROSSCOMPILING OFF)
  set(CMAKE_SYSTEM_PROCESSOR ${CMAKE_HOST_SYSTEM_PROCESSOR})
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(PARSEC_WITH_DEVEL_HEADERS ON CACHE BOOL "Install PaRSEC headers")
  set(BUILD_TOOLS OFF CACHE BOOL "Don't build the tools by default")

  FetchContent_Declare(
          PARSEC
          GIT_REPOSITORY      https://bitbucket.org/schuchart/parsec.git
          GIT_TAG             ${MADNESS_TRACKED_PARSEC_TAG}
  )
  FetchContent_MakeAvailable(PARSEC)
  FetchContent_GetProperties(PARSEC
          SOURCE_DIR PARSEC_SOURCE_DIR
          BINARY_DIR PARSEC_BINARY_DIR
          )

  # this is where PaRSECConfig.cmake will end up
  # must be in sync with the "install(FILES ...PaRSECConfig.cmake" statement in PaRSEC source
  set(PaRSEC_CONFIG "${CMAKE_INSTALL_PREFIX}/share/cmake/parsec/PaRSECConfig.cmake" CACHE INTERNAL "The location of installed PaRSECConfig.cmake file")
endif(TARGET PaRSEC::parsec)

# postcond check
if (NOT TARGET PaRSEC::parsec)
  message(FATAL_ERROR "FindOrFetchPARSEC could not make PaRSEC::parsec target available")
endif(NOT TARGET PaRSEC::parsec)

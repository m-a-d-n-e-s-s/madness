if (NOT TARGET PaRSEC::parsec)
  find_package(PaRSEC CONFIG QUIET COMPONENTS parsec HINTS ${PaRSEC_ROOT_DIR})
  if (TARGET PaRSEC::parsec)
    message(STATUS "Found PaRSEC CONFIG at ${PaRSEC_CONFIG}")
  endif (TARGET PaRSEC::parsec)
endif (NOT TARGET PaRSEC::parsec)

if (NOT TARGET PaRSEC::parsec)

  # configure PaRSEC
  set(SUPPORT_FORTRAN OFF CACHE BOOL "Disable Fortran support in PaRSEC")
  set(CMAKE_CROSSCOMPILING OFF)
  set(CMAKE_SYSTEM_PROCESSOR ${CMAKE_HOST_SYSTEM_PROCESSOR})

  FetchContent_Declare(
          PARSEC
          GIT_REPOSITORY     https://github.com/therault/parsec.git
          GIT_TAG            ${MADNESS_TRACKED_PARSEC_TAG}
  )
  FetchContent_MakeAvailable(PARSEC)
  FetchContent_GetProperties(PARSEC
          SOURCE_DIR PARSEC_SOURCE_DIR
          BINARY_DIR PARSEC_BINARY_DIR
          )

  # disable unity builds for parsec libs
  if (CMAKE_UNITY_BUILD)
    set(_parsec_libs_to_not_unity_build parsec-base-obj parsec-ptgpp parsec)
    foreach (_parsec_lib IN LISTS _parsec_libs_to_not_unity_build)
      if(TARGET ${_parsec_lib})
        set_target_properties(${_parsec_lib} PROPERTIES UNITY_BUILD OFF)
        message(STATUS "Will disable unity-build for target ${_parsec_lib}")
      else()
        message(FATAL_ERROR "FindOrFetchPaRSEC missing target ${_parsec_lib}")
      endif()
    endforeach ()
  endif(CMAKE_UNITY_BUILD)

  # this is where PaRSECConfig.cmake will end up
  # must be in sync with the "install(FILES ...PaRSECConfig.cmake" statement in PaRSEC source
  set(PaRSEC_CONFIG "${CMAKE_INSTALL_PREFIX}/share/cmake/parsec/PaRSECConfig.cmake" CACHE INTERNAL "The location of installed PaRSECConfig.cmake file")

  # export parsec targets from the build tree for the same to be possible for madness targets
  export(EXPORT parsec-targets FILE "${PROJECT_BINARY_DIR}/parsec-targets.cmake")

endif(NOT TARGET PaRSEC::parsec)

# postcond check
if (NOT TARGET PaRSEC::parsec)
  message(FATAL_ERROR "FindOrFetchPARSEC could not make PaRSEC::parsec target available")
endif(NOT TARGET PaRSEC::parsec)

if (NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json ${MADNESS_TRACKED_NLOHMANN_JSON_VERSION} QUIET CONFIG)
  if (TARGET nlohmann_json::nlohmann_json)
    message(STATUS "Found nlohmann_json CONFIG at ${nlohmann_json_CONFIG}")
  endif (TARGET nlohmann_json::nlohmann_json)
endif (NOT TARGET nlohmann_json::nlohmann_json)

if (NOT TARGET nlohmann_json::nlohmann_json)

  cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
  include(FetchContent)
  FetchContent_Declare(
          nlohmann_json
          GIT_REPOSITORY     https://github.com/nlohmann/json.git
          GIT_TAG            ${MADNESS_TRACKED_NLOHMANN_JSON_TAG}
          GIT_SHALLOW        TRUE
  )

  # configure nlohmann_json: header-only, no tests, but install+export so that
  # MADNESS targets that depend on it can be exported/installed in turn
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  set(JSON_Install ON CACHE INTERNAL "")

  FetchContent_MakeAvailable(nlohmann_json)

  # this is where nlohmann_jsonConfig.cmake will end up once installed;
  # must stay in sync with nlohmann_json's own install rules
  set(nlohmann_json_CONFIG "${CMAKE_INSTALL_PREFIX}/share/cmake/nlohmann_json/nlohmann_jsonConfig.cmake" CACHE INTERNAL "The location of installed nlohmann_jsonConfig.cmake file")

  # N.B. nlohmann_json already exports its own targets from the build tree
  # (nlohmann_jsonTargets, courtesy of JSON_Install) so that MADNESS targets
  # depending on it can be exported in turn -- do NOT re-export them here, or
  # the nlohmann_json target lands in two export sets and CMake refuses to
  # export MADworld ("cannot depend upon another target ... exported in more
  # than one export set").

endif (NOT TARGET nlohmann_json::nlohmann_json)

# postcond check
if (TARGET nlohmann_json::nlohmann_json)
  set(MADNESS_HAS_NLOHMANN_JSON ON CACHE BOOL "MADNESS has access to nlohmann_json")
else (TARGET nlohmann_json::nlohmann_json)
  message(FATAL_ERROR "FindOrFetchNlohmannJson could not make nlohmann_json::nlohmann_json target available")
endif (TARGET nlohmann_json::nlohmann_json)

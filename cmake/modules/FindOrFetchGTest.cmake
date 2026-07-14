# Provides a MADgtest INTERFACE target (gtest, no gtest_main) for MADNESS unit tests.
#
# Discover-or-fetch, per the MADNESS FindOrFetch* convention
# (cf. FindOrFetchCereal.cmake, FindOrFetchPARSEC.cmake):
#   1. reuse a GTest target an enclosing project already defined
#   2. find_package(GTest) on the system (conda, distro, installed)
#   3. FetchContent a pinned googletest as a last resort
#
# MADgtest deliberately does NOT pull in gtest_main: MADNESS test drivers
# provide their own main().
if (NOT TARGET MADgtest)

  # 1./2. already provided by a parent project, or found on the system
  if (NOT (TARGET GTest::gtest OR TARGET GTest::GTest))
    find_package(GTest CONFIG QUIET)
    if (NOT (TARGET GTest::gtest OR TARGET GTest::GTest))
      find_package(GTest MODULE QUIET)
    endif()
  endif()

  # 3. fetch if still not found
  if (NOT (TARGET GTest::gtest OR TARGET GTest::GTest))
    cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
    include(FetchContent)
    set(MADNESS_FETCH_GTEST_TAG "v1.14.0"
        CACHE STRING "googletest git tag fetched when GTest is not found")
    message(STATUS "GTest not found; fetching googletest ${MADNESS_FETCH_GTEST_TAG}")
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)  # never install into MADNESS's tree
    set(BUILD_GMOCK   OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        ${MADNESS_FETCH_GTEST_TAG})
    FetchContent_MakeAvailable(googletest)  # gtest built on demand, not part of `all`
  endif()

  # normalize to a single stable name the rest of the tree links against
  add_library(MADgtest INTERFACE)
  if (TARGET GTest::gtest)         # CMake >= 3.20 casing / fetched googletest
    target_link_libraries(MADgtest INTERFACE GTest::gtest)
  elseif (TARGET GTest::GTest)     # CMake < 3.20 module casing
    target_link_libraries(MADgtest INTERFACE GTest::GTest)
  elseif (TARGET gtest)            # raw fetched target fallback
    target_link_libraries(MADgtest INTERFACE gtest)
  else()
    message(FATAL_ERROR "FindOrFetchGTest: no usable gtest target after discover/fetch")
  endif()

endif()

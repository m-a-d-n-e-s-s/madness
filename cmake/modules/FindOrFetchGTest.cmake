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
    set(MADNESS_GTEST_TARGET GTest::gtest)
  elseif (TARGET GTest::GTest)     # CMake < 3.20 module casing
    set(MADNESS_GTEST_TARGET GTest::GTest)
  elseif (TARGET gtest)            # raw fetched target fallback
    set(MADNESS_GTEST_TARGET gtest)
  else()
    message(FATAL_ERROR "FindOrFetchGTest: no usable gtest target after discover/fetch")
  endif()
  target_link_libraries(MADgtest INTERFACE ${MADNESS_GTEST_TARGET})

  # Carry an rpath to a preinstalled shared gtest.
  #
  # CMake decides whether a consumer needs an rpath entry from the imported
  # target's IMPORTED_SONAME, not from the dylib/so on disk. Some packagers
  # relocate the library to @rpath/@loader_path after generating their export
  # file, leaving a bare soname behind (conda-forge's gtest is the case we hit:
  # install_name @rpath/libgtest.<v>.dylib, IMPORTED_SONAME libgtest.<v>.dylib).
  # CMake then emits no rpath, every test links cleanly, and every test dies at
  # startup with "Library not loaded: @rpath/libgtest...". Adding the directory
  # unconditionally is harmless when the rpath was not needed.
  get_target_property(_mad_gtest_imported ${MADNESS_GTEST_TARGET} IMPORTED)
  get_target_property(_mad_gtest_type ${MADNESS_GTEST_TARGET} TYPE)
  if (UNIX AND _mad_gtest_imported AND _mad_gtest_type STREQUAL "SHARED_LIBRARY")
    target_link_options(MADgtest INTERFACE
        "LINKER:-rpath,$<TARGET_FILE_DIR:${MADNESS_GTEST_TARGET}>")
  endif()
  unset(_mad_gtest_imported)
  unset(_mad_gtest_type)

endif()

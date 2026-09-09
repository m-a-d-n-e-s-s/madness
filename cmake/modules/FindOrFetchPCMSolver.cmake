# - Fetch and build PCMSolver, the polarizable continuum model library
#
# Included by external/pcm.cmake only after both discovery paths (the
# PCMSolverConfig package config, then the FindPCM module) have come up empty.
#
# On success:
#   PCM_FOUND          - TRUE
#   PCM_LIBRARIES      - PCMSolver::pcm, an ALIAS for the in-tree pcm-static target
#   PCM_INCLUDE_DIRS   - empty: the include directories ride on that target
#   MADNESS_PCM_FETCHED- TRUE, so the rest of the harness can tell this path apart
#                        from a PCMSolver found on the system
#
# On a missing prerequisite this prints why and returns without setting
# PCM_FOUND: PCM is an optional feature, and a machine without a Fortran
# compiler should still get a working MADNESS out of the default configure.

if (TARGET PCMSolver::pcm)
  return()
endif ()

cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
include(FetchContent)
include(CheckLanguage)

set(_pcm_missing)

# PEDRA, PCMSolver's cavity generator, is Fortran. Probe rather than
# enable_language() outright so a machine without gfortran degrades to "no PCM"
# instead of a failed configure -- and clear the cache entry check_language()
# leaves behind on failure, which would otherwise poison a parent project
# (TiledArray, MPQC) that enables Fortran after pulling MADNESS in.
if (NOT DEFINED CMAKE_Fortran_COMPILER)
  check_language(Fortran)
  if (NOT CMAKE_Fortran_COMPILER)
    unset(CMAKE_Fortran_COMPILER CACHE)
    list(APPEND _pcm_missing "a Fortran compiler")
  endif ()
endif ()

# Boost headers (odeint, math) are a hard requirement of PCMSolver's Green's
# functions. Probed by header rather than find_package(Boost) to stay clear of
# CMP0167 -- MADNESS's own Boost hook is off by default (ENABLE_BOOST), so
# there is usually no Boost target to reuse here.
find_path(PCM_BOOST_INCLUDE_DIR NAMES boost/version.hpp)
mark_as_advanced(PCM_BOOST_INCLUDE_DIR)
if (NOT PCM_BOOST_INCLUDE_DIR)
  list(APPEND _pcm_missing "Boost headers (>= 1.54)")
endif ()

find_package(ZLIB QUIET)
if (NOT ZLIB_FOUND)
  list(APPEND _pcm_missing "zlib")
endif ()

if (_pcm_missing)
  list(JOIN _pcm_missing ", " _pcm_missing)
  message(STATUS "PCMSolver not found and cannot be built here (missing: ${_pcm_missing}); "
                 "the `pcm` keyword will be unavailable. Install PCMSolver and point "
                 "-DPCM_ROOT_DIR at the prefix, supply the missing prerequisite, or "
                 "configure -DENABLE_PCM=OFF to stop looking.")
  return()
endif ()
unset(_pcm_missing)

# The Fortran runtime is what actually makes the static libpcm.a linkable from
# MADNESS's Fortran-free C++ targets: CMake does not carry Fortran's implicit
# link libraries across an OBJECT library, which is how PCMSolver assembles
# pcm-static, so the consumer link line comes up short of __gfortran_* symbols.
# Enabling the language here (rather than at the top of the project) keeps it
# out of every configure that does not fetch PCMSolver -- notably out of
# external/lapack.cmake's detection, which runs earlier and deliberately probes
# Fortran symbols from C.
enable_language(Fortran)

# PCMSolver v1.3.0 (2020) predates the toolchains it now has to build with; see
# cmake/patches/pcmsolver-v1.3.0.cmake for what each hunk fixes and why.
set(MADNESS_TRACKED_PCMSOLVER_TAG "${MADNESS_TRACKED_PCMSOLVER_TAG}" CACHE STRING
    "The tag/commit of the PCMSolver repository to track")
if (NOT MADNESS_TRACKED_PCMSOLVER_TAG)
  message(FATAL_ERROR "MADNESS_TRACKED_PCMSOLVER_TAG is empty; external/versions.cmake must set it")
endif ()

# PCMSolver's own options, pre-seeded so its option()/cmake_dependent_option()
# calls adopt them. Not FORCEd: a value already in the cache -- the user's, or a
# parent project's for the generically-named ones -- wins.
#   STATIC_LIBRARY_ONLY  a shared libpcm would need install-tree RPATH handling
#                        that MADNESS, static by default, does not otherwise do
#   BUILD_STANDALONE     the run_pcm driver; MADNESS uses the C API only
#   ENABLE_TESTS         PCMSolver's Catch suite would join MADNESS's CTest set
#   TEST_Fortran_API     ditto, for its Fortran 90 binding tests
set(STATIC_LIBRARY_ONLY ON CACHE BOOL "Build PCMSolver as a static library only")
set(BUILD_STANDALONE OFF CACHE BOOL "Build the PCMSolver standalone executables")
set(ENABLE_TESTS OFF CACHE BOOL "Build the PCMSolver unit tests")
set(TEST_Fortran_API OFF CACHE BOOL "Build the PCMSolver Fortran 90 API tests")

# PCMSolver reaches for FindPythonInterp and FindBoost, both slated for removal.
# The default policy is what it wants; say so explicitly so its configure does
# not spray CMP0148/CMP0167 author warnings over a MADNESS build the user
# cannot act on.
set(CMAKE_POLICY_DEFAULT_CMP0148 OLD)
set(CMAKE_POLICY_DEFAULT_CMP0167 OLD)

FetchContent_Declare(
    pcmsolver
    GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver.git
    GIT_TAG        ${MADNESS_TRACKED_PCMSOLVER_TAG}
    GIT_SHALLOW    TRUE
    PATCH_COMMAND  ${CMAKE_COMMAND} -DPCMSOLVER_SOURCE_DIR=<SOURCE_DIR>
                   -P ${PROJECT_SOURCE_DIR}/cmake/patches/pcmsolver-v1.3.0.cmake
)
FetchContent_MakeAvailable(pcmsolver)

unset(CMAKE_POLICY_DEFAULT_CMP0148)
unset(CMAKE_POLICY_DEFAULT_CMP0167)

if (NOT TARGET pcm-static)
  message(FATAL_ERROR "FindOrFetchPCMSolver: the fetched PCMSolver did not define pcm-static")
endif ()

# VersionInfo.hpp is generated by PCMSolver's version target, which its sources
# include but do not depend on -- an upstream race that only shows up under a
# parallel build. Wire it up.
add_dependencies(pcm-objlib pcmsolver-update-version)

# Hand the *installed* pcm-static the Fortran runtime, minus whatever the C++
# driver links anyway. Inside this build tree it is not needed -- Fortran being
# an enabled language, CMake sees Fortran in pcm-static's link closure and adds
# the runtime itself, and adding it here too only earns a "duplicate libraries"
# warning on every executable link. install(EXPORT) does record the languages
# (IMPORTED_LINK_INTERFACE_LANGUAGES), but a consumer acts on those only if it
# has Fortran enabled, which a pure C++ consumer of MADNESS has no reason to do.
set(_pcm_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
if (CMAKE_CXX_IMPLICIT_LINK_LIBRARIES)
  list(REMOVE_ITEM _pcm_fortran_libs ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
endif ()
if (_pcm_fortran_libs)
  set_property(TARGET pcm-static APPEND PROPERTY
      INTERFACE_LINK_LIBRARIES "$<INSTALL_INTERFACE:${_pcm_fortran_libs}>")
  set_property(TARGET pcm-static APPEND PROPERTY
      INTERFACE_LINK_DIRECTORIES
      "$<INSTALL_INTERFACE:${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}>")
endif ()
unset(_pcm_fortran_libs)

add_library(PCMSolver::pcm ALIAS pcm-static)

# PCMSolver only arranges its headers under a PCMSolver/ directory as part of
# install(); its build-tree include interface points at the flat source dirs
# instead. chem/pcm.h asks for <PCMSolver/pcmsolver.h>, so stage the three
# headers that entails. All three exist by now: the first two are sources, and
# PCMSolverExport.h is written by generate_export_header() at configure time.
set(_pcm_staged_include "${PROJECT_BINARY_DIR}/external/pcmsolver-include")
file(COPY "${pcmsolver_SOURCE_DIR}/api/pcmsolver.h"
          "${pcmsolver_SOURCE_DIR}/api/PCMInput.h"
          "${pcmsolver_BINARY_DIR}/include/PCMSolverExport.h"
     DESTINATION "${_pcm_staged_include}/PCMSolver")

# ... and point the target at that instead of at the flat source dirs, which
# would otherwise be added PUBLIC to every MADNESS target that links PCM and
# put PCMSolver's own src/ and include/ trees on their header search path.
# The install interface is already correct and is preserved verbatim.
set_property(TARGET pcm-static PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    "$<BUILD_INTERFACE:${_pcm_staged_include}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

# Export PCMSolver's targets from the build tree as well as the install tree
# (it only does the latter itself), so that MADNESS's own
# export(EXPORT madness ...) can name a target that pcm-static resolves to.
export(EXPORT PCMSolverTargets-static
       NAMESPACE PCMSolver::
       FILE "${PROJECT_BINARY_DIR}/PCMSolver-targets.cmake")

# find_path/find_library from a failed FindPCM run leave -NOTFOUND cache
# entries behind; drop them rather than shadowing them with normal variables.
unset(PCM_INCLUDE_DIRS CACHE)
unset(PCM_LIBRARIES CACHE)

set(PCM_FOUND TRUE)
set(PCM_LIBRARIES PCMSolver::pcm)
set(PCM_INCLUDE_DIRS "")  # carried by the target, in both build and install trees
set(MADNESS_PCM_FETCHED TRUE)
unset(_pcm_staged_include)

message(STATUS "PCMSolver not found; building ${MADNESS_TRACKED_PCMSOLVER_TAG} from source")

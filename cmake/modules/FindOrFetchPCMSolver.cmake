# - Fetch and build PCMSolver, the polarizable continuum model library
#
# Included by external/pcm.cmake only after both discovery paths (the
# PCMSolverConfig package config, then the FindPCM module) have come up empty.
#
# On success:
#   PCMSolver::pcm      - an ALIAS for the in-tree pcm-static target, carrying
#                         the include directories for both build and install tree
#   MADNESS_PCM_FETCHED - TRUE, so the rest of the harness can tell this path
#                         apart from a PCMSolver found on the system
# external/pcm.cmake turns that target into PCM_FOUND / PCM_LIBRARIES /
# PCM_INCLUDE_DIRS for the rest of the build.
#
# Boost headers, PCMSolver's only non-toolchain dependency, are fetched too when
# the host has none -- so on a machine with a Fortran compiler and zlib this
# needs nothing preinstalled.
#
# On a missing prerequisite that cannot be fetched (a Fortran compiler, zlib)
# this sets MADNESS_PCM_UNAVAILABLE_REASON and returns without defining the
# target: PCM is an optional feature, and a machine without gfortran should
# still get a working MADNESS out of the default configure. external/pcm.cmake
# turns that reason into the one warning the user sees.

if (TARGET PCMSolver::pcm)
  return()
endif ()

cmake_minimum_required(VERSION 3.18.0)  # for file(ARCHIVE_EXTRACT)
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

# Boost headers are a hard requirement of PCMSolver's Green's functions
# (SphericalDiffuse.cpp pulls in boost/numeric/odeint.hpp; Getkw wants
# boost/any.hpp). No compiled Boost library is needed -- PCMSolver sets
# BOOST_COMPONENTS_REQUIRED to the empty string.
#
# Probed by header rather than find_package(Boost) to stay clear of CMP0167 --
# MADNESS's own Boost hook is off by default (ENABLE_BOOST), so there is
# usually no Boost target to reuse here. Coming up empty is not fatal: the
# headers are fetched below, once the prerequisites that cannot be fetched have
# been cleared.
find_path(PCM_BOOST_INCLUDE_DIR NAMES boost/version.hpp)
mark_as_advanced(PCM_BOOST_INCLUDE_DIR)

find_package(ZLIB QUIET)
if (NOT ZLIB_FOUND)
  list(APPEND _pcm_missing "zlib")
endif ()

if (_pcm_missing)
  # Hand the reason back rather than reporting it here: external/pcm.cmake
  # issues one warning for the whole search, so a machine that simply cannot
  # build PCMSolver says why in the same breath as "you asked for PCM and are
  # not getting it", instead of dribbling a STATUS line into the configure log
  # a few hundred lines above the summary.
  list(JOIN _pcm_missing ", " _pcm_missing)
  set(MADNESS_PCM_UNAVAILABLE_REASON
      "no installed PCMSolver was found, and building one here would need: ${_pcm_missing}")
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

# Nothing installed: fetch the headers. Deliberately after the checks above, so
# a host that was going to fail on a missing Fortran compiler fails before
# spending a download on it, and after enable_language(Fortran), which is the
# last thing that can still go wrong cheaply.
#
# Done with file(DOWNLOAD) + file(ARCHIVE_EXTRACT) rather than FetchContent,
# which the other FindOrFetch* modules use, because Boost is wanted here as a
# header tree and not as a project: there is nothing to add_subdirectory(), and
# the primitives allow extracting *only* the header tree. Measured on the b2
# tarball: extracting everything takes 71 s and deleting the 310 MB that are
# never read takes another 15 s, against 34 s to extract boost/ alone.
if (NOT PCM_BOOST_INCLUDE_DIR)
  if (NOT MADNESS_TRACKED_BOOST_VERSION OR NOT MADNESS_TRACKED_BOOST_URL_HASH)
    message(FATAL_ERROR "MADNESS_TRACKED_BOOST_{VERSION,URL_HASH} are empty; "
                        "external/versions.cmake must set them")
  endif ()

  set(_pcm_boost_stage "${PROJECT_BINARY_DIR}/external/boost")
  set(_pcm_boost_root "${_pcm_boost_stage}/boost-${MADNESS_TRACKED_BOOST_VERSION}")

  # The extracted headers are their own stamp: no marker file to fall out of
  # step with the tree it describes, and a reconfigure costs one EXISTS check.
  if (NOT EXISTS "${_pcm_boost_root}/boost/version.hpp")
    set(_pcm_boost_tarball "${_pcm_boost_stage}/boost-${MADNESS_TRACKED_BOOST_VERSION}.tar.xz")
    message(STATUS "Boost headers not found; fetching ${MADNESS_TRACKED_BOOST_VERSION} for the PCMSolver build")
    file(DOWNLOAD
         "https://github.com/boostorg/boost/releases/download/boost-${MADNESS_TRACKED_BOOST_VERSION}/boost-${MADNESS_TRACKED_BOOST_VERSION}-b2-nodocs.tar.xz"
         "${_pcm_boost_tarball}"
         EXPECTED_HASH ${MADNESS_TRACKED_BOOST_URL_HASH}
         STATUS _pcm_boost_dl)
    list(GET _pcm_boost_dl 1 _pcm_boost_dl_msg)
    list(GET _pcm_boost_dl 0 _pcm_boost_dl_code)
    if (NOT _pcm_boost_dl_code EQUAL 0)
      # Same policy as a missing Fortran compiler: an optional feature that
      # cannot be assembled here is a warning and no PCM, not a failed
      # configure. Drop the partial file so a later run re-downloads.
      file(REMOVE "${_pcm_boost_tarball}")
      set(MADNESS_PCM_UNAVAILABLE_REASON
          "no installed PCMSolver or Boost headers were found, and the Boost "
          "${MADNESS_TRACKED_BOOST_VERSION} headers could not be downloaded (${_pcm_boost_dl_msg})")
      return()
    endif ()

    # PATTERNS keeps the 276 MB of libs/ and 31 MB of tools/ in the archive.
    # The member paths carry the boost-<version>/ prefix, which is why
    # _pcm_boost_root sits one level below the staging directory.
    file(ARCHIVE_EXTRACT
         INPUT "${_pcm_boost_tarball}"
         DESTINATION "${_pcm_boost_stage}"
         PATTERNS "boost-${MADNESS_TRACKED_BOOST_VERSION}/boost/*")
    if (NOT EXISTS "${_pcm_boost_root}/boost/version.hpp")
      message(FATAL_ERROR
          "FindOrFetchPCMSolver: the fetched Boost ${MADNESS_TRACKED_BOOST_VERSION} "
          "yielded no boost/version.hpp -- the release artifact's layout has "
          "changed; revisit the URL and PATTERNS in this file.")
    endif ()
    # 51 MB that is only needed again if the headers go away, in which case the
    # download is a second and a half.
    file(REMOVE "${_pcm_boost_tarball}")
  endif ()

  set(PCM_BOOST_INCLUDE_DIR "${_pcm_boost_root}" CACHE PATH
      "Boost include directory used by the PCMSolver source build" FORCE)
  unset(_pcm_boost_stage)
  unset(_pcm_boost_root)
  unset(_pcm_boost_tarball)
  unset(_pcm_boost_dl)
  unset(_pcm_boost_dl_msg)
  unset(_pcm_boost_dl_code)
endif ()

# Hand PCMSolver the headers through Boost_INCLUDE_DIR rather than
# BOOST_INCLUDEDIR: its cmake/downloaded/autocmake_boost.cmake overwrites the
# latter with `set(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIR})` before calling
# find_package(Boost) -- so the cache entry is the only hint that survives.
#
# This is load-bearing well beyond the no-Boost case. Left to itself, that same
# module treats a failed find_package(Boost) as a cue to build its own: it
# downloads boost_1_54_0.zip from a 2013 SourceForge URL and unpacks 116 MB of
# 2013-era headers into the build tree, silently, on every host where CMake
# does not find Boost on the default search path. Pinning Boost_INCLUDE_DIR
# keeps that path unreachable, which is also what makes hunk 2 of
# cmake/patches/pcmsolver-v1.3.0.cmake matter -- PCMSolver now compiles against
# a modern Boost, whose odeint needs C++14.
set(Boost_INCLUDE_DIR "${PCM_BOOST_INCLUDE_DIR}" CACHE PATH
    "Boost include directory")

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

# PCM_FOUND / PCM_LIBRARIES / PCM_INCLUDE_DIRS are normalized by
# external/pcm.cmake once this returns, off the PCMSolver::pcm target -- which
# includes clearing the -NOTFOUND cache entries a failed FindPCM run leaves
# behind. All this has to do is announce the target.
set(MADNESS_PCM_FETCHED TRUE)
unset(_pcm_staged_include)

message(STATUS "PCMSolver not found; building ${MADNESS_TRACKED_PCMSOLVER_TAG} from source")

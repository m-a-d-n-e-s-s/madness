# Patch script for the PCMSolver source tree fetched by FindOrFetchPCMSolver.cmake.
#
# Run as `cmake -DPCMSOLVER_SOURCE_DIR=<dir> -P pcmsolver-v1.3.0.cmake`, from
# FetchContent's PATCH_COMMAND. It is idempotent: a hunk already applied is
# skipped, and a hunk whose context is missing is a hard error -- so bumping
# MADNESS_TRACKED_PCMSOLVER_TAG past the point where upstream moves any of this
# fails loudly instead of silently building something untested.
#
# PCMSolver v1.3.0 is from 2020 and receives no maintenance; every hunk here is
# a toolchain-compatibility fix, not a change of behaviour.

cmake_minimum_required(VERSION 3.12.0)

if (NOT PCMSOLVER_SOURCE_DIR)
  message(FATAL_ERROR "PCMSOLVER_SOURCE_DIR must be set")
endif ()

# Replace `old` with `new` in `relpath`, exactly once.
function(pcm_patch relpath old new what)
  set(_file "${PCMSOLVER_SOURCE_DIR}/${relpath}")
  if (NOT EXISTS "${_file}")
    message(FATAL_ERROR "pcmsolver patch: ${relpath} does not exist -- "
                        "the fetched PCMSolver is not the pinned v1.3.0")
  endif ()
  file(READ "${_file}" _contents)
  string(FIND "${_contents}" "${new}" _already)
  if (NOT _already EQUAL -1)
    return()  # already patched
  endif ()
  string(FIND "${_contents}" "${old}" _found)
  if (_found EQUAL -1)
    message(FATAL_ERROR "pcmsolver patch: could not apply '${what}' to ${relpath} -- "
                        "upstream text has changed; revisit this patch script")
  endif ()
  string(REPLACE "${old}" "${new}" _contents "${_contents}")
  file(WRITE "${_file}" "${_contents}")
endfunction()

# 1. CMake >= 4.0 refuses a project whose cmake_minimum_required is below 3.5.
pcm_patch(CMakeLists.txt
    "cmake_minimum_required(VERSION 3.3 FATAL_ERROR)"
    "cmake_minimum_required(VERSION 3.5 FATAL_ERROR)"
    "raise cmake_minimum_required to 3.5")

# 2. Boost >= 1.87 requires C++14 (boost::numeric::odeint uses generic lambdas),
#    and SphericalDiffuse.cpp pulls odeint in.
pcm_patch(cmake/custom/compilers/CXXFlags.cmake
    "set(CMAKE_CXX_STANDARD 11)"
    "set(CMAKE_CXX_STANDARD 14)"
    "build PCMSolver as C++14")

# 3. The bundled Eigen 3.3.2 does not compile with clang >= 11: `trt` is a
#    Transpose<TranspositionsBase<...>>, which has no derived(). Fixed upstream
#    in Eigen 3.3.5; applied here because PCMSolver vendors the header.
pcm_patch(external/eigen3/include/eigen3/Eigen/src/Core/Transpositions.h
    "Product<OtherDerived, Transpose, AliasFreeProduct>(matrix.derived(), trt.derived())"
    "Product<OtherDerived, Transpose, AliasFreeProduct>(matrix.derived(), trt)"
    "Eigen 3.3.2 Transpositions.h clang fix")

# 4. `update_version` is far too generic a target name to add to a build tree
#    that MADNESS may itself be a subproject of (see the CMake harness notes in
#    CLAUDE.md). Nothing else in PCMSolver refers to it.
pcm_patch(cmake/custom/pcmsolver.cmake
    "add_custom_target(update_version"
    "add_custom_target(pcmsolver-update-version"
    "rename the update_version target")

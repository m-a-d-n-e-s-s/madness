if(ENABLE_PCM)

  # PCMSolver >= 1.2 installs a CMake package config exporting PCMSolver::pcm,
  # which carries the transitive link deps a static build needs (zlib, the
  # Fortran runtime) -- worth trying before FindPCM, which only knows how to
  # look for one header and one library and leaves those deps to luck.
  find_package(PCMSolver CONFIG QUIET)

  # QUIET because coming up empty here is not yet news -- the fetch below is
  # the last word on whether MADNESS ends up with PCM.
  if(NOT TARGET PCMSolver::pcm)
    find_package(PCM COMPONENTS pcm QUIET)
  endif()

  # Nothing on the system: build our own copy, following the FindOrFetch*
  # pattern the other optional dependencies use.
  if(NOT TARGET PCMSolver::pcm AND NOT PCM_FOUND)
    if(MADNESS_FETCH_PCMSOLVER)
      include(${PROJECT_SOURCE_DIR}/cmake/modules/FindOrFetchPCMSolver.cmake)
    else()
      message(STATUS "PCMSolver not found and -DMADNESS_FETCH_PCMSOLVER=OFF; the `pcm` "
                     "keyword will be unavailable. Install PCMSolver and point "
                     "-DPCM_ROOT_DIR at the prefix, or configure -DENABLE_PCM=OFF "
                     "to stop looking.")
    endif()
  endif()

  # Whether the CONFIG flavour is the one we ended up USING -- the only thing
  # madness-config.cmake should act on, and false for both the FindPCM path
  # (baked into the exported target as an absolute library path) and the
  # fetched path (exported by PCMSolver's own install rules alongside us).
  if(TARGET PCMSolver::pcm AND NOT MADNESS_PCM_FETCHED)
    set(MADNESS_PCM_VIA_CONFIG ON)
  endif()

  # Normalize the three flavours onto the variables the rest of the build reads
  # (PCM_FOUND / PCM_LIBRARIES / PCM_INCLUDE_DIRS). For both target-based
  # flavours the include directories ride on the target, so PCM_INCLUDE_DIRS
  # stays empty and only the FindPCM path fills it in.
  #
  # Dropping FindPCM's cache entries is load-bearing, not tidiness: they are
  # find_path/find_library results, so a configure that took the FindPCM path
  # and failed leaves PCM_INCLUDE_DIRS-NOTFOUND behind. Reconfigure that build
  # tree after making a PCMSolver visible -- or after allowing the fetch -- and
  # find_package(PCM) is skipped, nothing rewrites the entry, and
  # src/madness{,/chem}/CMakeLists.txt feed a NOTFOUND to
  # target_include_directories: "variables are used in this project, but they
  # are set to NOTFOUND", generation aborted.
  if(TARGET PCMSolver::pcm)
    unset(PCM_INCLUDE_DIRS CACHE)
    unset(PCM_LIBRARIES CACHE)
    set(PCM_FOUND TRUE)
    set(PCM_LIBRARIES PCMSolver::pcm)
    set(PCM_INCLUDE_DIRS "")
  endif()

  # Set the output variables
  if(PCM_FOUND)
    set(MADNESS_HAS_PCM 1)
  endif()

endif()

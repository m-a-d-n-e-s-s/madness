if(ENABLE_PCM)

  # PCMSolver >= 1.2 installs a CMake package config exporting PCMSolver::pcm,
  # which carries the transitive link deps a static build needs (zlib, the
  # Fortran runtime) -- worth trying before FindPCM, which only knows how to
  # look for one header and one library and leaves those deps to luck.
  #
  # chem/pcm.cc needs the v1.3.0 API (pcmsolver_default_input(), declared in
  # PCMInput.h), so an older install must be passed over rather than accepted
  # and left to fail at compile time. PCMSolver's version file is
  # SameMinorVersion, so this admits 1.3.x only -- the release the fetch below
  # pins and patches, and the last one upstream made.
  find_package(PCMSolver 1.3 CONFIG QUIET)
  if(NOT TARGET PCMSolver::pcm AND PCMSolver_CONSIDERED_VERSIONS)
    message(STATUS "Ignoring PCMSolver ${PCMSolver_CONSIDERED_VERSIONS} "
        "(${PCMSolver_CONSIDERED_CONFIGS}): MADNESS needs 1.3")
  endif()

  # QUIET because coming up empty here is not yet news -- the fetch below is
  # the last word on whether MADNESS ends up with PCM.
  if(NOT TARGET PCMSolver::pcm)
    find_package(PCM COMPONENTS pcm QUIET)

    # FindPCM knows no version, so check for the API itself. Rejecting has to
    # drop its find_path/find_library cache entries too, or the next configure
    # skips the search and resurrects the old copy.
    if(PCM_FOUND)
      file(STRINGS "${PCM_INCLUDE_DIRS}/PCMSolver/PCMInput.h" _pcm_has_default_input
          REGEX "pcmsolver_default_input")
      if(NOT _pcm_has_default_input)
        message(STATUS "Ignoring PCMSolver at ${PCM_INCLUDE_DIRS}: its PCMInput.h "
            "lacks pcmsolver_default_input(), so it predates the 1.3 API MADNESS needs")
        set(PCM_FOUND FALSE)
        unset(PCM_INCLUDE_DIRS CACHE)
        unset(PCM_LIBRARIES CACHE)
      endif()
      unset(_pcm_has_default_input)
    endif()
  endif()

  # Nothing on the system: build our own copy, following the FindOrFetch*
  # pattern the other optional dependencies use.
  if(NOT TARGET PCMSolver::pcm AND NOT PCM_FOUND)
    if(MADNESS_FETCH_PCMSOLVER)
      # Sets MADNESS_PCM_UNAVAILABLE_REASON instead of defining the target if a
      # prerequisite of the source build is missing.
      include(${PROJECT_SOURCE_DIR}/cmake/modules/FindOrFetchPCMSolver.cmake)
    else()
      set(MADNESS_PCM_UNAVAILABLE_REASON
          "no installed PCMSolver was found, and -DMADNESS_FETCH_PCMSOLVER=OFF forbids building one")
    endif()
  endif()

  # Whether the CONFIG flavour is the one we ended up USING -- the only thing
  # madness-config.cmake should act on, and false for both the FindPCM path
  # (baked into the exported target as an absolute library path) and the
  # fetched path (exported by PCMSolver's own install rules alongside us).
  #
  # The PCMSolver_CONFIG test is not redundant with the target test: the target
  # can exist without our find_package having produced it. A parent project
  # that consumes MADNESS via add_subdirectory()/FetchContent may define
  # PCMSolver::pcm itself, in which case the target is here, PCMSolver_CONFIG
  # is empty, and claiming VIA_CONFIG would have madness-config.cmake hand
  # find_dependency() an empty PATHS with NO_DEFAULT_PATH -- a search of
  # nowhere, failing the consumer's configure even when a perfectly good
  # PCMSolver is installed on the default search path.
  if(TARGET PCMSolver::pcm AND NOT MADNESS_PCM_FETCHED AND PCMSolver_CONFIG)
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
  else()
    # ENABLE_PCM defaults to OFF, so reaching here means PCM was asked for by
    # name and is not being delivered -- a warning, not a STATUS line, because
    # the build otherwise succeeds and the loss only surfaces much later, when a
    # deck asking for `pcm` hits the MADNESS_EXCEPTION in the stub half of
    # chem/pcm.cc.
    if(NOT MADNESS_PCM_UNAVAILABLE_REASON)
      set(MADNESS_PCM_UNAVAILABLE_REASON "no installed PCMSolver was found")
    endif()
    message(WARNING
        "ENABLE_PCM is ON but PCM support could not be configured: "
        "${MADNESS_PCM_UNAVAILABLE_REASON}. MADNESS will build without it, and any "
        "calculation requesting the `pcm` solvation model will abort at runtime. "
        "Either point -DPCM_ROOT_DIR at an installed PCMSolver prefix, or install "
        "what the source build is missing, or configure -DENABLE_PCM=OFF to stop "
        "looking.")
  endif()

  # Report PCMSolver to FeatureSummary once, as a package, by the outcome that
  # actually matters -- rather than as a MADNESS feature, which would file an
  # unavailable *dependency* under "the following features have been disabled"
  # next to switches like GENTENSOR that are nobody's package.
  #
  # Two mechanics force the hand here. FeatureSummary skips any package whose
  # _CMAKE_<pkg>_QUIET global property is set, which find_package() sets for
  # every QUIET call -- so both probes above are invisible in the package lists
  # as they stand. And that property is keyed on the *name*, so it suppresses a
  # same-named add_feature_info() entry along with the package (measured: a
  # QUIET find_package(Foo) plus add_feature_info(Foo ...) is reported in
  # neither list). set_package_properties() does not clear it.
  #
  # So: drop both probe names -- which of the three paths delivered PCM is an
  # implementation detail of the search, not something to report -- and register
  # one un-QUIETed PCMSolver entry on the side the verdict belongs on.
  set(_pcm_probes PCMSolver PCM)
  if(NOT ENABLE_BOOST)
    # The vendored PCMSolver build resolves Boost too -- its own
    # find_package(Boost), fed the headers FindOrFetchPCMSolver settled on --
    # and that lands in the same global lists. MADNESS itself does not depend
    # on Boost (ENABLE_BOOST is off), so reporting it would advertise a
    # dependency the project does not have. When ENABLE_BOOST *is* on,
    # external/boost.cmake ran earlier and that entry is the user's: leave it.
    list(APPEND _pcm_probes Boost)
  endif()
  foreach(_pcm_probe IN LISTS _pcm_probes)
    foreach(_pcm_list PACKAGES_FOUND PACKAGES_NOT_FOUND)
      get_property(_pcm_pkgs GLOBAL PROPERTY ${_pcm_list})
      if(_pcm_pkgs)
        list(REMOVE_ITEM _pcm_pkgs ${_pcm_probe})
        set_property(GLOBAL PROPERTY ${_pcm_list} "${_pcm_pkgs}")
      endif()
    endforeach()
  endforeach()
  unset(_pcm_probes)
  unset(_pcm_probe)
  unset(_pcm_list)
  unset(_pcm_pkgs)

  set_property(GLOBAL PROPERTY _CMAKE_PCMSolver_QUIET FALSE)
  set_package_properties(PCMSolver PROPERTIES
      TYPE OPTIONAL
      DESCRIPTION "polarizable continuum model of solvation"
      PURPOSE "provides the `pcm` solvation model"
      URL "https://pcmsolver.readthedocs.io/")
  if(PCM_FOUND)
    set_property(GLOBAL APPEND PROPERTY PACKAGES_FOUND PCMSolver)
  else()
    set_property(GLOBAL APPEND PROPERTY PACKAGES_NOT_FOUND PCMSolver)
  endif()

endif()

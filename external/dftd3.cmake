if(ENABLE_DFTD3)

  # A source build of simple-dftd3 installs a CMake package config exporting
  # s-dftd3::s-dftd3, which carries the transitive link deps a static build needs
  # (gfortran, blas) -- worth trying first even though our own use is one header
  # and one library. Older conda-forge builds ship pkg-config only, which is what
  # FindDFTD3 handles.
  #
  # Newer conda-forge builds (1.5.0) do ship a config, but a broken one: it
  # find_dependency()s mctc-lib, whose own config reports NOT FOUND because
  # `tomlf` is missing from the package. That prints a CMake Warning and leaves
  # the target undefined -- harmless, since FindDFTD3 then picks the library up,
  # but do not go looking for the cause in MADNESS. QUIET does not suppress a
  # message(WARNING) raised inside someone else's config file.
  find_package(s-dftd3 CONFIG QUIET)
  if(TARGET s-dftd3::s-dftd3 AND DEFINED s-dftd3_DIR)
    # Record the config used so the installed madness-config.cmake can re-find the dependency.
    set(s-dftd3_CONFIG "${s-dftd3_DIR}/s-dftd3Config.cmake" CACHE INTERNAL "s-dftd3 package config used by MADNESS")
  endif()

  if(NOT TARGET s-dftd3::s-dftd3)
    find_package(DFTD3)
  endif()

  # Set the output variables
  if(TARGET s-dftd3::s-dftd3 OR DFTD3_FOUND)
    set(MADNESS_HAS_DFTD3 1)

    # Resolve the header to an absolute path so chem/dispersion.cc can #include
    # it without an -I -- see the comment there and in chem/CMakeLists.txt.
    if(TARGET s-dftd3::s-dftd3)
      get_target_property(_sdftd3_incdirs s-dftd3::s-dftd3 INTERFACE_INCLUDE_DIRECTORIES)
    else()
      set(_sdftd3_incdirs ${DFTD3_INCLUDE_DIRS})
    endif()
    find_file(DFTD3_HEADER NAMES s-dftd3.h HINTS ${_sdftd3_incdirs})
    mark_as_advanced(DFTD3_HEADER)
    unset(_sdftd3_incdirs)
  endif()

endif()

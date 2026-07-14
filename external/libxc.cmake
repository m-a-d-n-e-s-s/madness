if(ENABLE_LIBXC)

  # Preferred: a proper CMake package config (libxc built with -DENABLE_CMAKE
  # ships lib/cmake/Libxc/LibxcConfig.cmake).
  find_package(Libxc CONFIG QUIET)

  # Fallback: most HPC `libxc` environment modules are autotools installs that
  # ship no CMake config — only headers + libs exposed via C_INCLUDE_PATH /
  # LIBRARY_PATH (and a pkg-config .pc). Discover those and synthesize the
  # Libxc::xc target so `-DENABLE_LIBXC=ON` + `module load libxc` just works.
  if(NOT TARGET Libxc::xc)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
      pkg_check_modules(PC_LIBXC QUIET libxc)
    endif()
    find_path(LIBXC_INCLUDE_DIR NAMES xc.h
      HINTS ${PC_LIBXC_INCLUDE_DIRS} ${LIBXC_ROOT} ENV LIBXC_ROOT
      PATH_SUFFIXES include
      PATHS ENV C_INCLUDE_PATH ENV CPATH)
    find_library(LIBXC_LIBRARY NAMES xc
      HINTS ${PC_LIBXC_LIBRARY_DIRS} ${LIBXC_ROOT} ENV LIBXC_ROOT
      PATH_SUFFIXES lib lib64
      PATHS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    if(LIBXC_INCLUDE_DIR AND LIBXC_LIBRARY)
      add_library(Libxc::xc UNKNOWN IMPORTED)
      set_target_properties(Libxc::xc PROPERTIES
        IMPORTED_LOCATION "${LIBXC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIR}")
      message(STATUS "Libxc: no CMake config; using module/pkg-config install "
                     "(lib=${LIBXC_LIBRARY}, include=${LIBXC_INCLUDE_DIR})")
    endif()
  endif()

  # Set the output variables
  if(TARGET Libxc::xc)
    set(MADNESS_HAS_LIBXC 1)
    message(STATUS "MADNESS_HAS_LIBXC=1 (Libxc::xc available)")
  else()
    message(STATUS "ENABLE_LIBXC=ON but no libxc found (CMake config, pkg-config, "
                   "or C_INCLUDE_PATH/LIBRARY_PATH) — falling back to LDA-only XC")
  endif()

endif()

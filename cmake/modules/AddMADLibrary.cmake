macro(add_mad_library _name _source_files _header_files _dep_mad_comp _include_dir)

  if (NOT ${_source_files}) # no sources given? use add_mad_hdr_library
    message (FATAL_ERROR "add_mad_library: no sources given; probably want add_mad_hdr_library instead")
  endif()

  add_library(MAD${_name}-obj OBJECT ${${_source_files}} ${${_header_files}})
  # make library target(s)
  # if building shared library, build static as well using the same objects
  add_library(MAD${_name} $<TARGET_OBJECTS:MAD${_name}-obj>)
  add_dependencies(madness-libraries MAD${_name})
  if(BUILD_SHARED_LIBS)
    if (NOT DEFINED CMAKE_POSITION_INDEPENDENT_CODE)
      set_target_properties(MAD${_name}-obj PROPERTIES POSITION_INDEPENDENT_CODE TRUE)  # this is the default anyway, but produce a warning just in case
      message(WARNING "building shared libraries, setting default for POSITION_INDEPENDENT_CODE to true (set CMAKE_POSITION_INDEPENDENT_CODE to change the default)")
    endif()
  endif(BUILD_SHARED_LIBS)

  # Pass the private MAD${_name} compile flags to MAD${_name}-obj
  target_compile_definitions(MAD${_name}-obj PRIVATE 
      $<TARGET_PROPERTY:MAD${_name},COMPILE_DEFINITIONS>)
  target_include_directories(MAD${_name}-obj PRIVATE 
      $<TARGET_PROPERTY:MAD${_name},INCLUDE_DIRECTORIES>)
  target_compile_options(MAD${_name}-obj PRIVATE
      $<TARGET_PROPERTY:MAD${_name},COMPILE_OPTIONS>)
  # PRIVATE: warning flags must not propagate to consumers of the installed MAD${_name}.
  target_link_libraries(MAD${_name}-obj PRIVATE madness_internal_warnings)

  # target-common setup
  add_custom_target(install-madness-${_name}
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${_name} -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing ${_name} library components"
      USES_TERMINAL)
  add_dependencies(install-madness-${_name} install-madness-common)
  add_dependencies(install-madness-libraries install-madness-${_name})
  foreach(_dep ${_dep_mad_comp})
    if(TARGET install-madness-${_dep})
      add_dependencies(install-madness-${_name} install-madness-${_dep})
    endif()
  endforeach()

  # configure each target
    set(targetname MAD${_name})

    target_include_directories(${targetname} PUBLIC
        $<INSTALL_INTERFACE:${MADNESS_INSTALL_INCLUDEDIR}>)
    set_target_properties(${targetname} PROPERTIES PUBLIC_HEADER "${${_header_files}}")

    # Add library to the list of installed components
    install(TARGETS ${targetname} EXPORT madness
      COMPONENT ${_name}
      PUBLIC_HEADER DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/${_include_dir}"
      LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")
  
    # Create a target to install the component
    add_dependencies(install-madness-${_name} ${targetname})

    set(LINK_FLAGS "")
    foreach(_dep ${_dep_mad_comp})
      # NB header-only components also live under the MAD${_dep} name
      set(deptargetname MAD${_dep})

      if(TARGET ${deptargetname})
        # Propagate the dependency's usage requirements through the link alone.
        #
        # Do NOT also copy the dependency's INTERFACE_{COMPILE_DEFINITIONS,INCLUDE_DIRECTORIES,
        # COMPILE_OPTIONS} into ${targetname}'s own properties. The link below already gives
        # both ${targetname}'s consumers and MAD${_name}-obj (which forwards
        # $<TARGET_PROPERTY:${targetname},INCLUDE_DIRECTORIES> & co. above) the dependency's
        # full, transitively-resolved usage requirements, so the copy adds nothing to the
        # build graph -- but it does leave a $<TARGET_PROPERTY:${deptargetname},...> genex
        # sitting in ${targetname}'s own property. Downstream tooling that reads those
        # properties raw (e.g. TiledArray's DetectMADNESSConfig.cmake, which feeds
        # MADworld's INTERFACE_INCLUDE_DIRECTORIES into a standalone try_compile that knows
        # nothing about MADmisc) then dies with 'Target "MADmisc" not found.'
        # Copying the *resolved* values instead is no better: get_target_property() does not
        # traverse the dependency's own link interface, so the MPI/json/PaRSEC/LAPACK
        # requirements that reach us through it would be silently dropped, and the snapshot
        # would capture only what the dependency happens to carry at macro-call time.
        # PUBLIC (not INTERFACE) even when ${deptargetname} is header-only: MAD${_name}-obj
        # compiles ${targetname}'s own sources against $<TARGET_PROPERTY:${targetname},
        # INCLUDE_DIRECTORIES> & co., and only a PUBLIC link puts the dependency's usage
        # requirements into those properties. An INTERFACE link would serve consumers but
        # leave our own translation units without the dependency's include dirs and macros.
        target_link_libraries(${targetname} PUBLIC ${deptargetname})

        # import LINK_FLAGS from dependent
        get_property(deptargetname_LINK_FLAGS_SET TARGET ${deptargetname} PROPERTY LINK_FLAGS SET)
        if (deptargetname_LINK_FLAGS_SET)
          get_property(deptargetname_LINK_FLAGS TARGET ${deptargetname} PROPERTY LINK_FLAGS)
          set(LINK_FLAGS "${LINK_FLAGS} ${deptargetname_LINK_FLAGS}")
        endif ()
        
      endif()
    endforeach(_dep ${_dep_mad_comp})
    set_target_properties(${targetname} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
    target_compile_features(${targetname} INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

endmacro()


macro(add_mad_hdr_library _name _header_files _dep_mad_comp _include_dir)

  message (STATUS "in add_mad_hdr_library(${_name})")

  # make INTERFACE library
  add_library(MAD${_name} INTERFACE)
  
  # Add target dependencies
  add_dependencies(madness-libraries MAD${_name})
  
  target_include_directories(MAD${_name} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:${MADNESS_INSTALL_INCLUDEDIR}>
  )
  
  # Add library to the list of installed components
  install(TARGETS MAD${_name} EXPORT madness
      COMPONENT ${_name})
  
  # Create a target to install the component
  add_custom_target(install-madness-${_name}
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${_name} -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing ${_name} library components"
      USES_TERMINAL)
  add_dependencies(install-madness-${_name} MAD${_name})
  add_dependencies(install-madness-libraries install-madness-${_name})

  foreach(_dep ${_dep_mad_comp})
    if(TARGET install-madness-${_dep})
      add_dependencies(install-madness-${_name} install-madness-${_dep})
    endif()
    # NB the dependency target is MAD${_dep}, not ${_dep}; see add_mad_library() for why
    # the link interface alone carries the dependency's usage requirements.
    # MAD${_name} is a plain INTERFACE library (no sources), so INTERFACE is the only
    # scope CMake accepts here.
    if(TARGET MAD${_dep})
      target_link_libraries(MAD${_name} INTERFACE MAD${_dep})
    endif()
  endforeach()
  
  target_compile_features(MAD${_name} INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")
  if (CMAKE_CXX_STANDARD GREATER_EQUAL 20)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      target_compile_options(MAD${_name} INTERFACE "-Wno-deprecated-volatile")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      target_compile_options(MAD${_name} INTERFACE "-Wno-volatile")
    endif()
  endif()

endmacro()

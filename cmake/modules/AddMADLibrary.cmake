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
      if (${_dep}_is_mad_hdr_lib)
        set(deptargetname MAD${_dep})
      else(${_dep}_is_mad_hdr_lib)
        set(deptargetname MAD${_dep})
      endif(${_dep}_is_mad_hdr_lib)

      if(TARGET ${deptargetname})
        target_compile_definitions(${targetname} PUBLIC
            $<TARGET_PROPERTY:${deptargetname},INTERFACE_COMPILE_DEFINITIONS>)
        target_include_directories(${targetname} PUBLIC
            $<TARGET_PROPERTY:${deptargetname},INTERFACE_INCLUDE_DIRECTORIES>)
        target_compile_options(${targetname} PUBLIC
            $<TARGET_PROPERTY:${deptargetname},INTERFACE_COMPILE_OPTIONS>)
        if (${_dep}_is_mad_hdr_lib)
          target_link_libraries(${targetname} INTERFACE ${_dep})
        else()
          target_link_libraries(${targetname} PUBLIC ${deptargetname})
        endif()

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
  add_dependencies(libraries-madness MAD${_name})
  
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
    if(TARGET ${_dep})
        target_compile_definitions(MAD${_name} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_DEFINITIONS>)
        target_include_directories(MAD${_name} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_INCLUDE_DIRECTORIES>)
        target_compile_options(MAD${_name} PUBLIC 
          $<TARGET_PROPERTY:${_dep},INTERFACE_COMPILE_OPTIONS>)
      if (${_dep}_is_mad_hdr_lib)
        target_link_libraries(MAD${_name} INTERFACE ${_dep})
      else()
        target_link_libraries(MAD${_name} PUBLIC ${_dep})
      endif()
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

  set(${_name}_is_mad_hdr_lib TRUE)
endmacro()

macro(add_mad_library _name _source_files _header_files _dep_mad_comp _include_dir)

  # Create the MADNESS library
  add_library(MAD${_name}-obj OBJECT ${${_source_files}} ${${_header_files}})
  add_library(MAD${_name} $<TARGET_OBJECTS:MAD${_name}-obj>)
  set_target_properties(MAD${_name} PROPERTIES PUBLIC_HEADER "${${_header_files}}")
  
  # Add target dependencies
  add_library(${_name} ALIAS MAD${_name})
  add_dependencies(libraries MAD${_name})
  
  # Add library to the list of installed components
  install(TARGETS MAD${_name} EXPORT madness
      COMPONENT ${_name}
      PUBLIC_HEADER DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/${_include_dir}"
      LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")
  
  # Create a target to install the component
  add_custom_target(install-${_name}
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${_name} -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing ${_name} library components")
  add_dependencies(install-${_name} MAD${_name})
  add_dependencies(install-libraries install-${_name})

  foreach(_dep ${_dep_mad_comp})
    if(TARGET install-${_dep})
      add_dependencies(install-${_name} install-${_dep})
    endif()
    if(TARGET MAD${_dep})
      target_link_libraries(MAD${_name} PUBLIC MAD${_dep})
      append_target_properties(MAD${_dep} MAD${_name} 
          "INTERFACE_INCLUDE_DIRECTORIES;INTERFACE_COMPILE_DEFINITIONS;INTERFACE_COMPILE_OPTIONS")
    elseif(TARGET ${_dep})
      target_link_libraries(MAD${_name} PUBLIC ${_dep})
      append_target_properties(${_dep} MAD${_name}
          "INTERFACE_INCLUDE_DIRECTORIES;INTERFACE_COMPILE_DEFINITIONS;INTERFACE_COMPILE_OPTIONS")
    endif()
  endforeach()
  
  append_target_properties(MAD${_name} MAD${_name}-obj 
      "INTERFACE_INCLUDE_DIRECTORIES;INTERFACE_COMPILE_DEFINITIONS;INTERFACE_COMPILE_OPTIONS")
  
  # Add compile and linker flags to library
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    target_link_libraries(MAD${_name} PUBLIC "-Wl,-no_pie")
  endif()
  
endmacro()
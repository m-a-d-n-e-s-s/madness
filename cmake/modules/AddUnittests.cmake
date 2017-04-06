macro(add_unittests _component _sources _libs)

  # Add targets and for world_unittests
  add_custom_target(${_component}_unittests)
  add_dependencies(unittests ${_component}_unittests)
    
  # Add a test that builds the unit tests
  add_test(build_${_component}_unittests
      "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${_component}_unittests)
  
  foreach(_source ${${_sources}})
    # Get the test name (the file name of the first source)
    string(REGEX MATCH "[A-Za-z_][A-Za-z0-9_]*\\.cc" _test_source "${_source}")
    string(REGEX MATCHALL "[A-Za-z0-9_\\.\\$<:>]+" _source_list "${_source}")
    get_filename_component(_test "${_test_source}" NAME_WE)
    
    # Create test executables
    add_executable(${_test} EXCLUDE_FROM_ALL ${_source_list})
    target_link_libraries(${_test} ${_libs})
    
    # Add the test and set dependencies
    add_test(NAME ${_component}-${_test} COMMAND ${_test})
    add_dependencies(${_component}_unittests ${_test})
    set_tests_properties(${_component}-${_test} 
        PROPERTIES DEPENDS build_${_component}_unittests)
 
  endforeach()

endmacro()
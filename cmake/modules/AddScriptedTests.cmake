# macro will build (i.e. copy the test script to the build directory) and execute the script

# example
# CMakeLists.txt: add_scripted_tests(nemo_test1.py nemo)


macro(add_scripted_tests _testcase_in _binary _labels)

  # convert test.py.in to test.py
#  get_filename_component(_testcase "${_testcase_in}" NAME_WLE)
  # possibly path/to/src/test.py.in and path/to/build/test.py
  set(_testcase ${_testcase_in})

  # Add targets and for scripted tests
  add_custom_target_subproject(madness ${_testcase}_${_binary}_scripted_tests)
  add_dependencies(scripted_tests-madness ${_testcase}_${_binary}_scripted_tests-madness)

#   Add a test that builds the binary
  add_test(madness/test/${_binary}/build
          "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${_binary})
  set_tests_properties(madness/test/${_binary}/build PROPERTIES DEPENDS ${_binary})

  # make sure that the build step has all labels
  set_property(TEST madness/test/${_binary}/build PROPERTY LABELS "${_labels}" APPEND)

#   Add a test that copies the test scripts and replaces the variable to the source directory
#   containing the reference json outputs
  if (0)
  add_test(madness/test/scripted_tests/${_binary}/${_testcase}/copy
          COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/${_testcase}
          ${CMAKE_CURRENT_BINARY_DIR}/${_testcase})
  set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/copy
          PROPERTIES DEPENDS ${_binary} LABELS "${_labels}")
  endif()
  #  copy the test scripts and replaces the variable to the source directory containing the reference json outputs
  set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})
  set(BINARY ${_binary})
  set(TESTCASE ${_testcase})
  configure_file(${_testcase} ${_testcase} @ONLY)

#  message(STATUS "testcase: " ${_testcase})
#  message(STATUS "binary:   " ${_binary})
#  message(STATUS "labels:   " ${_labels})
#  message(STATUS "sourcedir "  ${CMAKE_CURRENT_SOURCE_DIR})
#  message(STATUS "binarydir "  ${CMAKE_CURRENT_BINARY_DIR})

  # Add the tests (execution and result) and set dependencies
  add_test(NAME madness/test/scripted_tests/${_binary}/${_testcase}/run COMMAND ${_testcase} --reference_directory=${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/run
          PROPERTIES LABELS "${_labels}")

  # if "labels" is "verylong", potentially skip the test on small machines
    if ("${_labels}" MATCHES "verylong")
        set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/run
                PROPERTIES SKIP_RETURN_CODE 77)
        set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/run
                PROPERTIES SKIP_MESSAGE "Test skipped on small machines")
    endif()


endmacro()
# macro will build (i.e. copy the test script to the build directory) and execute the script

# example
# CMakeLists.txt: add_scripted_tests(nemo_test1.py nemo)


macro(add_scripted_tests _testcase _binary _labels)

  # Add targets and for scripted tests
  add_custom_target_subproject(madness ${_testcase}_scripted_tests)
  add_dependencies(scripted_tests-madness ${_testcase}_scripted_tests-madness)

#   Add a test that builds the binary
  add_test(madness/test/${_binary}/build
          "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${_binary})
  set_tests_properties(madness/test/${_binary}/build
          PROPERTIES DEPENDS ${_binary} LABELS scripted_tests)

#   Add a test that copies the test scripts
  add_test(madness/test/scripted_tests/${_binary}/${_testcase}/copy
          COMMAND ${CMAKE_COMMAND} -E copy
          ${CMAKE_CURRENT_SOURCE_DIR}/${_testcase}
          ${CMAKE_CURRENT_BINARY_DIR}/${_testcase})
  set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/copy
          PROPERTIES DEPENDS ${_binary} LABELS scripted_tests)

  message(STATUS "testcase: " ${_testcase})
  message(STATUS "binary:   " ${_binary})
  message(STATUS "labels:   " ${_labels})
  message(STATUS "sourcedir "  ${CMAKE_CURRENT_SOURCE_DIR})
  message(STATUS "binarydir "  ${CMAKE_CURRENT_BINARY_DIR})

  # Add the tests (execution and result) and set dependencies
  add_test(NAME madness/test/scripted_tests/${_binary}/${_testcase}/run COMMAND ${_testcase})
  set_tests_properties(madness/test/scripted_tests/${_binary}/${_testcase}/run
          PROPERTIES LABELS scripted_tests)

endmacro()
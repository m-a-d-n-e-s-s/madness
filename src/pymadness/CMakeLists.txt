if(NOT MADNESS_HAS_PYTHON)
  return()
endif()

pybind11_add_module(_pymadness
    src/_pymadness/module.cpp
    src/_pymadness/bind_world.cpp
    src/_pymadness/bind_defaults.cpp
    src/_pymadness/bind_function.cpp
    src/_pymadness/bind_operators.cpp
    src/_pymadness/bind_tensor.cpp
)

target_link_libraries(_pymadness PRIVATE madness)
target_compile_features(_pymadness PRIVATE cxx_std_17)
target_include_directories(_pymadness PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_BINARY_DIR}/src
)

# Install the compiled extension next to the Python package
install(TARGETS _pymadness
    LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}/pymadness"
    COMPONENT python)

# Install the pure-Python package
install(DIRECTORY src/pymadness/
    DESTINATION "${MADNESS_INSTALL_LIBDIR}/pymadness"
    COMPONENT python
    FILES_MATCHING PATTERN "*.py")

# --- Tests ---
if(BUILD_TESTING)
  # Discover the module output directory so tests can import it
  set(PYMADNESS_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Copy pure-python package next to the .so for testing (re-copies on every build)
  add_custom_command(TARGET _pymadness POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
              "${CMAKE_CURRENT_SOURCE_DIR}/src/pymadness"
              "${CMAKE_CURRENT_BINARY_DIR}/pymadness"
      COMMENT "Syncing pymadness Python package to build directory")

  add_test(NAME pymadness/test_basic
      COMMAND ${Python3_EXECUTABLE} -m pytest -xvs
              "${CMAKE_CURRENT_SOURCE_DIR}/tests/"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(pymadness/test_basic PROPERTIES
      ENVIRONMENT "PYTHONPATH=${PYMADNESS_MODULE_DIR}:${PYMADNESS_MODULE_DIR}/pymadness")
endif()

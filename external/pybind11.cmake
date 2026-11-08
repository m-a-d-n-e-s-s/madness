if(ENABLE_PYTHON)

  find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
  # Build against the Python3 found above. As a subproject pybind11 otherwise
  # defaults to its classic FindPythonLibsNew mode, which ignores Python3_* and
  # resolves `python3` from PATH on its own -- possibly a different minor
  # version (e.g. Homebrew's), yielding an extension that will not import.
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 2.11 CONFIG)

  if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found via CONFIG, trying FetchContent...")
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)
  endif()

  if(TARGET pybind11::module)
    set(MADNESS_HAS_PYTHON 1)
    message(STATUS "Python bindings enabled (pybind11 ${pybind11_VERSION})")
  else()
    message(WARNING "pybind11 targets not available, Python bindings disabled")
  endif()

endif()

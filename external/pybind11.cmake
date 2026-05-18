if(ENABLE_PYTHON)

  find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
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

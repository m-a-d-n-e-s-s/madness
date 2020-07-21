# Find MPI

if(ENABLE_MPI)

  # find threads first
  if (NOT TARGET Threads::Threads)
    message(FATAL_ERROR "include external/pthread.cmake BEFORE external/mpi.cmake")
  endif()

  # Try to find MPI
  find_package(MPI REQUIRED)
  cmake_minimum_required(VERSION 3.10) # FindMPI needs to provide MPI_<lang>_HEADER_DIR

  # Set the variables 
  if(MPI_C_FOUND)
    set(MPI_FOUND         ${MPI_C_FOUND})
    string(STRIP "${MPI_C_COMPILE_FLAGS}" MPI_COMPILE_FLAGS)
    set(MPI_INCLUDE_PATH  ${MPI_C_INCLUDE_PATH})
    string(STRIP "${MPI_C_LINK_FLAGS}" MPI_LINK_FLAGS)
    set(MPI_LIBRARIES     ${MPI_C_LIBRARIES})
    set(MPI_HEADER_DIR    ${MPI_C_HEADER_DIR})
  elseif(MPI_CXX_FOUND)
    set(MPI_FOUND         ${MPI_CXX_FOUND})
    string(STRIP "${MPI_CXX_COMPILE_FLAGS}" MPI_COMPILE_FLAGS)
    set(MPI_INCLUDE_PATH  ${MPI_CXX_INCLUDE_PATH})
    string(STRIP "${MPI_CXX_LINK_FLAGS}" MPI_LINK_FLAGS)
    set(MPI_LIBRARIES     ${MPI_CXX_LIBRARIES})
    set(MPI_HEADER_DIR    ${MPI_CXX_HEADER_DIR})
  else()
    message(FATAL_ERROR "No suitable MPI compiler was not found.")
  endif()
  # will hardwire to particular mpi.h to avoid the issues with CPATH defined by Intel MPI envvar script overriding mpi.h in system dirs
  if (EXISTS "${MPI_HEADER_DIR}")
    get_filename_component(MADNESS_MPI_HEADER "${MPI_HEADER_DIR}/mpi.h" ABSOLUTE CACHE)
  elseif (EXISTS "${MPI_INCLUDE_PATH}")
    get_filename_component(MADNESS_MPI_HEADER "${MPI_INCLUDE_PATH}/mpi.h" ABSOLUTE CACHE)
  else()
    set(MADNESS_MPI_HEADER "mpi.h" CACHE STRING "Path to the main MPI header")
  endif()
  message(STATUS "MPI main header: ${MADNESS_MPI_HEADER}")

  # filter out -pthread from COMPILE and LINK flags, use Threads::Threads instead
  # this is to avoid issues later consuming madness targets in codes using CUDA
  # see https://gitlab.kitware.com/cmake/cmake/merge_requests/2512
  string(REGEX REPLACE "[ ]?-pthread" "" MPI_COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  string(REGEX REPLACE "[ ]?-pthread" "" MPI_LINK_FLAGS "${MPI_LINK_FLAGS}")
  set(MPI_LIBRARIES ${MPI_LIBRARIES} Threads::Threads)

  # MPIEXEC was deprecated in CMake 3.10
  if (CMAKE_VERSION VERSION_LESS 3.10.0)
    set(MPIEXEC_EXECUTABLE ${MPIEXEC})
  endif()

else()
  # Disable MPI via config.h
  set(STUBOUTMPI 1)
endif()

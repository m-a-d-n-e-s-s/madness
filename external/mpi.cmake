# Find MPI

if(ENABLE_MPI)

  # find threads first
  if (NOT TARGET Threads::Threads)
    message(FATAL_ERROR "include external/pthread.cmake BEFORE external/mpi.cmake")
  endif()

  # Try to find MPI
  find_package(MPI REQUIRED)
  
  # Set the variables 
  if(MPI_C_FOUND)
    set(MPI_FOUND         ${MPI_C_FOUND})
    string(STRIP "${MPI_C_COMPILE_FLAGS}" MPI_COMPILE_FLAGS)
    set(MPI_INCLUDE_PATH  ${MPI_C_INCLUDE_PATH})
    string(STRIP "${MPI_C_LINK_FLAGS}" MPI_LINK_FLAGS)
    set(MPI_LIBRARIES     ${MPI_C_LIBRARIES})
  elseif(MPI_CXX_FOUND)
    set(MPI_FOUND         ${MPI_CXX_FOUND})
    string(STRIP "${MPI_CXX_COMPILE_FLAGS}" MPI_COMPILE_FLAGS)
    set(MPI_INCLUDE_PATH  ${MPI_CXX_INCLUDE_PATH})
    string(STRIP "${MPI_CXX_LINK_FLAGS}" MPI_LINK_FLAGS)
    set(MPI_LIBRARIES     ${MPI_CXX_LIBRARIES})
  else()
    message(FATAL_ERROR "No suitable MPI compiler was not found.")
  endif()

  # filter out -pthread from COMPILE and LINK flags, use Threads::Threads instead
  # this is to avoid issues later consuming madness targets in codes using CUDA
  # see https://gitlab.kitware.com/cmake/cmake/merge_requests/2512
  string(REGEX REPLACE "[ ]?-pthread" "" MPI_COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  string(REGEX REPLACE "[ ]?-pthread" "" MPI_LINK_FLAGS "${MPI_LINK_FLAGS}")
  set(MPI_LIBRARIES ${MPI_LIBRARIES} Threads::Threads)

else()
  # Disable MPI via config.h
  set(STUBOUTMPI 1)
endif()
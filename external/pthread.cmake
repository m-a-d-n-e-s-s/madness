# Find Threads
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

# Check that pthreads was found
if(NOT CMAKE_USE_PTHREADS_INIT)
  message(FATAL_ERROR "MADNESS requires pthreads support.")
endif()

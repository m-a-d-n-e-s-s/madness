# Find Threads
if (NOT TARGET Threads::Threads)
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  find_package(Threads REQUIRED)
endif()

# Check that pthreads was found
if(NOT TARGET Threads::Threads)
  message(FATAL_ERROR "MADNESS requires pthreads support.")
endif()

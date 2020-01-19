# Find Threads
if (NOT TARGET Threads::Threads)
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  find_package(Threads REQUIRED)
endif()

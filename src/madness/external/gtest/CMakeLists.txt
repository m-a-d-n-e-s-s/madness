# src/madness/external/gtest

if(BUILD_TESTING)
  add_library(MADgtest gtest-all.cc)
  target_include_directories(MADgtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(MADgtest Threads::Threads)
endif()

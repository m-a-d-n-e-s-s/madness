# src/madness/external/gtest

if(ENABLE_UNITTESTS)
  set(gtest_sources gtest-all.cc)
  set(gtest_headers )
  add_mad_library(gtest gtest_sources gtest_headers "Threads::Threads" "madness/external/gtest")
  target_include_directories(MADgtest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif()

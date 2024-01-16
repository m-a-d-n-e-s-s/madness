cmake_minimum_required(VERSION 3.16)
project(madness-tutorial)

find_package(madness REQUIRED OPTIONAL_COMPONENTS chem)

if (TARGET MADchem)
  add_executable(simple_hf simple_hf.cpp)
  target_link_libraries(simple_hf MADchem)
endif()

add_executable(test_runtime test_runtime.cpp)
target_link_libraries(test_runtime MADworld)

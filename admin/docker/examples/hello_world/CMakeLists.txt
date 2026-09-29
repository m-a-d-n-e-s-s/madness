cmake_minimum_required(VERSION 3.20)
project(HelloMADNESS CXX)

# MADNESS requires C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate installed MADNESS package
find_package(madness REQUIRED CONFIG)

# Define executable and link against the MADNESS imported target
add_executable(hello_madness hello.cc)
target_link_libraries(hello_madness PRIVATE madness)

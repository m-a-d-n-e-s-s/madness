# src/madness/external/gtest

if(ENABLE_UNITTESTS)
  add_library(MADgtest gtest-all.cc)
  target_include_directories(MADgtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  
  target_compile_options(MADgtest PUBLIC "${CXX11_FLAG}")
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    target_link_libraries(MADgtest PUBLIC "-Wl,-no_pie")
  endif()
  if(LIBUNWIND_FOUND)
    target_include_directories(MADgtest PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
    target_link_libraries(MADgtest PUBLIC ${LIBUNWIND_LIBRARIES})
  endif()
  if(GPERFTOOLS_FOUND)
    target_include_directories(MADgtest PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
    target_link_libraries(MADgtest PUBLIC ${GPERFTOOLS_LIBRARIES})
  endif()
endif()
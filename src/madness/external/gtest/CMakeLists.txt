# src/madness/external/gtest

if(ENABLE_UNITTESTS)
  add_library(MADgtest gtest-all.cc)
  target_include_directories(MADgtest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
  
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    target_link_libraries(MADgtest INTERFACE "-Wl,-no_pie")
  endif()
  if(LIBUNWIND_FOUND)
    target_include_directories(MADgtest PRIVATE ${LIBUNWIND_INCLUDE_DIRS})
    target_link_libraries(MADgtest PRIVATE ${LIBUNWIND_LIBRARIES})
  endif()
  if(GPERFTOOLS_FOUND)
    target_include_directories(MADgtest PRIVATE ${GPERFTOOLS_INCLUDE_DIRS})
    target_link_libraries(MADgtest PRIVATE ${GPERFTOOLS_LIBRARIES})
  endif()
endif()
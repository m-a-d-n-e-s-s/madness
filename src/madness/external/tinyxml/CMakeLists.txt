# src/madness/external/tinyxml

set(MADTINYXML_HEADERS tinystr.h tinyxml.h)
set(MADTINYXML_SOURCES tinystr.cc tinyxml.cc tinyxmlerror.cc tinyxmlparser.cc)

# Create the MADtinyxml library
add_mad_library(tinyxml MADTINYXML_SOURCES MADTINYXML_HEADERS "dummy" "madness/external/tinyxml/")

target_compile_options(MADtinyxml PUBLIC "${CXX11_FLAG}")
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  target_link_libraries(MADtinyxml PUBLIC "-Wl,-no_pie")
endif()
if(GPERFTOOLS_FOUND)
  target_include_directories(MADtinyxml PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(MADtinyxml PUBLIC ${GPERFTOOLS_LIBRARIES})
endif()
if(LIBUNWIND_FOUND)
  target_include_directories(MADtinyxml PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
  target_link_libraries(MADtinyxml PUBLIC ${LIBUNWIND_LIBRARIES})
endif()
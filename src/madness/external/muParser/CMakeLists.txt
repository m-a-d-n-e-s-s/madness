# src/madness/external/muParser

set(MADMUPARSER_HEADERS 
    muParserBase.h muParserBytecode.h muParserCallback.h muParserDef.h
    muParserDLL.h muParserError.h muParserFixes.h muParser.h muParserInt.h 
    muParserStack.h muParserTest.h muParserToken.h muParserTokenReader.h)
set(MADMUPARSER_SOURCES
    muParserBase.cpp muParserBytecode.cpp muParserCallback.cpp muParser.cpp
    muParserError.cpp muParserTokenReader.cpp muParserBase.h muParserBytecode.h
    muParserCallback.h muParserDef.h muParserDLL.h muParserError.h
    muParserFixes.h muParser.h muParserInt.h muParserStack.h muParserTest.h 
    muParserToken.h muParserTokenReader.h)

# Create the MADmuparser library
add_mad_library(muparser MADMUPARSER_SOURCES MADMUPARSER_HEADERS "dummy" "madness/external/muParser/")

target_compile_options(MADmuparser PUBLIC "${CXX11_FLAG}")
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  target_link_libraries(MADmuparser PUBLIC "-Wl,-no_pie")
endif()
if(GPERFTOOLS_FOUND)
  target_include_directories(MADmuparser PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(MADmuparser PUBLIC ${GPERFTOOLS_LIBRARIES})
endif()
if(LIBUNWIND_FOUND)
  target_include_directories(MADmuparser PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
  target_link_libraries(MADmuparser PUBLIC ${LIBUNWIND_LIBRARIES})
endif()

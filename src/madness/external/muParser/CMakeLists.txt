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
add_mad_library(muparser MADMUPARSER_SOURCES MADMUPARSER_HEADERS "dummy" "madness/external/muParser")

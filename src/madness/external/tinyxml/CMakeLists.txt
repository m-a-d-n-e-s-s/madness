# src/madness/external/tinyxml

set(MADTINYXML_HEADERS tinystr.h tinyxml.h)
set(MADTINYXML_SOURCES tinystr.cc tinyxml.cc tinyxmlerror.cc tinyxmlparser.cc)

# Create the MADtinyxml library
add_mad_library(tinyxml MADTINYXML_SOURCES MADTINYXML_HEADERS "dummy" "madness/external/tinyxml/")

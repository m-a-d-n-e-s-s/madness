# src/madness/external/nlohmann_json

set(MADNLOHMANN_JSON_HEADERS json.hpp)
set(MADNLOHMANN_JSON_SOURCES json.cc )
#set(MADTINYXML_SOURCES tinystr.cc tinyxml.cc tinyxmlerror.cc tinyxmlparser.cc)

# Create the MADtinyxml library
add_mad_library(nlohmann_json MADNLOHMANN_JSON_SOURCES MADNLOHMANN_JSON_HEADERS "dummy" "madness/external/nlohmann_json/")
#add_mad_library(tinyxml MADTINYXML_SOURCES MADTINYXML_HEADERS "dummy" "madness/external/tinyxml/")
#add_mad_hdr_library(nlohmann_json MADNLOHMANN_JSON_HEADERS "dummy" "madness/external/nlohmann_json/")

#target_compile_options(MADtinyxml PRIVATE -O0)
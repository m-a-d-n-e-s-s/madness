# src/madness/misc

set(MADMISC_HEADERS misc.h ran.h phandler.h interpolation_1d.h cfft.h)
set(MADMISC_SOURCES
    checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc)

# Create the MADmisc library
add_mad_library(misc MADMISC_SOURCES MADMISC_HEADERS "world" "madness/misc/")

if(ENABLE_UNITTESTS)
  
  # The list of unit test source files
  set(MISC_TEST_SOURCES interp3.cc)
  
  add_unittests(misc MISC_TEST_SOURCES "MADmisc;MADgtest")

endif()
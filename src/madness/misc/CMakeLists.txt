# src/madness/misc

set(MADMISC_HEADERS misc.h ran.h phandler.h interpolation_1d.h cfft.h info.h gnuplot.h)
set(MADMISC_SOURCES
    checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc info.cc unique_filename.cc)
# retrieve git metadata
include(GetGitMetadata)
vgkit_cmake_git_metadata()
# feed MADNESS_GIT_{REVISION,DESCRIPTION} to info.cc only to avoid recompiling everything
set_source_files_properties(
        info.cc
        PROPERTIES
          COMPILE_DEFINITIONS
          "MADNESS_GIT_REVISION=\"${MADNESS_GIT_REVISION}\";MADNESS_GIT_DESCRIPTION=\"${MADNESS_GIT_DESCRIPTION}\""
)

# Create the MADmisc library
add_mad_library(misc MADMISC_SOURCES MADMISC_HEADERS "world" "madness/misc/")

if(BUILD_TESTING)
  # The list of unit test source files
  # test_gnuplot.cc breaks the CI since gnuplot missing so disable 
  set(MISC_TEST_SOURCES interp3.cc )
  
  add_unittests(misc "${MISC_TEST_SOURCES}" "MADmisc;MADgtest" "unittests;short")

endif()

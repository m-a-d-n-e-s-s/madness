# src/madness/misc

set(MADMISC_HEADERS misc.h ran.h phandler.h interpolation_1d.h cfft.h)
set(MADMISC_SOURCES
    checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc)

# Create the MADmisc library
add_mad_library(misc MADMISC_SOURCES MADMISC_HEADERS "world" "madness/misc/")

# Create the compile-time git library: git_version.cmake will write the gitinfo.cc file,
add_custom_command(
    OUTPUT  "${PROJECT_BINARY_DIR}/src/madness/misc/git_version.h"
           SYMBOLIC "${PROJECT_BINARY_DIR}/src/madness/misc/_git_version.h"
    COMMAND ${CMAKE_COMMAND}
      -Dlocal_dir="${PROJECT_SOURCE_DIR}/cmake/"
      -Doutput_dir="${PROJECT_BINARY_DIR}/src/madness/misc"
      -Dsource_dir="${PROJECT_SOURCE_DIR}/src/madness/misc"
      -P "${PROJECT_SOURCE_DIR}/cmake/git_version.cmake"
 	COMMENT "generating the gitversion file"
)
add_custom_target_subproject(madness mgit
    DEPENDS "${PROJECT_BINARY_DIR}/src/madness/misc/git_version.h"
)

set(MADGIT_HEADERS gitinfo.h)
set(MADGIT_SOURCES gitinfo.cc )

add_mad_library(git MADGIT_SOURCES MADGIT_HEADERS "world" "madness/misc/")
add_dependencies(MADgit-obj mgit-madness)

if(ENABLE_UNITTESTS)
  
  # The list of unit test source files
  set(MISC_TEST_SOURCES interp3.cc)
  
  add_unittests(misc "${MISC_TEST_SOURCES}" "MADmisc;MADgtest")

endif()
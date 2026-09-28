# src/madness/misc

set(MADMISC_HEADERS misc.h ran.h cfft.h info.h gnuplot.h array_of_bools.h kahan_accumulator.h)
set(MADMISC_SOURCES
    checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc info.cc numerics.cc unique_filename.cc)
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

# Create the MADmisc library.
# MADmisc is a leaf: it must not include anything from world/tensor/mra/chem, so that
# MADworld can depend on it (see world/CMakeLists.txt).
add_mad_library(misc MADMISC_SOURCES MADMISC_HEADERS "" "madness/misc/")
target_link_libraries(MADmisc PUBLIC Threads::Threads)

# No unit tests registered here: interp3.cc moved to src/madness/mra along with
# interpolation_1d.h, and test_gnuplot.cc stays disabled (gnuplot missing in CI).

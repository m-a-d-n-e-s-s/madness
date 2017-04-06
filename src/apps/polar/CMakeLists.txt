# src/apps/polar

# Add source files for moldft_polar
set(MOLDFT_POLAR_SOURCES 
    moldft.cc lda.cc  molecule.cc molecularbasis.cc atomutil.cc corepotential.cc
    distpm.cc)
if(LIBXC_FOUND)
  list(APPEND MOLDFT_POLAR_SOURCES xcfunctional_libxc.cc)
else()
  list(APPEND MOLDFT_POLAR_SOURCES xcfunctional_ldaonly.cc)
endif()

# Create the moldft_polar executable
add_executable(moldft_polar ${MOLDFT_POLAR_SOURCES})
target_link_libraries(moldft_polar MADmra)
if(LIBXC_FOUND)
  target_link_libraries(moldft_polar libxc)
endif()


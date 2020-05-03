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
if(LIBXC_FOUND)
  add_mad_executable(moldft_polar ${MOLDFT_POLAR_SOURCES} "MADmra;libxc")
else(LIBXC_FOUND)
  add_mad_executable(moldft_polar ${MOLDFT_POLAR_SOURCES} MADmra)
endif(LIBXC_FOUND)


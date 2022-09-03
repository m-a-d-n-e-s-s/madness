# src/apps/molresponse

add_subdirectory(testing)

# Set the CHEM sources and header files
set(MOLRESPONSE_SOURCES
        density.cc
        property.cc
        TDDFT.cc
        global_functions.cc
        calc_runner.cc
        timer.cc
        basic_operators.cc
        Plot_VTK.cc
        iter_freq2.cc
        iterate_excited.cc
        Lambda_X.cc
        Theta_X.cc
        response_parameters.cpp
        iterate_gamma.cc
        )
set(MOLRESPONSE_HEADERS
        density.h
        property.h
        TDDFT.h
        global_functions.h
        timer.h
        basic_operators.h
        Plot_VTK.h
        response_parameters.h
        )
# Create the MADchem library
add_mad_library(response MOLRESPONSE_SOURCES MOLRESPONSE_HEADERS "chem" "madness/chem/molresponse/")

set(targetname MADresponse)

add_dependencies(applications-madness MADresponse)
install(TARGETS MADresponse DESTINATION "${MADNESS_INSTALL_BINDIR}")
#



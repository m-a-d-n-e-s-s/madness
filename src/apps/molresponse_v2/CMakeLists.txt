# src/apps/molresponse_v2/CMakeLists.txt

if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
    message(STATUS "Using C++17 or higher standard")
else()
    message(FATAL_ERROR "C++17 or higher is required.")
endif()

set(MOLRESPONSE_SOURCES2
    ResponsePreliminariesCalculator.cpp
    ResponseManager.cpp
)

set(MOLRESPONSE_HEADERS2
    ResponsePreliminariesCalculator.hpp
    ResponsePreliminaries.hpp
    GroundStateData.hpp
    AccuracyProtocols.hpp
    MolecularProperties.hpp
    ResponseManager.hpp
)

# Create the MADresponse2 library
add_mad_library(response2 MOLRESPONSE_SOURCES2 MOLRESPONSE_HEADERS2 "chem" "madness/chem/molresponse2/")
set(targetname MADresponse2)
add_dependencies(applications-madness test_preliminaries)
add_dependencies(applications-madness molresponse2)
add_mad_executable(test_preliminaries "test_preliminaries.cpp" "MADresponse2;MADchem")
add_mad_executable(molresponse2 "molresponse2.cpp" "MADresponse2;MADchem")

# Add install targets for executables
install(TARGETS test_preliminaries molresponse2 DESTINATION "${MADNESS_INSTALL_BINDIR}")


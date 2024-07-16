# src/apps/molresponse/testing


set(MY_EXECUTABLES
        mad-dft
)


foreach (EXECUTABLE ${MY_EXECUTABLES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADall_response;MADchem;")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()


# For testing script add MADGtest library



# The list of unit test source files
set(RESPONSE_TEST_SOURCES test_coordinator test_manager
)

foreach (EXECUTABLE ${RESPONSE_TEST_SOURCES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADall_response;MADchem")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()

# Copy the test resources found in resources directory
file(COPY ${CMAKE_SRC_DIR}resources DESTINATION ${CMAKE_BINARY_DIR}/src/apps/molresponse/response_coordinator2)





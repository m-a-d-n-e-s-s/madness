# src/apps/moldft

set(MY_EXECUTABLES
	madqc
)

foreach (EXECUTABLE ${MY_EXECUTABLES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADall_response;MADchem;")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()


install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")
# For testing script



# The list of unit test source files
set(RESPONSE_TEST_SOURCES test_managers 
)

foreach (EXECUTABLE ${RESPONSE_TEST_SOURCES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADall_response;MADchem")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()




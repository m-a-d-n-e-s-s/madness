# src/apps/moldft

set(MY_EXECUTABLES
	madqc2
)

foreach (EXECUTABLE ${MY_EXECUTABLES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADresponse2;MADchem;")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()


install(TARGETS madqc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
# For testing script



# src/apps/molresponse/testing


set(MY_EXECUTABLES
        mad-freq
        mad-excited
        frequency_calc
        excited_state_calc
        full_frequency_test
        full_excited_test
        full_ground_test
        generate_excited_data
        generate_frequency_data
        test_development
        x_space_testing
        )

add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem")

foreach(EXECUTABLE ${MY_EXECUTABLES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADall_response;MADchem")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()

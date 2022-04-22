# src/apps/molresponse/testing

add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)


# Uses [mol] [xc] [operator] as input
add_mad_executable(mad-freq "mad-freq.cpp" "MADall_response;MADchem")
add_mad_executable(mad-excited "mad-excited.cpp" "MADall_response;MADchem")


# Uses "response.in" file as input
add_mad_executable(frequency_calc "frequency_calc.cpp" "MADall_response;MADchem")
add_mad_executable(excited_calc "excited_state_calc.cpp" "MADall_response;MADchem")

# Uses [xc] [operator] to create the full test suite
add_mad_executable(full_frequency_test "full_frequency_test.cpp" "MADall_response;MADchem")
add_mad_executable(full_excited_test "full_excited_test.cpp" "MADchem;MADresponse_base")
add_mad_executable(full_ground_test.cpp "full_ground_test.cpp" "MADchem;MADresponse_base")



add_mad_executable(generate_data "generate_frequency_data.cpp" "MADall_response;MADchem")
add_mad_executable(database-test "database_testing.cpp" "MADall_response;MADchem")
# Add dependencies for MADchem

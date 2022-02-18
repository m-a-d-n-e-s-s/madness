# src/apps/molresponse/testing

add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)




add_mad_executable(testing_libxc "testing_libxc_response.cpp" "MADall_response;MADchem")
add_mad_executable(generate_data "generate_frequency_answers.cpp" "MADall_response;MADchem")
add_mad_executable(single_molecule_response_test "single_molecule_test.cpp" "MADall_response;MADchem")
add_mad_executable(test_frequency_response "large_response_test.cpp" "MADall_response;MADchem")
add_mad_executable(test_response_running "test_response_running.cpp" "MADall_response;MADchem")
add_mad_executable(test_frequency_base "test_frequency_base.cpp" "MADall_response;MADchem")
add_mad_executable(test_excited "test_response_excited.cpp" "MADchem;MADresponse_base")
# Add dependencies for MADchem

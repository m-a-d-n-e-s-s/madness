# src/apps/molresponse/testing

add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)




add_mad_executable(test_response_running "test_response_running.cpp" "MADall_response;MADchem")
add_mad_executable(test_excited "test_response_excited.cpp" "MADchem;MADresponse_base")
# Add dependencies for MADchem

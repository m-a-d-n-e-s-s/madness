add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)




add_mad_executable(test_response_base "test_response_base.cpp" "MADchem;MADresponse_base")
add_mad_executable(test_excited_response "test_excited_response.cpp" "MADresponse_base;MADchem")
# Add dependencies for MADchem

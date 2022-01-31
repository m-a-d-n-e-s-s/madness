add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)
add_mad_executable(test_response_parameters "madness_catch_main.cc;test_response_parameters.cpp" "MADresponse;MADchem")




add_mad_executable(response_development "madness_catch_main_init.cc;test_response_base.cpp" "MADchem;MADresponse_base")
# Add dependencies for MADchem

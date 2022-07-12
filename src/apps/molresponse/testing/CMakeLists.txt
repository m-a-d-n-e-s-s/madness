add_mad_executable(test_schema_json "madness_catch_main.cc;qcschema_json_testing.cpp" "MADresponse;MADchem" ,)
add_mad_executable(test_response_parameters "madness_catch_main.cc;test_response_parameters.cpp" "MADresponse;MADchem")
add_mad_executable(test_response_functions "madness_catch_main.cc;test_response_functions.cpp" "MADresponse;MADchem")
add_mad_executable(test_density_vector "madness_catch_main.cc;test_density_vector.cpp;../ResponseBase.cpp;../densityVector.cpp;../global_functions.cc;../response_parameters.cpp;../timer.cc;../basic_operators.cc;../Plot_VTK.cc;../property.cc;" "MADchem")



# src/apps/molresponse


add_executable(molresponse
molresponse.cc
density.cc
property.cc
TDDFT.cc
global_functions.cc
calc_runner.cc
timer.cc
basic_operators.cc
Plot_VTK.cc
iter_freq2.cc
iterate_excited.cc
Lambda_X.cc
Theta_X.cc
iterate_gamma.cc
)

add_executable(test_schema_json
        madness_catch_main.cc
        qcschema_json_testing.cpp
        density.cc
        property.cc
        TDDFT.cc
        global_functions.cc
        calc_runner.cc
        timer.cc
        basic_operators.cc
        Plot_VTK.cc
        iter_freq2.cc
        iterate_excited.cc
        Lambda_X.cc
        Theta_X.cc
        iterate_gamma.cc
        )
target_link_libraries(molresponse -lstdc++fs MADchem  )
target_link_libraries(test_schema_json -lstdc++fs MADchem  )


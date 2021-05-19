# src/apps/molresponse

add_executable(test_resp
test_response.cc
density.cc
property.cc
TDDFT.cc
global_functions.cc
iterate_excited.cc
)
target_link_libraries(test_resp -lstdc++fs MADchem  )
add_executable(molresponse
molresponse.cc
)
#density.cc
# property.cc
# TDDFT.cc
# iterate_freq2.cc
# iterate_excited2.cc
# iterate_xy.cc
# iterate_gamma.cc
# Lambda_X.cc
# Theta_X.cc
#global_functions.cc
	#  basic_operators
	#  Plot_VTK.cc timer.cc)
target_link_libraries(molresponse -lstdc++fs MADchem  )


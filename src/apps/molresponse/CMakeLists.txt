# src/apps/molresponse

add_executable(molresponse molresponse.cc density.cc property.cc TDDFT.cc iterate_freq.cc iterate_xy.cc global_functions.cc basic_operators Plot_VTK.cc)
target_link_libraries(molresponse -lstdc++fs MADchem  )

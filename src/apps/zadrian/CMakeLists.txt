# src/apps/adrian

add_executable(densityTest densityTest.cc density.cc property.cc TDDFT.cc iterate_freq.cc iterate_xy.cc global_functions.cc TDHF_Basic_Operators2.cc Plot_VTK.cc)
target_link_libraries(densityTest -lstdc++fs MADchem  )

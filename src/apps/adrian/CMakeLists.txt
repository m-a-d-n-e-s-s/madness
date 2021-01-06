# src/apps/adrian

add_executable(densityTest densityTest.cc density.cc property.cc TDDFT.cc global_functions.cc TDHF_Basic_Operators2.cc Plot_VTK.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(densityTest -lstdc++fs MADchem  )

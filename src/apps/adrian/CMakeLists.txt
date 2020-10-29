# src/apps/adrian

add_executable(adrianXY adrianXY.cc TDDFT.cc TDHF_Basic_Operators2.cc Plot_VTK.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(adrianXY MADchem)

add_executable(densityTest densityTest.cc density.cc TDDFT.cc  property_operators.cc property_functions.cc global_functions.cc TDHF_Basic_Operators2.cc Plot_VTK.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(densityTest MADchem)

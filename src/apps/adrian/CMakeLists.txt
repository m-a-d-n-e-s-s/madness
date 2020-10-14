# src/apps/moldft

add_executable(adrianXY adrianXY.cc TDDFT.cc TDHF_Basic_Operators2.cc Plot_VTK.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(adrianXY MADchem)


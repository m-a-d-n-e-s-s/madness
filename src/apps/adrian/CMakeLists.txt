# src/apps/moldft

add_executable(adrianXY adrianXY.cc TDHF2.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(adrianXY MADchem)


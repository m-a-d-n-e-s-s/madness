# src/apps/moldft

add_executable(adrianXY adrianXY.cc TDDFT.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(adrianXY MADchem)


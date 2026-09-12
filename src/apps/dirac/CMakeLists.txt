# src/apps/dirac

add_mad_executable(DFdriver "DFdriver.cc;DF.cc;fcwf.cc" MADchem)

if(BUILD_TESTING)
  add_unittests(dirac "test_dirac.cc" "MADchem;MADgtest" "unittests;short")
endif()

#target_link_libraries(DFdriver MADchem MADmra)

add_dependencies(applications-madness DFdriver)

install(TARGETS DFdriver DESTINATION "${MADNESS_INSTALL_BINDIR}")


#add_executable(rk rk.cc)
#target_link_libraries(rk MADchem MADmra -lgmp -lmpfr)

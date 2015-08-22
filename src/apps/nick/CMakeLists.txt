# src/apps/nick

add_executable(projPsi projPsi.cc extra.cc wavef.cc hyp.cc mpreal.cc)
target_link_libraries(projPsi MADmra)

add_executable(lineplot lineplot.cc)
target_link_libraries(lineplot MADmra)

add_executable(test_nick test.cc wavef.cc hyp.cc mpreal.cc)
target_link_libraries(test_nick MADmra)

add_executable(wfSlice wfSlice.cc)
target_link_libraries(wfSlice MADmra)

add_executable(toDX toDX.cc wavef.cc hyp.cc mpreal.cc)
target_link_libraries(toDX MADmra)

add_executable(hello hello.cc wavef.cc hyp.cc)
target_link_libraries(hello MADmra)

#install(TARGETS projPsi test hello lineplot wfSlice toDX
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")
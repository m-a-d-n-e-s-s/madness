# src/apps/nick

add_mad_executable(projPsi projPsi.cc extra.cc wavef.cc hyp.cc mpreal.cc MADmra)

add_mad_executable(lineplot lineplot.cc MADmra)

add_mad_executable(test_nick test.cc wavef.cc hyp.cc mpreal.cc MADmra)

add_mad_executable(wfSlice wfSlice.cc MADmra)

add_mad_executable(toDX toDX.cc wavef.cc hyp.cc mpreal.cc MADmra)

add_mad_executable(hello hello.cc wavef.cc hyp.cc MADmra)

#install(TARGETS projPsi test hello lineplot wfSlice toDX
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")

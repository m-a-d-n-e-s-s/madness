# src/apps

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(plot)
add_subdirectory(tdse)
add_subdirectory(moldft)
# molresponse v1, the frozen legacy baseline, and the v2 engine were REMOVED
# (M1: v3 is standalone; parity + full-tensor Raman references archived in
# madness_studies/refs/parity_baseline_v2 before deletion).
add_subdirectory(molresponse)
add_subdirectory(madqc)
add_subdirectory(pno)
add_subdirectory(nemo)
add_subdirectory(cis)
add_subdirectory(dirac)
add_subdirectory(mp2)
add_subdirectory(cc2)
add_subdirectory(oep)
add_subdirectory(znemo)
add_subdirectory(zcis)
#add_subdirectory(hf)
#add_subdirectory(interior_bc)

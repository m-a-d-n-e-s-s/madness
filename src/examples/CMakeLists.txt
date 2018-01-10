# src/examples

set(EXAMPLE_SOURCES
    h2dft hedft hello hatom_energy h2 he tdse_example heat heat2 csqrt 
    sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio 
    dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic testsolver
    testspectralprop dielectric_external_field mp2 tiny oep h2dynamic newsolver
    cc2 nemo helium_exact density_smoothing siam_example lrccs ac_corr)
 
if(LIBXC_FOUND)
  list(APPEND EXAMPLE_SOURCES hefxc)
endif()

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_executable(${example} ${example}.cc)
  target_link_libraries(${example} MADchem)
endforeach()

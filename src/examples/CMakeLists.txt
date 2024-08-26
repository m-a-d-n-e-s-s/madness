# src/examples


add_subdirectory(writecoeffs)

set(EXAMPLE_SOURCES
    madinfo h2dft hedft hello hatom_energy h2 he tdse_example heat heat2 csqrt hatom hatom_sf_dirac
    sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio 
    dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic testsolver
    testspectralprop dielectric_external_field tiny h2dynamic newsolver testcomplexfunctionsolver
    helium_exact density_smoothing siam_example ac_corr dirac-hatom
    derivatives array_worldobject writecoeff writecoeff2)
 
if(LIBXC_FOUND)
  list(APPEND EXAMPLE_SOURCES hefxc)
endif()

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_mad_executable(${example} ${example}.cc "MADchem")
  add_dependencies(numerical-examples-madness ${example})
endforeach()

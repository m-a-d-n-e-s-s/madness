# src/examples


add_subdirectory(writecoeffs)

set(EXAMPLE_SOURCES
    madinfo h2dft hedft hello hatom_energy h2 he tdse_example heat heat2 csqrt hatom hatom_sf_dirac
    sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio 
    dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic
    testspectralprop dielectric_external_field tiny h2dynamic newsolver testcomplexfunctionsolver
    helium_exact density_smoothing siam_example ac_corr dirac-hatom
    derivatives array_worldobject pcr)
 
if(LIBXC_FOUND)
  list(APPEND EXAMPLE_SOURCES hefxc)
endif()

# MRA dimensions beyond 3 each example needs; see madness_mra_filter_targets()
foreach(_ex ac_corr csqrt hatom_1d nonlinschro sininteg tdse1d testcomplexfunctionsolver testspectralprop)
  set(mra_dims_${_ex} 1)
endforeach()
foreach(_ex he helium_exact tiny)
  set(mra_dims_${_ex} 6)
endforeach()
madness_mra_filter_targets(EXAMPLE_SOURCES)

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_mad_executable(${example} ${example}.cc "MADchem")
  add_dependencies(numerical-examples-madness ${example})
endforeach()

if(NOT (MADNESS_BUILD_MADWORLD_ONLY OR MADNESS_BUILD_LIBRARIES_ONLY) AND LAPACK_FOUND)
	add_subdirectory(periodic)
	# quantum-chemistry example decks (label `qctest`); needs the madqc target
	# from src/apps, which src/CMakeLists.txt processes before this directory
	add_subdirectory(qc)
endif()


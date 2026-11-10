# src/examples

set(EXAMPLE_SOURCES
    writecoeff writecoeff2 writecoeff3 writecoeff_json h2_write h2_write_json h2_write_json_2 writecoeff_hdf5)
# MRA dimensions beyond 3 each example needs; see madness_mra_filter_targets()
foreach(_ex writecoeff writecoeff2 writecoeff3 writecoeff_json writecoeff_hdf5)
  set(mra_dims_${_ex} 2)
endforeach()
madness_mra_filter_targets(EXAMPLE_SOURCES)

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_mad_executable(${example} ${example}.cc "MADchem")
  add_dependencies(numerical-examples-madness ${example})
endforeach()

# src/examples

set(EXAMPLE_SOURCES
    writecoeff writecoeff2 writecoeff3 h2_write)

# Create executables for example applications
foreach(example ${EXAMPLE_SOURCES})
  add_mad_executable(${example} ${example}.cc "MADchem")
  add_dependencies(numerical-examples-madness ${example})
endforeach()


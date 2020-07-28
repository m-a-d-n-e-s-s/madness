# src/madness/external

if(ENABLE_UNITTESTS)
  add_subdirectory(gtest)
endif()

if (NOT MADNESS_BUILD_MADWORLD_ONLY)
  add_subdirectory(muParser)
  add_subdirectory(tinyxml)
endif (NOT MADNESS_BUILD_MADWORLD_ONLY)

# if want embedded elemental ...
if(ENABLE_ELEMENTAL AND ENABLE_ELEMENTAL_EMBEDDED)

  # Note: Here we are changing the build configuration of for Elemental since it
  # requires [Hybrid|Pure][Debug|Release] as the build configuratino. This
  # should not be required if we upgrade Elemental to v0.85.
  
  # Store the build configuratin
  set(PUSH_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  
  # Set the build type used by Elemental. We use Pure (MPI only/no OpenMP) since
  # MADNESS threads do not play nice with Elemental threads.
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "PureDebug" 
        CACHE STRING "Specifies the build type on single-configuration generators" FORCE)
  else()
    set(CMAKE_BUILD_TYPE "PureRelease" 
        CACHE STRING "Specifies the build type on single-configuration generators" FORCE)
  endif()
  
  # Set the configuration variables used by elemental
  set(MATH_LIBS "${LAPACK_LIBRARIES}")
  if(((NOT ENABLE_SPINLOCKS) OR ENABLE_NEVER_SPIN) OR (CMAKE_SYSTEM_NAME MATCHES "Darwin"))
    set(HAVE_SPINLOCKS OFF CACHE BOOL "Enable if pthread lib supports spinlocks")
  else()
    set(HAVE_SPINLOCKS ON CACHE BOOL "Enable if pthread lib supports spinlocks")
  endif()

  message(STATUS "Begin Elemental Configure")
  add_subdirectory(elemental)
  message(STATUS "End Elemental Configure")
  set(ELEMENTAL_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/elemental/include
      ${CMAKE_CURRENT_BINARY_DIR}/elemental/include PARENT_SCOPE)
  set(ELEMENTAL_LIBRARIES elemental PARENT_SCOPE)
  
  # Restore the build configuration.
  set(CMAKE_BUILD_TYPE ${PUSH_CMAKE_BUILD_TYPE}
      CACHE STRING "Specifies the build type on single-configuration generators" FORCE) 

  set(ELEMENTAL_FOUND 1 PARENT_SCOPE)
  set(ELEMENTAL_PACKAGE_NAME elemental CACHE INTERNAL "Elemental package name (component name is same)" FORCE)
  set(HAVE_ELEMENTAL_H ON CACHE INTERNAL "Have elemental.hpp" FORCE)
  
  # Add build-elemental target
  add_custom_target(build-elemental ALL
      COMMAND ${CMAKE_COMMAND} --build . --target elemental
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      COMMENT "Building elemental"
      USES_TERMINAL)
  
  # Add a component install target for elemental 
  add_custom_target(install-elemental
      COMMAND ${CMAKE_COMMAND} -DCOMPONENT=elemental -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing elemental library components"
      USES_TERMINAL)
  add_dependencies(install-elemental elemental)
  
  # elemental target has been added to the list of exported targets in hacked CMakeFiles.txt in elemental subdir
  
endif(ENABLE_ELEMENTAL AND ENABLE_ELEMENTAL_EMBEDDED)

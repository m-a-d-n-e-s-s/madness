# src/madness

### MADWorld runtime
add_subdirectory(world)
# may need to build elemental as part of MADWorld
add_subdirectory(external)
if (ENABLE_ELEMENTAL)
	if (MADNESS_HAS_ELEMENTAL_EMBEDDED)
		include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/elemental/include
				${CMAKE_CURRENT_BINARY_DIR}/external/elemental/include)
	elseif (MADNESS_HAS_ELEMENTAL)
		include_directories(${ELEMENTAL_SOURCE_DIR}/include
				${ELEMENTAL_SOURCE_DIR}/external/suite_sparse/include
				${ELEMENTAL_BINARY_DIR}/include)
	endif ()
endif(ENABLE_ELEMENTAL)
if(MADNESS_BUILD_MADWORLD_ONLY)  # madness library includes just madworld
	add_library(madness $<TARGET_OBJECTS:MADworld-obj>)
endif(MADNESS_BUILD_MADWORLD_ONLY)

### CLAPACK will be build if BLAS/LAPACK is available
if (LAPACK_FOUND)
  set(MADCLAPACK_HEADERS tensor/cblas.h tensor/cblas_types.h
		tensor/clapack.h tensor/clapack_fortran.h
		tensor/lapacke_types.h) # this part of MADlinalg is purely independent of MADtensor
  # add custom target for "installing" MADclapack headers
  install(FILES ${MADCLAPACK_HEADERS}
		DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness/tensor"
		COMPONENT clapack)
  add_custom_target(install-madness-clapack
		COMMAND ${CMAKE_COMMAND} -DCOMPONENT=clapack -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
		COMMENT "Installing clapack headers"
		USES_TERMINAL)
endif(LAPACK_FOUND)

### numerics + apps
if (NOT MADNESS_BUILD_MADWORLD_ONLY)
	add_subdirectory(misc)
	add_subdirectory(tensor)
	add_subdirectory(mra)
	add_subdirectory(chem)
    # Create the super library that contains all MADNESS libraries
	add_library(madness
    	$<TARGET_OBJECTS:MADworld-obj> $<TARGET_OBJECTS:MADmisc-obj>
    	$<TARGET_OBJECTS:MADtensor-obj> ${MADTENSOR_MTXMQ_OBJ}
    	$<TARGET_OBJECTS:MADlinalg-obj> $<TARGET_OBJECTS:MADmuparser-obj>
    	$<TARGET_OBJECTS:MADtinyxml-obj> $<TARGET_OBJECTS:MADmra-obj>
		$<TARGET_OBJECTS:MADchem-obj>)
endif(NOT MADNESS_BUILD_MADWORLD_ONLY)

# Set the dependencies for the super-library
target_compile_features(madness INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")
# DISABLEPIE flag can break linking of dependent libraries (e.g. on Linux using gcc6)
# instead for each dependent executable target T do:
# target_link_libraries(T MADNESS_DISABLEPIE_LINKER_FLAGS)
#
#if(LINKER_HAS_DISABLEPIE_SUPPORT)
#  target_link_libraries(madness INTERFACE "${DISABLEPIE_LINKER_FLAG}")
#endif()
if(GPERFTOOLS_FOUND)
  target_include_directories(madness PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${GPERFTOOLS_LIBRARIES})
endif()
if(LIBUNWIND_FOUND AND NOT GPERFTOOLS_LIBRARIES MATCHES "unwind")
  target_include_directories(madness PUBLIC ${LIBUNWIND_INCLUDE_DIR})
  target_link_libraries(madness PUBLIC ${LIBUNWIND_LIBRARIES})
endif()
if(ELEMENTAL_FOUND)
  target_link_libraries(madness PUBLIC ${ELEMENTAL_PACKAGE_NAME})
endif()
if(PAPI_FOUND)
  target_include_directories(madness PUBLIC ${PAPI_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${PAPI_LIBRARIES})
endif()
if(TBB_FOUND AND TARGET TBB::tbb)
  target_link_libraries(madness PUBLIC TBB::tbb)
  if (TBB_COMPILE_FLAGS)
    target_compile_definitions(madness PUBLIC "${TBB_COMPILE_FLAGS}")
  endif(TBB_COMPILE_FLAGS)
endif()
if(TARGET PaRSEC::parsec)
  target_link_libraries(madness PUBLIC PaRSEC::parsec)
endif()
if(MPI_FOUND)
  target_include_directories(madness PUBLIC ${MPI_INCLUDE_PATH})
  target_compile_definitions(madness PUBLIC -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1)
  target_compile_options(madness PUBLIC "${MPI_COMPILE_FLAGS}")
  target_link_libraries(madness PUBLIC ${MPI_LINK_FLAGS})
  target_link_libraries(madness PUBLIC ${MPI_LIBRARIES})
endif()
if(HAVE_IBMBGQ OR HAVE_IBMBGP)
  target_link_libraries(madness PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests/libMADMTXM.a)
endif()
if (LAPACK_FOUND)

	# if LAPACK_LIBRARIES includes IMPORTED targets make sure they are available
	if (LAPACK_LIBRARIES MATCHES OpenMP::OpenMP_C AND NOT TARGET OpenMP::OpenMP_C)
		find_package(OpenMP REQUIRED COMPONENTS C)
	endif()

	target_link_libraries(madness PUBLIC ${LAPACK_LIBRARIES})
	if (LAPACK_INCLUDE_DIRS)
		target_include_directories(madness PUBLIC ${LAPACK_INCLUDE_DIRS})
	endif(LAPACK_INCLUDE_DIRS)
	if (LAPACK_COMPILE_OPTIONS)
		target_compile_options(madness PUBLIC ${LAPACK_COMPILE_OPTIONS})
	endif(LAPACK_COMPILE_OPTIONS)
	if (LAPACK_COMPILE_DEFINITIONS)
		target_compile_definitions(madness PUBLIC ${LAPACK_COMPILE_DEFINITIONS})
	endif(LAPACK_COMPILE_DEFINITIONS)
endif(LAPACK_FOUND)
target_link_libraries(madness PUBLIC Threads::Threads)
if (WORLD_GET_DEFAULT_DISABLED)
  target_compile_definitions(madness PUBLIC -DMADNESS_DISABLE_WORLD_GET_DEFAULT=1)
endif (WORLD_GET_DEFAULT_DISABLED)
# MADchem dependencies
if (LIBXC_FOUND)
	target_include_directories(madness PUBLIC ${LIBXC_INCLUDE_DIRS})
	target_link_libraries(madness PUBLIC ${LIBXC_LIBRARIES})
endif ()
if (PCM_FOUND)
	target_include_directories(madness PUBLIC ${PCM_INCLUDE_DIRS})
	target_link_libraries(madness PUBLIC ${PCM_LIBRARIES})
endif ()
if (BOOST_FOUND)
	target_include_directories(madness PUBLIC ${BOOST_INCLUDE_DIRS})
	target_link_libraries(madness PUBLIC ${BOOST_LIBRARIES})
endif ()

# Add common header files and madness library to the install list
install(FILES constants.h fortran_ctypes.h madness_config.h ${CMAKE_CURRENT_BINARY_DIR}/config.h
    DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness"
    COMPONENT common)
add_custom_target(install-madness-common
     COMMAND ${CMAKE_COMMAND} -DCOMPONENT=common -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
     COMMENT "Installing common library components"
     USES_TERMINAL)
     
install(TARGETS madness EXPORT madness
        COMPONENT madness
        LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")
add_custom_target(install-madness
        COMMAND ${CMAKE_COMMAND} -DCOMPONENT=madness -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
        COMMENT "Installing madness library components"
        USES_TERMINAL)
#add_dependencies(install-madness madness)
add_dependencies(install-madness install-madness-common)
add_dependencies(install-madness-libraries install-madness)

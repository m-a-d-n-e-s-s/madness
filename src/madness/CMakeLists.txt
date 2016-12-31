# src/madness

# Add source subdirectories
add_subdirectory(external)

if (MADNESS_HAS_ELEMENTAL_EMBEDDED)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/elemental/include
      ${CMAKE_CURRENT_BINARY_DIR}/external/elemental/include)
elseif (MADNESS_HAS_ELEMENTAL)
  include_directories(${ELEMENTAL_SOURCE_DIR}/include
      ${ELEMENTAL_SOURCE_DIR}/external/suite_sparse/include
      ${ELEMENTAL_BINARY_DIR}/include)
endif ()

add_subdirectory(world)
add_subdirectory(misc)
add_subdirectory(tensor)
add_subdirectory(mra)

# Create the super library that contains all MADNESS libraries
add_library(madness
    $<TARGET_OBJECTS:MADworld-obj> $<TARGET_OBJECTS:MADmisc-obj>
    $<TARGET_OBJECTS:MADtensor-obj> ${MADTENSOR_MTXMQ_OBJ}
    $<TARGET_OBJECTS:MADlinalg-obj> $<TARGET_OBJECTS:MADmuparser-obj>
    $<TARGET_OBJECTS:MADtinyxml-obj> $<TARGET_OBJECTS:MADmra-obj>)

# Set the dependencies for the super-library
if(CXX11_COMPILE_FLAG)
  target_compile_options(madness INTERFACE $<INSTALL_INTERFACE:${CXX11_COMPILE_FLAG}>)
endif()
if(LINKER_HAS_DISABLEPIE_SUPPORT)
  target_link_libraries(madness INTERFACE "${DISABLEPIE_LINKER_FLAG}")
endif()
if(GPERFTOOLS_FOUND)
  target_include_directories(madness PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${GPERFTOOLS_LIBRARIES})
endif()
if(LIBUNWIND_FOUND AND NOT GPERFTOOLS_LIBRARIES MATCHES "unwind")
  target_include_directories(madness PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${LIBUNWIND_LIBRARIES})
endif()
if(ELEMENTAL_FOUND)
  target_link_libraries(madness PUBLIC ${ELEMENTAL_PACKAGE_NAME})
endif()
if(PAPI_FOUND)
  target_include_directories(madness PUBLIC ${PAPI_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${PAPI_LIBRARIES})
endif()
if(TBB_FOUND)
  target_include_directories(madness PUBLIC ${TBB_INCLUDE_DIRS})
  if(TBB_USE_DEBUG AND TBB_LIBRARIES_DEBUG)
    target_compile_definitions(madness PUBLIC "${TBB_COMPILE_FLAGS_DEBUG}")
    target_link_libraries(madness PUBLIC ${TBB_LIBRARIES_DEBUG})
  else()
    target_link_libraries(madness PUBLIC ${TBB_LIBRARIES_RELEASE})
  endif()
endif()
if(PARSEC_FOUND)
  target_include_directories(madness PUBLIC ${PARSEC_INCLUDE_DIRS})
  target_link_libraries(madness PUBLIC ${PARSEC_LIBRARIES} ${PARSEC_LIBRARIES_DEP}  ${PARSEC_LIBRARY_DIRS_DEP_FLAGS})
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
target_link_libraries(madness PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(madness PUBLIC ${CMAKE_THREAD_LIBS_INIT})
if (WORLD_GET_DEFAULT_DISABLED)
  target_compile_definitions(madness PUBLIC -DMADNESS_DISABLE_WORLD_GET_DEFAULT=1)
endif (WORLD_GET_DEFAULT_DISABLED)

# Add common header files and madness library to the install list
install(TARGETS madness EXPORT madness
    COMPONENT madness
    LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")
install(FILES constants.h fortran_ctypes.h madness_config.h ${CMAKE_CURRENT_BINARY_DIR}/config.h
    DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness"
    COMPONENT common)
add_custom_target(install-common
     COMMAND ${CMAKE_COMMAND} -DCOMPONENT=common -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
     COMMENT "Installing common library components")


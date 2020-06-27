# src/madness/world

set(MADWORLD_HEADERS 
    info.h archive.h print.h worldam.h future.h worldmpi.h
    world_task_queue.h array_addons.h stack.h vector.h worldgop.h 
    world_object.h buffer_archive.h nodefaults.h dependency_interface.h 
    worldhash.h worldref.h worldtypes.h dqueue.h parallel_archive.h 
    vector_archive.h madness_exception.h worldmem.h thread.h worldrmi.h 
    safempi.h worldpapi.h worldmutex.h print_seq.h worldhashmap.h range.h 
    atomicint.h posixmem.h worldptr.h deferred_cleanup.h MADworld.h world.h 
    uniqueid.h worldprofile.h timers.h binary_fstream_archive.h mpi_archive.h 
    text_fstream_archive.h worlddc.h mem_func_wrapper.h taskfn.h group.h 
    dist_cache.h distributed_id.h type_traits.h function_traits.h stubmpi.h 
    bgq_atomics.h binsorter.h parsec.h meta.h worldinit.h)
set(MADWORLD_SOURCES
    madness_exception.cc world.cc timers.cc future.cc redirectio.cc
    archive_type_names.cc info.cc debug.cc print.cc worldmem.cc worldrmi.cc
    safempi.cc worldpapi.cc worldref.cc worldam.cc worldprofile.cc thread.cc 
    world_task_queue.cc worldgop.cc deferred_cleanup.cc worldmutex.cc
    binary_fstream_archive.cc text_fstream_archive.cc lookup3.c worldmpi.cc 
    group.cc parsec.cc archive.cc)

# Create the MADworld-obj and MADworld library targets
add_mad_library(world MADWORLD_SOURCES MADWORLD_HEADERS "common;${ELEMENTAL_PACKAGE_NAME}" "madness/world")

# DISABLEPIE flag can break linking of dependent libraries (e.g. on Linux using gcc6)
# instead for each dependent executable target T do:
# target_link_libraries(T ${MADNESS_DISABLEPIE_LINKER_FLAGS})
#
#if(LINKER_HAS_DISABLEPIE_SUPPORT)
#  target_link_libraries(MADworld INTERFACE "${DISABLEPIE_LINKER_FLAG}")
#endif()

  set(targetname MADworld)

  # create an additional, MADworld-only, export
  install(TARGETS ${targetname} EXPORT madworld
      COMPONENT world
      PUBLIC_HEADER DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness/world"
      LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")

  # Set the dependencies for the MADworld library
  #target_compile_definitions(${targetname} PUBLIC "-DMADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\"")
  target_include_directories(${targetname} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src/>)
  if (GPERFTOOLS_FOUND)
    target_include_directories(${targetname} PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${GPERFTOOLS_LIBRARIES})
  endif ()
  if (LIBUNWIND_FOUND AND NOT GPERFTOOLS_LIBRARIES MATCHES "unwind")
    target_include_directories(${targetname} PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${LIBUNWIND_LIBRARIES})
  endif ()
  if (ELEMENTAL_FOUND)
    target_link_libraries(${targetname} PUBLIC ${ELEMENTAL_PACKAGE_NAME})
  endif ()
  if (PAPI_FOUND)
    target_include_directories(${targetname} PUBLIC ${PAPI_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${PAPI_LIBRARIES})
  endif ()
  if (TBB_FOUND)
    target_include_directories(${targetname} PUBLIC ${TBB_INCLUDE_DIRS})
    if (TBB_USE_DEBUG AND TBB_LIBRARIES_DEBUG)
      target_compile_definitions(${targetname} PUBLIC "${TBB_COMPILE_FLAGS_DEBUG}")
      target_link_libraries(${targetname} PUBLIC ${TBB_LIBRARIES_DEBUG})
    else ()
      target_link_libraries(${targetname} PUBLIC ${TBB_LIBRARIES_RELEASE})
    endif ()
    if (TBB_COMPILE_FLAGS)
      target_compile_definitions(${targetname} PUBLIC "${TBB_COMPILE_FLAGS}")
    endif(TBB_COMPILE_FLAGS)
  endif ()
  if (TARGET PaRSEC::parsec)
    target_link_libraries(${targetname} PUBLIC PaRSEC::parsec)
  endif ()
  if (MPI_FOUND)
    target_include_directories(${targetname} PUBLIC ${MPI_INCLUDE_PATH})
    target_compile_definitions(${targetname} PUBLIC -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1 -DMADNESS_MPI_HEADER=$<1:"${MADNESS_MPI_HEADER}">)
    target_compile_options(${targetname} PUBLIC ${MPI_COMPILE_FLAGS})
    target_link_libraries(${targetname} PUBLIC ${MPI_LINK_FLAGS})
    target_link_libraries(${targetname} PUBLIC ${MPI_LIBRARIES})
  endif ()
  target_link_libraries(${targetname} PUBLIC Threads::Threads)
  if (WORLD_GET_DEFAULT_DISABLED)
    target_compile_definitions(${targetname} PUBLIC -DMADNESS_DISABLE_WORLD_GET_DEFAULT=1)
  endif (WORLD_GET_DEFAULT_DISABLED)

if (TARGET PaRSEC::parsec)
  add_subproject_parsec_target_to_export_set("madworld" "world")
endif()

if(ENABLE_UNITTESTS)

  # The list of unit test source files
  set(WORLD_TEST_SOURCES test_prof.cc test_ar.cc test_hashdc.cc test_hello.cc
      test_atomicint.cc test_future.cc test_future2.cc test_future3.cc 
      test_dc.cc test_hashthreaded.cc test_queue.cc test_world.cc 
      test_worldprofile.cc test_binsorter.cc test_vector.cc test_worldptr.cc 
      test_worldref.cc test_stack.cc test_googletest.cc test_tree.cc)

  add_unittests(world "${WORLD_TEST_SOURCES}" "MADworld;MADgtest")    

  find_package(CUDA)
  if (TARGET PaRSEC::parsec AND CUDA_FOUND)
    CMAKE_PUSH_CHECK_STATE()
    STRING(REPLACE "-std=c++1z" "" _tmp "${CMAKE_CXX_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${_tmp}")
    cuda_add_library(MADtest_cuda hello_world.cu EXCLUDE_FROM_ALL)
    CMAKE_POP_CHECK_STATE()
    target_link_libraries(test_world PRIVATE MADtest_cuda)
  endif ()

  set_tests_properties(world-test_googletest PROPERTIES WILL_FAIL TRUE)
  
endif()


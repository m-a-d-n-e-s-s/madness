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
    bgq_atomics.h binsorter.h parsec.h)
set(MADWORLD_SOURCES
    madness_exception.cc world.cc timers.cc future.cc redirectio.cc
    archive_type_names.cc info.cc debug.cc print.cc worldmem.cc worldrmi.cc
    safempi.cc worldpapi.cc worldref.cc worldam.cc worldprofile.cc thread.cc 
    world_task_queue.cc worldgop.cc deferred_cleanup.cc worldmutex.cc
    binary_fstream_archive.cc text_fstream_archive.cc lookup3.c worldmpi.cc 
    group.cc parsec.cc)

# Create the MADworld-obj and MADworld library targets
add_mad_library(world MADWORLD_SOURCES MADWORLD_HEADERS "common;${ELEMENTAL_PACKAGE_NAME}" "madness/world")

if(LINKER_HAS_DISABLEPIE_SUPPORT)
  target_link_libraries(MADworld INTERFACE "${DISABLEPIE_LINKER_FLAG}")
endif()

# create an additional, MADworld-only, export
install(TARGETS MADworld EXPORT madworld
      COMPONENT world
      PUBLIC_HEADER DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness/world"
      LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
      INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")

# Set the dependencies for the MADworld library
#target_compile_definitions(MADworld PUBLIC "-DMADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\"")
target_include_directories(MADworld PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src/>)
if(GPERFTOOLS_FOUND)
  target_include_directories(MADworld PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(MADworld PUBLIC ${GPERFTOOLS_LIBRARIES})
endif()
if(LIBUNWIND_FOUND AND NOT GPERFTOOLS_LIBRARIES MATCHES "unwind")
  target_include_directories(MADworld PUBLIC ${LIBUNWIND_INCLUDE_DIRS})
  target_link_libraries(MADworld PUBLIC ${LIBUNWIND_LIBRARIES})
endif()
if(ELEMENTAL_FOUND)
  target_link_libraries(MADworld PUBLIC ${ELEMENTAL_PACKAGE_NAME})
endif()
if(PAPI_FOUND)
  target_include_directories(MADworld PUBLIC ${PAPI_INCLUDE_DIRS})
  target_link_libraries(MADworld PUBLIC ${PAPI_LIBRARIES})
endif()
if(TBB_FOUND)
  target_include_directories(MADworld PUBLIC ${TBB_INCLUDE_DIRS})
  if(TBB_USE_DEBUG AND TBB_LIBRARIES_DEBUG)
    target_compile_definitions(MADworld PUBLIC "${TBB_COMPILE_FLAGS_DEBUG}")
    target_link_libraries(MADworld PUBLIC ${TBB_LIBRARIES_DEBUG})
  else()
    target_link_libraries(MADworld PUBLIC ${TBB_LIBRARIES_RELEASE})
  endif()
endif()
if(PARSEC_FOUND)
  target_include_directories(MADworld PUBLIC ${PARSEC_INCLUDE_DIRS})
  target_link_libraries(MADworld PUBLIC ${PARSEC_LIBRARIES} ${PARSEC_LIBRARIES_DEP})
  # convert library directories list to linker flags
  set(_parsec_libdir_flags "")
  foreach(_dir ${PARSEC_LIBRARY_DIRS_DEP})
    list(APPEND _parsec_libdir_flags "-L${_dir}")
  endforeach()
  set(PARSEC_LIBRARY_DIRS_DEP_FLAGS "${_parsec_libdir_flags}" CACHE INTERNAL "PARSEC(+deps) library directories as linker flag list")
  target_link_libraries(MADworld PUBLIC ${PARSEC_LIBRARY_DIRS_DEP_FLAGS})
endif()
if(MPI_FOUND)
  target_include_directories(MADworld PUBLIC ${MPI_INCLUDE_PATH})
  target_compile_definitions(MADworld PUBLIC -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1)
  target_compile_options(MADworld PUBLIC ${MPI_COMPILE_FLAGS})
  target_link_libraries(MADworld PUBLIC ${MPI_LINK_FLAGS})
  target_link_libraries(MADworld PUBLIC ${MPI_LIBRARIES})
endif()
target_link_libraries(MADworld PUBLIC ${CMAKE_THREAD_LIBS_INIT})
if (WORLD_GET_DEFAULT_DISABLED)
  target_compile_definitions(MADworld PUBLIC -DMADNESS_DISABLE_WORLD_GET_DEFAULT=1)
endif (WORLD_GET_DEFAULT_DISABLED)

if(ENABLE_UNITTESTS)

  # The list of unit test source files
  set(WORLD_TEST_SOURCES test_prof.cc test_ar.cc test_hashdc.cc test_hello.cc
      test_atomicint.cc test_future.cc test_future2.cc test_future3.cc 
      test_dc.cc test_hashthreaded.cc test_queue.cc test_world.cc 
      test_worldprofile.cc test_binsorter.cc test_vector.cc test_worldptr.cc 
      test_worldref.cc test_stack.cc test_googletest.cc test_tree.cc)


  add_unittests(world WORLD_TEST_SOURCES "MADworld;MADgtest")

  set_tests_properties(world-test_googletest PROPERTIES WILL_FAIL TRUE)

  if (ENABLE_PARSEC)
    find_package(CUDA)
    if (CUDA_FOUND) # no way to make sure PARSEC has CUDA
                    # so just look for it and hope for the best
      cuda_add_library(MADtest_cuda hello_world.cu EXCLUDE_FROM_ALL)  
      target_link_libraries(test_world MADtest_cuda)
    endif(CUDA_FOUND)
  endif (ENABLE_PARSEC)
  
endif()


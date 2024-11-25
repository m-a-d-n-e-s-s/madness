# src/madness/world

set(MADWORLD_HEADERS 
    archive.h print.h worldam.h future.h worldmpi.h
    world_task_queue.h array_addons.h stack.h vector.h worldgop.h 
    world_object.h buffer_archive.h nodefaults.h dependency_interface.h 
    worldhash.h worldref.h worldtypes.h dqueue.h parallel_archive.h parallel_dc_archive.h
    vector_archive.h madness_exception.h worldmem.h thread.h worldrmi.h 
    safempi.h worldpapi.h worldmutex.h print_seq.h worldhashmap.h range.h 
    atomicint.h posixmem.h worldptr.h deferred_cleanup.h MADworld.h world.h 
    uniqueid.h worldprofile.h timers.h binary_fstream_archive.h mpi_archive.h 
    text_fstream_archive.h worlddc.h mem_func_wrapper.h taskfn.h group.h 
    dist_cache.h distributed_id.h type_traits.h function_traits.h stubmpi.h 
    bgq_atomics.h binsorter.h parsec.h meta.h worldinit.h thread_info.h
    cloud.h test_utilities.h timing_utilities.h units.h)
set(MADWORLD_SOURCES
    madness_exception.cc world.cc timers.cc future.cc redirectio.cc
    archive_type_names.cc debug.cc print.cc worldmem.cc worldrmi.cc
    safempi.cc worldpapi.cc worldref.cc worldam.cc worldprofile.cc thread.cc 
    world_task_queue.cc worldgop.cc deferred_cleanup.cc worldmutex.cc
    binary_fstream_archive.cc text_fstream_archive.cc lookup3.c worldmpi.cc 
    group.cc parsec.cc archive.cc units.cc)

if(MADNESS_ENABLE_CEREAL)
    set(MADWORLD_HEADERS ${MADWORLD_HEADERS} "cereal_archive.h")
endif()

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

# N.B. If cereal is fetched it will be part of cereal export set
if(MADNESS_ENABLE_CEREAL)
  include(${PROJECT_SOURCE_DIR}/cmake/modules/FindOrFetchCereal.cmake)
  target_link_libraries(${targetname} PUBLIC cereal)
  target_compile_definitions(${targetname} PUBLIC MADNESS_HAS_CEREAL)
endif()

# Set the dependencies for the MADworld library
target_include_directories(${targetname} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src/>)
if (GPERFTOOLS_FOUND)
  target_include_directories(${targetname} PUBLIC ${GPERFTOOLS_INCLUDE_DIRS})
  target_link_libraries(${targetname} PUBLIC ${GPERFTOOLS_LIBRARIES})
endif ()
if (LIBUNWIND_FOUND AND NOT GPERFTOOLS_LIBRARIES MATCHES "unwind")
  target_include_directories(${targetname} PUBLIC ${LIBUNWIND_INCLUDE_DIR})
  target_link_libraries(${targetname} PUBLIC ${LIBUNWIND_LIBRARIES})
endif ()
if (ELEMENTAL_FOUND)
  target_link_libraries(${targetname} PUBLIC ${ELEMENTAL_PACKAGE_NAME})
endif ()
if (PAPI_FOUND)
  target_include_directories(${targetname} PUBLIC ${PAPI_INCLUDE_DIRS})
  target_link_libraries(${targetname} PUBLIC ${PAPI_LIBRARIES})
endif ()
if (TBB_FOUND AND TARGET TBB::tbb)  # make sure TBB was discovered by MADNESS
  target_link_libraries(${targetname} PUBLIC TBB::tbb)
  if (TBB_COMPILE_FLAGS)
    target_compile_definitions(${targetname} PUBLIC "${TBB_COMPILE_FLAGS}")
  endif(TBB_COMPILE_FLAGS)
endif ()
if (TARGET PaRSEC::parsec)
  target_link_libraries(${targetname} PUBLIC PaRSEC::parsec)
endif ()
if (MPI_FOUND)
  target_compile_definitions(${targetname} PUBLIC -DMADNESS_MPI_HEADER=$<1:"${MADNESS_MPI_HEADER}">)
  target_link_libraries(${targetname} PUBLIC MPI::MPI_CXX)
endif ()
target_link_libraries(${targetname} PUBLIC Threads::Threads)
if (WORLD_GET_DEFAULT_DISABLED)
  target_compile_definitions(${targetname} PUBLIC -DMADNESS_DISABLE_WORLD_GET_DEFAULT=1)
endif (WORLD_GET_DEFAULT_DISABLED)

# record the list of dependencies to be used in madness-config.cmake
get_target_property(MADNESS_MADWORLD_INTERFACE_LINK_LIBRARIES MADworld INTERFACE_LINK_LIBRARIES)
set(MADNESS_MADWORLD_INTERFACE_LINK_LIBRARIES "${MADNESS_MADWORLD_INTERFACE_LINK_LIBRARIES}" CACHE INTERNAL "public library dependencies of MADworld" FORCE)

if(BUILD_TESTING)

  # The list of unit test source files
  set(WORLD_TEST_SOURCES test_prof.cc test_ar.cc test_hashdc.cc test_hello.cc
      test_atomicint.cc test_future.cc test_future2.cc test_future3.cc 
      test_dc.cc test_hashthreaded.cc test_queue.cc test_world.cc 
      test_worldprofile.cc test_binsorter.cc test_vector.cc test_worldptr.cc 
      test_worldref.cc test_stack.cc test_googletest.cc test_tree.cc
          )

  add_unittests(world "${WORLD_TEST_SOURCES}" "MADworld;MADgtest" "unittests;short")

  if (TARGET PaRSEC::parsec AND PARSEC_HAVE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
      # cmake 3.17 decouples C++ and CUDA standards, see https://gitlab.kitware.com/cmake/cmake/issues/19123
      # cmake 3.18 knows that CUDA 11 provides cuda_std_17
      cmake_minimum_required(VERSION 3.18.0)
      set(CMAKE_CUDA_STANDARD 17)
      set(CMAKE_CUDA_EXTENSIONS OFF)
      set(CMAKE_CUDA_STANDARD_REQUIRED ON)
      set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
      # N.B. need relaxed constexpr for std::complex
      # see https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#constexpr-functions%5B/url%5D:
      if (DEFINED CMAKE_CUDA_FLAGS)
        set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr ${CMAKE_CUDA_FLAGS}")
      else()
        set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr")
      endif()

      enable_language(CUDA)
      add_library(MADtest_cuda hello_world.cu)
      target_link_libraries(test_world PRIVATE MADtest_cuda)
    endif(CMAKE_CUDA_COMPILER)
  endif ()

  set_tests_properties(madness/test/world/test_googletest/run PROPERTIES WILL_FAIL TRUE)
  
endif()


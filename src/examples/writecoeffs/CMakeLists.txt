cmake_minimum_required(VERSION 3.0.0)
project(h5cpp-examples LANGUAGES CXX C)

get_directory_property(IS_TOP_LEVEL PARENT_DIRECTORY)
if (NOT IS_TOP_LEVEL)
  find_package(h5cpp REQUIRED)
  if (H5CPP_MPI_ENABLED)
    find_package(MPI REQUIRED)
    include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  endif ()
endif()

if (NOT H5CPP_BOOST_ENABLED)
  set(CMAKE_CXX_STANDARD 17)
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

add_executable(basic_files basic_files.cpp)
target_link_libraries(basic_files PRIVATE h5cpp::h5cpp)

add_executable(selecting_datasets selecting_datasets.cpp)
target_link_libraries(selecting_datasets PRIVATE h5cpp::h5cpp)

add_executable(std_vector_io std_vector_io.cpp)
target_link_libraries(std_vector_io PRIVATE h5cpp::h5cpp)

add_executable(complex_io complex_io.cpp)
target_link_libraries(complex_io PRIVATE h5cpp::h5cpp)

add_executable(writing_image writing_image.cpp
    rgbpixel.cpp)
target_link_libraries(writing_image PRIVATE h5cpp::h5cpp)

add_executable(write_single_vector write_single_vector.cpp)
target_link_libraries(write_single_vector PRIVATE h5cpp::h5cpp)

add_executable(read_dataset read_dataset.cpp)
target_link_libraries(read_dataset PRIVATE h5cpp::h5cpp)

add_executable(append_scalar_data append_scalar_data.cpp)
target_link_libraries(append_scalar_data PRIVATE h5cpp::h5cpp)

add_executable(write_vector_list write_vector_list.cpp)
target_link_libraries(write_vector_list PRIVATE h5cpp::h5cpp)

add_executable(append_vector_data append_vector_data.cpp)
target_link_libraries(append_vector_data PRIVATE h5cpp::h5cpp)


add_custom_target(examples)
add_dependencies(examples
  basic_files
  selecting_datasets
  std_vector_io
  complex_io
  writing_image
  append_scalar_data
  append_vector_data
  )

#
# build this only for new versions of HDF5
#

add_subdirectory(attributes)

if (H5CPP_SWMR_ENABLED)
  add_subdirectory(swmr)
endif ()

if (H5CPP_MPI_ENABLED)
  add_subdirectory(mpi)
endif ()

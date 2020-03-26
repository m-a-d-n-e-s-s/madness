# -*- mode: cmake -*-

######################
# Find Elemental
######################

if(ENABLE_ELEMENTAL AND DEFINED ELEMENTAL_TAG)

  include(ExternalProject)
  include(ConvertIncludesListToCompilerArgs)
  include(ConvertLibrariesListToCompilerArgs)

  find_package(Git REQUIRED)
  
  if(NOT DEFINED ELEMENTAL_URL)
    set(ELEMENTAL_URL https://github.com/elemental/Elemental.git)
  endif()
  message(STATUS "Will pull Elemental from ${ELEMENTAL_URL}")
  
  # Create a cache entry for Elemental build variables.
  # Note: This will not overwrite user specified values.
  set(ELEMENTAL_SOURCE_DIR "${PROJECT_BINARY_DIR}/external/source/elemental" CACHE PATH 
        "Path to the Elemental source directory")
  set(ELEMENTAL_BINARY_DIR "${PROJECT_BINARY_DIR}/external/build/elemental" CACHE PATH 
        "Path to the Elemental build directory")
  
  # Set Elemental compile flags
  append_flags(ELEMENTAL_CFLAGS "${CMAKE_C_FLAGS}")
  append_flags(ELEMENTAL_CXXFLAGS "${CMAKE_CXX_FLAGS}")
  append_flags(ELEMENTAL_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  
  if(CMAKE_BUILD_TYPE)
    string(TOLOWER ELEMENTAL_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
    append_flags(ELEMENTAL_CFLAGS "${CMAKE_C_FLAGS_${ELEMENTAL_BUILD_TYPE}}")
    append_flags(ELEMENTAL_CXXFLAGS "${CMAKE_CXX_FLAGS_${ELEMENTAL_BUILD_TYPE}}")
  endif()

  # Use same build type for Elemental. Assume this is newer Elemental.
  # Older El required decorated build types (PureRelease, HybridDebug, etc.)
  if (NOT DEFINED ELEMENTAL_CMAKE_BUILD_TYPE)
    set(ELEMENTAL_CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
  endif (NOT DEFINED ELEMENTAL_CMAKE_BUILD_TYPE)
  if (${ELEMENTAL_CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message(WARNING "Compiling Elemental with in Debug mode will make many calls non-reentrant! Set ELEMENTAL_CMAKE_BUILD_TYPE=Release to avoid")
  endif()
  
  # Set the configuration variables used by elemental
  if((ENABLE_SPINLOCKS OR NOT ENABLE_NEVER_SPIN) AND NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(ELEMENTAL_HAVE_SPINLOCKS ON)
  else()
    set(ELEMENTAL_HAVE_SPINLOCKS OFF)
  endif()

  # Override BLAS+LAPACK selection by Elemental
  # unless ELEMENTAL_MATH_LIBS is given by the user, use LAPACK_LIBRARIES as the default value for Elemental's MATH_LIBS
  if (NOT ELEMENTAL_MATH_LIBS AND NOT ("${ELEMENTAL_CMAKE_EXTRA_ARGS}" MATCHES "-DMATH_LIBS"))
    # see lapack.cmake
    # process LAPACK_LIBRARIES for CMAKE_REQUIRED_LIBRARIES (this is likely only to work with Makefile generator):
    # 1. get rid of the surrounding quotes
    string(REGEX REPLACE "\"" "" PROCESSED_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
    # 2. convert a space-separated string of libs into a list
    string(REGEX REPLACE " " ";" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
    # 3. restore (and protect!) the space in "-framework X"
    string(REGEX REPLACE "-framework;(.*)" "-framework \\1" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
    # 4. list -> string
    string(REGEX REPLACE ";" "\ " PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
    set(ELEMENTAL_MATH_LIBS "${PROCESSED_LAPACK_LIBRARIES}")
  endif()

  #
  # Obtain Elemental source **only** if needed (that's why not using ExternalProject)
  #
  # make directory
  message(STATUS "Checking Elemental source directory: ${ELEMENTAL_SOURCE_DIR}")
  if(EXISTS "${ELEMENTAL_SOURCE_DIR}")
    message(STATUS "Checking Elemental source directory: ${ELEMENTAL_SOURCE_DIR} - found")
  else(EXISTS "${ELEMENTAL_SOURCE_DIR}")
    # Create the external source directory
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/external/source)
      set(error_code 1)
      execute_process(
          COMMAND "${CMAKE_COMMAND}" -E make_directory "${PROJECT_BINARY_DIR}/external/source"
          RESULT_VARIABLE error_code)
      if(error_code)
        message(FATAL_ERROR "Failed to create the external source directory in build tree.")
      endif()
    endif()
  endif()
  # checkout if needed
  if(NOT EXISTS ${ELEMENTAL_SOURCE_DIR}/.git)
    message(STATUS "Pulling Elemental from: ${ELEMENTAL_URL}")
    set(error_code 1)
    set(number_of_tries 0)
    while(error_code AND number_of_tries LESS 3)
      execute_process(
          COMMAND ${GIT_EXECUTABLE}
                  clone ${ELEMENTAL_URL} elemental
          WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/external/source
          OUTPUT_FILE git.log
          ERROR_FILE git.log
          RESULT_VARIABLE error_code)
      math(EXPR number_of_tries "${number_of_tries} + 1")
    endwhile()
    if(number_of_tries GREATER 1)
      message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
    endif()
    if(error_code)
      message(FATAL_ERROR "Failed to clone repository: '${ELEMENTAL_URL}'")
    endif()
  endif()
  # reset to the desired Elemental tag
  if(EXISTS ${ELEMENTAL_SOURCE_DIR}/.git)
    set(error_code 1)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" fetch
      COMMAND "${GIT_EXECUTABLE}" checkout ${ELEMENTAL_TAG}
      WORKING_DIRECTORY "${ELEMENTAL_SOURCE_DIR}"
      RESULT_VARIABLE error_code)
    if(error_code)
      message(FATAL_ERROR "Failed to checkout tag: '${ELEMENTAL_TAG}'")
    endif()        
  endif()
  
  # Create or clean the build directory
  if(EXISTS "${ELEMENTAL_BINARY_DIR}")
    set(error_code 1)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E remove -f "./*"
        WORKING_DIRECTORY ${ELEMENTAL_BINARY_DIR}
        RESULT_VARIABLE error_code)
    if(error_code)
      message(FATAL_ERROR "Failed to delete existing files in the Elemental build directory.")
    endif()
  else()
    set(error_code 1)
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${ELEMENTAL_BINARY_DIR}"
        RESULT_VARIABLE error_code)
    if(error_code)
      message(FATAL_ERROR "Failed to create the Elemental build directory.")
    endif()
  endif()

  # since 0.85 package name is 'El', before that it was 'elemental'
  # detect the version by searching the main header
  message(STATUS "Looking for the top Elemental header")
  if(EXISTS ${ELEMENTAL_SOURCE_DIR}/include/El.hpp)
    message(STATUS "Looking for the top Elemental header - found El.hpp")
    set(HAVE_EL_H ON CACHE INTERNAL "Have El.hpp" FORCE)
    set(ELEMENTAL_PACKAGE_NAME El CACHE INTERNAL "Elemental package name" FORCE)
    set(ELEMENTAL_CONFIG_NAME Elemental CACHE INTERNAL "Elemental configure file prefix" FORCE)
  elseif(EXISTS ${ELEMENTAL_SOURCE_DIR}/include/elemental.hpp)
    message(STATUS "Looking for the top Elemental header - found elemental.hpp")
    set(HAVE_ELEMENTAL_H ON CACHE INTERNAL "Have elemental.hpp" FORCE)
    set(ELEMENTAL_PACKAGE_NAME elemental CACHE INTERNAL "Elemental package name (component name is same)" FORCE)
    # EFV: fairly sure the assumptions about which components Elemental builds by default are valid for
    # recent versions only 
    message (FATAL_ERROR "ELEMENTAL_TAG is too old to handle")
  else()
    message(FATAL_ERROR "Looking for the top Elemental header - not found")
  endif()

  # if CMAKE_POSITION_INDEPENDENT_CODE is defined, pass it to El
  if (DEFINED CMAKE_POSITION_INDEPENDENT_CODE)
    list(APPEND ELEMENTAL_CMAKE_EXTRA_ARGS -DCMAKE_POSITION_INDEPENDENT_CODE=${CMAKE_POSITION_INDEPENDENT_CODE})
  endif ()

  # concat all arguments
  set(ELEMENTAL_ARGS       -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
          -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
          -DCMAKE_BUILD_TYPE=${ELEMENTAL_CMAKE_BUILD_TYPE}
          -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
          -DMPI_C_COMPILER=${MPI_C_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
          "-DC_FLAGS=${ELEMENTAL_CFLAGS}"
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          "-DCXX_FLAGS=${ELEMENTAL_CXXFLAGS}"
          -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
          -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
          -DCMAKE_EXE_LINKER_FLAGS=${ELEMENTAL_LDFLAGS}
          -DMATH_LIBS=${ELEMENTAL_MATH_LIBS}
          -DHAVE_SPINLOCKS=${ELEMENTAL_HAVE_SPINLOCKS}
          ${ELEMENTAL_CMAKE_EXTRA_ARGS})
  message(STATUS "Elemental CMake arguments: ${ELEMENTAL_ARGS}")

  set(error_code 1)
  message (STATUS "** Configuring Elemental")
  execute_process(
      COMMAND ${CMAKE_COMMAND}
      ${ELEMENTAL_SOURCE_DIR}
      -G "${CMAKE_GENERATOR}"
      ${ELEMENTAL_ARGS}
      WORKING_DIRECTORY "${ELEMENTAL_BINARY_DIR}"
      RESULT_VARIABLE error_code)
  if(error_code)
    message(FATAL_ERROR "The Elemental cmake configuration failed.")
  else(error_code)
    message (STATUS "** Done configuring Elemental")
  endif(error_code)

  file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/CMake)

  find_package(${ELEMENTAL_CONFIG_NAME} REQUIRED PATHS ${CMAKE_INSTALL_PREFIX}
      COMPONENTS REQUIRED ${ELEMENTAL_PACKAGE_NAME} pmrrr ElSuiteSparse)
  set(ELEMENTAL_FOUND 1)

  # Add update-elemental target that will pull updates to the Elemental source
  # from the git repository. This is done outside ExternalProject_add to prevent
  # Elemental from doing a full pull, configure, and build everytime the project
  # is built.
  add_custom_target(update-elemental
    COMMAND ${GIT_EXECUTABLE} pull --rebase origin master
    COMMAND ${CMAKE_COMMAND} -E touch_nocreate ${ELEMENTAL_BINARY_DIR}/stamp/elemental-configure
    WORKING_DIRECTORY ${ELEMENTAL_SOURCE_DIR}
    COMMENT "Updating source for 'elemental' from ${ELEMENTAL_URL}"
    USES_TERMINAL)

  add_custom_target(build-elemental ALL
      COMMAND ${CMAKE_COMMAND} --build . --target El
      WORKING_DIRECTORY ${ELEMENTAL_BINARY_DIR}
      COMMENT "Building elemental"
      USES_TERMINAL)

  # TODO update Elemental installation when fixed upstream (installing El does not install everything it needs)
  add_custom_target(install-elemental
      COMMAND ${CMAKE_COMMAND} -P ${ELEMENTAL_BINARY_DIR}/cmake_install.cmake
      COMMENT "Installing elemental"
      USES_TERMINAL)
  add_dependencies(install-elemental build-elemental)

  # Add clean-elemental target that will delete files generated by Elemental build.
  add_custom_target(clean-elemental
    COMMAND $(MAKE) clean
    COMMAND ${CMAKE_COMMAND} -E touch_nocreate ${ELEMENTAL_BINARY_DIR}/stamp/elemental-configure
    WORKING_DIRECTORY ${ELEMENTAL_BINARY_DIR}
    COMMENT Cleaning build directory for 'elemental'
    USES_TERMINAL)

  # Since 'install-elemental' target cannot be linked to the 'install' target,
  # we will do it manually here.
  install(CODE
      "
      execute_process(
          COMMAND \"${CMAKE_MAKE_PROGRAM}\" \"install\" 
          WORKING_DIRECTORY \"${ELEMENTAL_BINARY_DIR}\"
          RESULT_VARIABLE error_code)
      if(error_code)
        message(FATAL_ERROR \"Failed to install 'elemental'\")
      endif()
      "
      )
  
  # Set build dependencies and compiler arguments
  # WARNING: should only need to satisfy build-elemental to use El, however its cmake config file
  # is inadequate for using from build tree, hence MUST install even for build tree. Solving here will
  # take far too much effort. Will file an issue.
  add_dependencies(${ELEMENTAL_PACKAGE_NAME} install-elemental)
  
  # These build variables are not used anyway
  set(ELEMENTAL_INCLUDE_DIRS
      "${ELEMENTAL_BINARY_DIR}/include"
     )
  set(ELEMENTAL_LIBRARIES 
      "${ELEMENTAL_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${ELEMENTAL_PACKAGE_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${ELEMENTAL_BINARY_DIR}/external/pmrrr/${CMAKE_STATIC_LIBRARY_PREFIX}pmrrr${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${ELEMENTAL_BINARY_DIR}/external/suite_sparse/${CMAKE_STATIC_LIBRARY_PREFIX}ElSuiteSparse${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${ELEMENTAL_BINARY_DIR}/download/parmetis/build/libparmetis/${CMAKE_STATIC_LIBRARY_PREFIX}parmetis${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${ELEMENTAL_BINARY_DIR}/download/parmetis/build/metis/libmetis/${CMAKE_STATIC_LIBRARY_PREFIX}parmetis${CMAKE_STATIC_LIBRARY_SUFFIX}"
      "${ELEMENTAL_BINARY_DIR}/download/scalapack/build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}scalapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
     )

endif(ENABLE_ELEMENTAL AND DEFINED ELEMENTAL_TAG)

# Find BLAS and LAPACK.

set(lapack_is_optional 0)
set(missing_lapack_message_level "FATAL_ERROR")

# if building just the runtime, included this by mistake, warn and make optional
if (MADNESS_BUILD_MADWORLD_ONLY)
  message(WARNING "MADNESS_BUILD_MADWORLD_ONLY=ON, but included external/lapack.cmake; must be error in CMakeLists.txt")
  set(lapack_is_optional 1)
  set(missing_lapack_message_level "STATUS")
endif (MADNESS_BUILD_MADWORLD_ONLY)

include(CheckCFortranFunctionExists)
include(CMakePushCheckState)
include(CheckCSourceRuns)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckFunctionExists)

# Helper function: Print unified, detailed error message with installation recommendations across OS/architectures
function(madness_report_lapack_error _reason)
  set(_help_msg "\n================================================================================")
  string(APPEND _help_msg "\nMADNESS Linear Algebra / LAPACK Configuration Error:")
  string(APPEND _help_msg "  ${_reason}")
  string(APPEND _help_msg "\n================================================================================\n")

  string(APPEND _help_msg "MADNESS relies heavily on fine-grained task-based multithreading with lock-free\n")
  string(APPEND _help_msg "task queues and work stealing. OpenBLAS contains internal threading locks and\n")
  string(APPEND _help_msg "thread-local state that cause severe concurrency bugs, memory corruption, and\n")
  string(APPEND _help_msg "deadlocks under MADNESS multi-threaded execution (especially on ARM architectures).\n")
  string(APPEND _help_msg "Therefore, OpenBLAS is rejected by default.\n")

  string(APPEND _help_msg "RECOMMENDED LINEAR ALGEBRA PACKAGES BY ARCHITECTURE:")
  string(APPEND _help_msg "\n--------------------------------------------------------------------------------\n")

  string(APPEND _help_msg "[1] ARM Architecture (AArch64 / Apple Silicon):\n")
  string(APPEND _help_msg "    Preference Hierarchy: (1) ARMPL  ->  (2) NVPL  ->  (3) BLIS (serial) + Netlib LAPACK\n")
  string(APPEND _help_msg "    * macOS (Apple Silicon M1/M2/M3/M4):\n")
  string(APPEND _help_msg "      - Apple Accelerate is built-in and supported automatically:\n")
  string(APPEND _help_msg "          (No installation needed; CMake detects -framework Accelerate)\n")
  string(APPEND _help_msg "      - Or install BLIS & Netlib LAPACK via Homebrew:\n")
  string(APPEND _help_msg "          brew install blis lapack\n")
  string(APPEND _help_msg "    * Ubuntu / Debian (ARM64):\n")
  string(APPEND _help_msg "      - (1) Arm Performance Libraries (ARMPL, Best Performance):\n")
  string(APPEND _help_msg "          Download DEB from: https://developer.arm.com/downloads/-/arm-performance-libraries\n")
  string(APPEND _help_msg "          Install: sudo dpkg -i arm-performance-libraries_*.deb\n")
  string(APPEND _help_msg "          (Default location: /opt/arm/arm-performance-libraries)\n")
  string(APPEND _help_msg "      - (2) NVIDIA Performance Libraries (NVPL, for Grace/ARM):\n")
  string(APPEND _help_msg "          sudo apt-get install nvpl-blas nvpl-lapack\n")
  string(APPEND _help_msg "          (Or available via NVIDIA HPC SDK)\n")
  string(APPEND _help_msg "      - (3) Serial BLIS + Netlib LAPACK (Safe Fallback):\n")
  string(APPEND _help_msg "          sudo apt-get install libblis4-serial liblapack-dev liblapack3\n")
  string(APPEND _help_msg "          Note: Avoid installing libopenblas-dev as it overrides system BLAS alternatives.\n")
  string(APPEND _help_msg "    * Fedora / RHEL / CentOS / Rocky Linux (AArch64):\n")
  string(APPEND _help_msg "      - (1) Arm Performance Libraries (ARMPL, Best Performance):\n")
  string(APPEND _help_msg "          Download RPM from: https://developer.arm.com/downloads/-/arm-performance-libraries\n")
  string(APPEND _help_msg "          Install: sudo dnf install arm-performance-libraries_*.rpm\n")
  string(APPEND _help_msg "      - (2) NVIDIA Performance Libraries (NVPL):\n")
  string(APPEND _help_msg "          sudo dnf install nvpl-blas nvpl-lapack\n")
  string(APPEND _help_msg "      - (3) Serial BLIS + Netlib LAPACK (Safe Fallback):\n")
  string(APPEND _help_msg "          sudo dnf install blis-devel lapack-devel\n")

  string(APPEND _help_msg "[2] x86_64 Architecture (Intel / AMD):\n")
  string(APPEND _help_msg "    Preference Hierarchy: (1) Intel oneAPI MKL (Strongly Preferred)  ->  (2) BLIS (serial) + Netlib LAPACK\n")
  string(APPEND _help_msg "    * Ubuntu / Debian (x86_64):\n")
  string(APPEND _help_msg "      - (1) Intel oneAPI MKL (Strongly Preferred for performance):\n")
  string(APPEND _help_msg "          sudo apt-get install intel-oneapi-mkl-devel\n")
  string(APPEND _help_msg "          (Or install Intel oneAPI Base Toolkit, or: conda install mkl-devel)\n")
  string(APPEND _help_msg "      - (2) Serial BLIS + Netlib LAPACK (Safe Fallback):\n")
  string(APPEND _help_msg "          sudo apt-get install libblis4-serial liblapack-dev liblapack3\n")
  string(APPEND _help_msg "    * Fedora / RHEL / CentOS / Rocky Linux (x86_64):\n")
  string(APPEND _help_msg "      - (1) Intel oneAPI MKL (Strongly Preferred for performance):\n")
  string(APPEND _help_msg "          sudo dnf install intel-oneapi-mkl-devel\n")
  string(APPEND _help_msg "      - (2) Serial BLIS + Netlib LAPACK (Safe Fallback):\n")
  string(APPEND _help_msg "          sudo dnf install blis-devel lapack-devel\n")
  string(APPEND _help_msg "    * macOS (Intel x86_64):\n")
  string(APPEND _help_msg "      - Apple Accelerate framework (built-in, automatically detected)\n")
  string(APPEND _help_msg "      - Or Intel MKL via Conda (conda install mkl-devel) or Intel oneAPI\n")

  string(APPEND _help_msg "MANUAL CONFIGURATION / OVERRIDE:\n")
  string(APPEND _help_msg "--------------------------------------------------------------------------------\n")
  string(APPEND _help_msg "You can always manually specify your BLAS and LAPACK libraries by passing\n")
  string(APPEND _help_msg "-DLAPACK_LIBRARIES (and optionally -DLAPACK_INCLUDE_DIRS) to CMake:\n")
  string(APPEND _help_msg "  cmake .. -DLAPACK_LIBRARIES=\"/path/to/liblapack.a;/path/to/libblas.a\" \\\n")
  string(APPEND _help_msg "           -DLAPACK_INCLUDE_DIRS=\"/path/to/include\"\n\n")
  string(APPEND _help_msg "To explicitly override OpenBLAS safety checks and force OpenBLAS at your own risk:\n")
  string(APPEND _help_msg "  cmake .. -DENABLE_OPENBLAS=ON\n")
  string(APPEND _help_msg "================================================================================\n")

  message("${missing_lapack_message_level}" "${_help_msg}")
endfunction()

# Helper function: Check if a library or list of libraries resolves to or dynamically loads OpenBLAS
function(madness_check_is_openblas _lib_list _out_var)
  set(${_out_var} FALSE PARENT_SCOPE)
  if(NOT _lib_list)
    return()
  endif()

  # 1. Path check (including resolved symlinks)
  foreach(_lib ${_lib_list})
    if(EXISTS "${_lib}")
      get_filename_component(_real_path "${_lib}" REALPATH)
      if(_real_path MATCHES "openblas" OR _lib MATCHES "openblas")
        set(${_out_var} TRUE PARENT_SCOPE)
        return()
      endif()
      # On Linux, inspect ELF dynamic dependencies (DT_NEEDED) using readelf if available
      if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND _lib MATCHES "\\.so")
        find_program(READELF_EXECUTABLE NAMES readelf)
        if(READELF_EXECUTABLE)
          execute_process(
            COMMAND ${READELF_EXECUTABLE} -d "${_lib}"
            OUTPUT_VARIABLE _readelf_out
            ERROR_QUIET
          )
          if(_readelf_out MATCHES "NEEDED.*openblas")
            set(${_out_var} TRUE PARENT_SCOPE)
            return()
          endif()
          if(_readelf_out MATCHES "NEEDED.*libblas\\.so")
            find_file(_sys_libblas NAMES libblas.so.3 libblas.so
              PATHS /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE} /usr/lib /lib/${CMAKE_LIBRARY_ARCHITECTURE} /lib
              NO_DEFAULT_PATH
            )
            if(_sys_libblas)
              get_filename_component(_sys_libblas_real "${_sys_libblas}" REALPATH)
              if(_sys_libblas_real MATCHES "openblas")
                set(${_out_var} TRUE PARENT_SCOPE)
                return()
              endif()
            endif()
          endif()
        endif()
      endif()
    elseif(_lib MATCHES "openblas")
      set(${_out_var} TRUE PARENT_SCOPE)
      return()
    endif()
  endforeach()

  # 2. Static symbol check by test compilation/linking
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES ${_lib_list})
  unset(_has_openblas_symbol CACHE)
  check_function_exists(openblas_get_config _has_openblas_symbol)
  cmake_pop_check_state()
  if(_has_openblas_symbol)
    set(${_out_var} TRUE PARENT_SCOPE)
    return()
  endif()
  unset(_has_openblas_symbol CACHE)

  # 3. Dynamic runtime check (checks if OpenBLAS was loaded into memory at runtime)
  if(NOT CMAKE_CROSSCOMPILING)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES ${_lib_list})
    find_package(Threads QUIET)
    if(TARGET Threads::Threads)
      list(APPEND CMAKE_REQUIRED_LIBRARIES Threads::Threads)
    endif()
    find_library(_m_chk_lib NAMES m)
    if(_m_chk_lib)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${_m_chk_lib}")
    endif()
    find_library(_gfort_chk_lib NAMES gfortran)
    if(_gfort_chk_lib)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${_gfort_chk_lib}")
    endif()
    unset(_CHECK_DL_NO_OPENBLAS CACHE)
    check_c_source_runs(
      "
      #define _GNU_SOURCE
      #include <stdio.h>
      #include <string.h>
      #if defined(__linux__)
      #include <link.h>
      static int callback(struct dl_phdr_info *info, size_t size, void *data) {
          if (info->dlpi_name != NULL && info->dlpi_name[0] != 0) {
              if (strstr(info->dlpi_name, \"openblas\") != NULL) {
                  *((int*)data) = 1;
              }
          }
          return 0;
      }
      #endif
      extern void sgemm_(const char*, const char*, const int*, const int*, const int*, const float*, const float*, const int*, const float*, const int*, const float*, float*, const int*);
      int main(void) {
          int has_openblas = 0;
      #if defined(__linux__)
          char trans = 'N';
          int m = 0, n = 0, k = 0, lda = 1, ldb = 1, ldc = 1;
          float alpha = 1.0f, beta = 0.0f, a[1] = {0}, b[1] = {0}, c[1] = {0};
          sgemm_(&trans, &trans, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
          dl_iterate_phdr(callback, &has_openblas);
      #endif
          return (has_openblas == 0) ? 0 : 1;
      }
      " _CHECK_DL_NO_OPENBLAS
    )
    cmake_pop_check_state()
    if(NOT _CHECK_DL_NO_OPENBLAS)
      set(${_out_var} TRUE PARENT_SCOPE)
      return()
    endif()
    unset(_CHECK_DL_NO_OPENBLAS CACHE)
  endif()
endfunction()

set(LAPACK_FOUND FALSE)
set(HAVE_INTEL_MKL 0)
set(HAVE_ARMPL 0)
set(HAVE_NVPL 0)
set(HAVE_BLIS 0)
set(HAVE_ACML 0)

if(NOT LAPACK_LIBRARIES)
  set(USER_LAPACK_LIBRARIES FALSE)

  # 1. Search for Intel MKL
  if(ENABLE_MKL)
    find_package(MKL)
    if(MKL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
      set(HAVE_INTEL_MKL 1)
      set(LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE)
      set(LAPACK_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
    endif()
  endif()

  # 2. Search for Arm Performance Libraries (ARMPL) on ARM
  if(ENABLE_ARMPL AND NOT LAPACK_FOUND)
    find_package(ARMPL)
    if(ARMPL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ARMPL_LIBRARIES})
      set(LAPACK_INCLUDE_DIRS ${ARMPL_INCLUDE_DIRS})
      set(HAVE_ARMPL 1)
    endif()
  endif()

  # 3. Search for NVIDIA Performance Libraries (NVPL) on ARM
  if(ENABLE_NVPL AND NOT LAPACK_FOUND)
    find_package(NVPL)
    if(NVPL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${NVPL_LIBRARIES})
      set(LAPACK_INCLUDE_DIRS ${NVPL_INCLUDE_DIRS})
      set(HAVE_NVPL 1)
    endif()
  endif()

  # 4. Search for BLIS (serial) + LAPACK (Netlib)
  if(ENABLE_BLIS AND NOT LAPACK_FOUND)
    find_package(BLIS)
    if(BLIS_FOUND)
      # BLIS provides BLAS; find a compatible LAPACK library (e.g. Netlib LAPACK)
      # Search specifically in Netlib directories first to avoid picking up OpenBLAS symlinks
      set(_lapack_hints
        ${LAPACK_ROOT}
        ${LAPACK_DIR}
        ${LAPACK_ROOT_DIR}
        $ENV{LAPACK_DIR}
        $ENV{LAPACK_ROOT}
        ${BLIS_ROOT_DIR}
        ${BLIS_DIR}
        ${BLIS_ROOT}
        /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/lapack
        /usr/lib/aarch64-linux-gnu/lapack
        /usr/lib/x86_64-linux-gnu/lapack
        /usr/lib/lapack
        /usr/lib64/lapack
        /usr/lib64
        /usr/local/lib64
        /usr/local/lib/lapack
        /usr/local/lib
        /usr/lib
      )
      # Search for static liblapack.a first (prevents DT_NEEDED libblas.so.3 alternative poisoning on Ubuntu/Debian)
      find_library(BLIS_LAPACK_LIBRARY
        NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}" lapack
        HINTS ${_lapack_hints}
        PATH_SUFFIXES lapack
      )
      if(NOT BLIS_LAPACK_LIBRARY)
        find_library(BLIS_LAPACK_LIBRARY
          NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}" lapack
          HINTS ${_lapack_hints}
        )
      endif()

      # Check if pairing with BLIS_LAPACK_LIBRARY pulls in OpenBLAS
      if(BLIS_LAPACK_LIBRARY AND NOT ENABLE_OPENBLAS)
        madness_check_is_openblas("${BLIS_LAPACK_LIBRARY};${BLIS_LIBRARIES}" _blis_lapack_is_openblas)
        if(_blis_lapack_is_openblas)
          # Try falling back specifically to static liblapack.a
          unset(BLIS_LAPACK_LIBRARY CACHE)
          unset(BLIS_LAPACK_LIBRARY)
          find_library(BLIS_LAPACK_LIBRARY
            NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
            HINTS ${_lapack_hints}
            PATH_SUFFIXES lapack
          )
          if(BLIS_LAPACK_LIBRARY)
            madness_check_is_openblas("${BLIS_LAPACK_LIBRARY};${BLIS_LIBRARIES}" _blis_lapack_is_ob2)
            if(_blis_lapack_is_ob2)
              message(STATUS "Rejecting LAPACK library that pulls in OpenBLAS: ${BLIS_LAPACK_LIBRARY}")
              unset(BLIS_LAPACK_LIBRARY CACHE)
              unset(BLIS_LAPACK_LIBRARY)
            endif()
          endif()
        endif()
      endif()

      if(BLIS_LAPACK_LIBRARY)
        set(LAPACK_FOUND TRUE)
        set(_extra_lapack_libs)
        if(BLIS_LAPACK_LIBRARY MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}$")
          find_library(_gfortran_lib NAMES gfortran)
          if(_gfortran_lib)
            list(APPEND _extra_lapack_libs "${_gfortran_lib}")
          else()
            list(APPEND _extra_lapack_libs "-lgfortran")
          endif()
        endif()
        set(LAPACK_LIBRARIES ${BLIS_LAPACK_LIBRARY} ${BLIS_LIBRARIES} ${_extra_lapack_libs})
        set(LAPACK_INCLUDE_DIRS ${BLIS_INCLUDE_DIRS})
        set(HAVE_BLIS 1)
      else()
        message(STATUS "BLIS BLAS found (${BLIS_LIBRARY}), but accompanying compatible Netlib LAPACK library was not found.")
      endif()
    endif()
  endif()

  # 5. Search for ACML
  if(ENABLE_ACML AND NOT LAPACK_FOUND)
    find_package(ACML)
    if(ACML_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ACML_LIBRARIES})
      set(HAVE_ACML 1)
    endif()
  endif()

  # 6. Search for system specific BLAS/LAPACK checks (Darwin / Accelerate)
  if(NOT LAPACK_FOUND AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # Accelerate is always present, so no need to search
    set(LAPACK_LIBRARIES "-framework Accelerate")
    set(LAPACK_FOUND TRUE)
  endif()

  # 7. Search for Netlib lapack and blas libraries
  if(NOT LAPACK_FOUND)
    set(_netlib_hints
      ${LAPACK_ROOT}
      ${LAPACK_DIR}
      ${LAPACK_ROOT_DIR}
      $ENV{LAPACK_DIR}
      $ENV{LAPACK_ROOT}
      ${BLAS_ROOT}
      ${BLAS_DIR}
      ${BLAS_ROOT_DIR}
      $ENV{BLAS_DIR}
      $ENV{BLAS_ROOT}
      /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/lapack
      /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/blas
      /usr/lib/aarch64-linux-gnu/lapack
      /usr/lib/aarch64-linux-gnu/blas
      /usr/lib/x86_64-linux-gnu/lapack
      /usr/lib/x86_64-linux-gnu/blas
      /usr/lib/lapack
      /usr/lib/blas
      /usr/lib64/lapack
      /usr/lib64/blas
      /usr/lib64
      /usr/local/lib64
      /usr/local/lib/lapack
      /usr/local/lib/blas
      /usr/local/lib
      /usr/lib
    )
    find_library(LAPACK_lapack_LIBRARY
      NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}" lapack
      HINTS ${_netlib_hints}
      PATH_SUFFIXES lapack
    )
    if(NOT LAPACK_lapack_LIBRARY)
      find_library(LAPACK_lapack_LIBRARY
        NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}" lapack
        HINTS ${_netlib_hints}
      )
    endif()
    find_library(LAPACK_blas_LIBRARY
      NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}" blas
      HINTS ${_netlib_hints}
      PATH_SUFFIXES blas
    )
    if(NOT LAPACK_blas_LIBRARY)
      find_library(LAPACK_blas_LIBRARY
        NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}" blas
        HINTS ${_netlib_hints}
      )
    endif()

    if(LAPACK_lapack_LIBRARY AND NOT ENABLE_OPENBLAS)
      madness_check_is_openblas("${LAPACK_lapack_LIBRARY}" _lapack_is_ob)
      if(_lapack_is_ob)
        unset(LAPACK_lapack_LIBRARY CACHE)
        unset(LAPACK_lapack_LIBRARY)
        find_library(LAPACK_lapack_LIBRARY
          NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}"
          HINTS ${_netlib_hints}
          PATH_SUFFIXES lapack
        )
        if(LAPACK_lapack_LIBRARY)
          madness_check_is_openblas("${LAPACK_lapack_LIBRARY}" _lapack_is_ob2)
          if(_lapack_is_ob2)
            message(STATUS "Rejecting OpenBLAS library found for LAPACK: ${LAPACK_lapack_LIBRARY}")
            unset(LAPACK_lapack_LIBRARY CACHE)
            unset(LAPACK_lapack_LIBRARY)
          endif()
        endif()
      endif()
    endif()

    if(LAPACK_blas_LIBRARY AND NOT ENABLE_OPENBLAS)
      madness_check_is_openblas("${LAPACK_blas_LIBRARY}" _blas_is_ob)
      if(_blas_is_ob)
        unset(LAPACK_blas_LIBRARY CACHE)
        unset(LAPACK_blas_LIBRARY)
        find_library(LAPACK_blas_LIBRARY
          NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}"
          HINTS ${_netlib_hints}
          PATH_SUFFIXES blas
        )
        if(LAPACK_blas_LIBRARY)
          madness_check_is_openblas("${LAPACK_blas_LIBRARY}" _blas_is_ob2)
          if(_blas_is_ob2)
            message(STATUS "Rejecting OpenBLAS library found for BLAS: ${LAPACK_blas_LIBRARY}")
            unset(LAPACK_blas_LIBRARY CACHE)
            unset(LAPACK_blas_LIBRARY)
          endif()
        endif()
      endif()
    endif()

    if(LAPACK_lapack_LIBRARY AND LAPACK_blas_LIBRARY)
      madness_check_is_openblas("${LAPACK_lapack_LIBRARY};${LAPACK_blas_LIBRARY}" _pair_is_ob)
      if(_pair_is_ob AND NOT ENABLE_OPENBLAS)
        message(STATUS "Rejecting LAPACK+BLAS combination that pulls in OpenBLAS: ${LAPACK_lapack_LIBRARY} ${LAPACK_blas_LIBRARY}")
      else()
        set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARY} ${LAPACK_blas_LIBRARY})
        if(LAPACK_lapack_LIBRARY MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}$" OR LAPACK_blas_LIBRARY MATCHES "${CMAKE_STATIC_LIBRARY_SUFFIX}$")
          find_library(_gfortran_lib NAMES gfortran)
          if(_gfortran_lib)
            list(APPEND LAPACK_LIBRARIES "${_gfortran_lib}")
          endif()
        endif()
        set(LAPACK_FOUND TRUE)
      endif()
    endif()
  endif()

  # 8. Fallback to OpenBLAS only if explicitly requested via ENABLE_OPENBLAS=ON
  if(ENABLE_OPENBLAS AND NOT LAPACK_FOUND)
    find_library(LAPACK_openblas_LIBRARY NAMES openblas)
    if(LAPACK_openblas_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_openblas_LIBRARY})
      set(LAPACK_FOUND TRUE)
    endif()
  endif()

else()
  set(USER_LAPACK_LIBRARIES TRUE)
endif()

if(NOT LAPACK_LIBRARIES)
  madness_report_lapack_error("No suitable BLAS and LAPACK libraries were found automatically on your system.")
  return()
endif()

cmake_push_check_state()

# process LAPACK_LIBRARIES for CMAKE_REQUIRED_LIBRARIES
string(REGEX REPLACE "\"" "" PROCESSED_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
string(REGEX REPLACE " " ";" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
string(REGEX REPLACE "-framework;(.*)" "-framework\\\\ \\1" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")

set(CMAKE_REQUIRED_LIBRARIES ${PROCESSED_LAPACK_LIBRARIES} ${CMAKE_REQUIRED_LIBRARIES}
        Threads::Threads)
if(LAPACK_INCLUDE_DIRS)
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS} ${CMAKE_REQUIRED_INCLUDES})
endif()

# Verify that we can link against BLAS
check_c_fortran_function_exists(sgemm BLAS_WORKS)

if(BLAS_WORKS)
  message(STATUS "A library with BLAS API found.")
else()
  madness_report_lapack_error("Unable to link against BLAS function 'sgemm' using specified libraries: ${LAPACK_LIBRARIES}")
endif()

# Verify that we can link against LAPACK
check_c_fortran_function_exists(cheev LAPACK_WORKS)

if(LAPACK_WORKS)
  message(STATUS "A library with LAPACK API found.")
else()
  madness_report_lapack_error("Unable to link against LAPACK function 'cheev' using specified libraries: ${LAPACK_LIBRARIES}")
endif()

# Check for OpenBLAS rejection/warning
if(NOT ENABLE_OPENBLAS)
  madness_check_is_openblas("${LAPACK_LIBRARIES}" LAPACK_IS_OPENBLAS)
  if(LAPACK_IS_OPENBLAS)
    if(USER_LAPACK_LIBRARIES)
      message(WARNING "OpenBLAS detected in user-specified LAPACK_LIBRARIES. Note: OpenBLAS is not safe for multithreaded concurrent calls on ARM.")
    else()
      madness_report_lapack_error("Auto-detected linear algebra library '${LAPACK_LIBRARIES}' is OpenBLAS (or dynamically depends on OpenBLAS), which is unsafe for multithreading on ARM.")
    endif()
  endif()
endif()

set(LAPACK_FOUND TRUE)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

# Introspect LAPACK_LIBRARIES (both user-specified and auto-detected)
unset(USER_LAPACK_LIBRARIES_IS_MKL CACHE)
check_function_exists(mkl_get_version USER_LAPACK_LIBRARIES_IS_MKL)
if(USER_LAPACK_LIBRARIES_IS_MKL)
  message(STATUS "LAPACK provides an MKL library")
  set(HAVE_INTEL_MKL 1)
  list(APPEND LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE)
  list(REMOVE_DUPLICATES LAPACK_COMPILE_DEFINITIONS)
endif()

unset(USER_LAPACK_LIBRARIES_IS_ARMPL CACHE)
check_function_exists(armplversion USER_LAPACK_LIBRARIES_IS_ARMPL)
if(USER_LAPACK_LIBRARIES_IS_ARMPL)
  message(STATUS "LAPACK provides an ARMPL library")
  set(HAVE_ARMPL 1)
endif()

unset(USER_LAPACK_LIBRARIES_IS_NVPL CACHE)
check_function_exists(nvpl_blas_get_version USER_LAPACK_LIBRARIES_IS_NVPL)
if(USER_LAPACK_LIBRARIES_IS_NVPL)
  message(STATUS "LAPACK provides an NVPL library")
  set(HAVE_NVPL 1)
endif()

unset(USER_LAPACK_LIBRARIES_IS_BLIS CACHE)
check_function_exists(bli_info_get_version_str USER_LAPACK_LIBRARIES_IS_BLIS)
if(USER_LAPACK_LIBRARIES_IS_BLIS)
  message(STATUS "LAPACK provides a BLIS library")
  set(HAVE_BLIS 1)
  message(WARNING
    "BLIS was selected as the BLAS library. Note that BLIS performance for multi-threaded small matrix contractions (critical in MADNESS MRA multiwavelet operations) is significantly slower than vendor-tuned libraries due to framework overhead and memory pool lock contention. For optimal performance, it is strongly recommended to install Intel MKL (on x86_64) or Arm Performance Libraries (ARMPL, on ARM64) if at all possible."
  )
endif()

unset(USER_LAPACK_LIBRARIES_IS_ACML CACHE)
check_function_exists(acmlversion USER_LAPACK_LIBRARIES_IS_ACML)
if(USER_LAPACK_LIBRARIES_IS_ACML)
  message(STATUS "LAPACK provides an ACML library")
  set(HAVE_ACML 1)
endif()

# LAPACK_LIBRARIES might include IMPORTED targets that are not globally available
if(LAPACK_LIBRARIES MATCHES OpenMP::OpenMP_C AND NOT TARGET OpenMP::OpenMP_C)
  find_package(OpenMP REQUIRED COMPONENTS C)
endif()
if(LAPACK_LIBRARIES MATCHES Threads::Threads AND NOT TARGET Threads::Threads)
  find_package(Threads REQUIRED)
endif()

cmake_pop_check_state()

# Set the fortran mangling scheme.
if(LAPACK_WORKS STREQUAL "cheev_")
  set(FORTRAN_LINKAGE_LCU 1)
elseif(LAPACK_WORKS STREQUAL "cheev")
  set(FORTRAN_LINKAGE_LC 1)
elseif(LAPACK_WORKS STREQUAL "cheev__")
  set(FORTRAN_LINKAGE_LCUU 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV")
  set(FORTRAN_LINKAGE_UC 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV_")
  set(FORTRAN_LINKAGE_UCU 1)
endif()

# unquote LAPACK_COMPILE_OPTIONS, LAPACK_INCLUDE_DIRS, and LAPACK_COMPILE_DEFINITIONS
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_OPTIONS "${LAPACK_COMPILE_OPTIONS}")
string(REGEX REPLACE "\"" "" LAPACK_INCLUDE_DIRS "${LAPACK_INCLUDE_DIRS}")
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_DEFINITIONS "${LAPACK_COMPILE_DEFINITIONS}")

# epilogue: final sanity checks
if(USER_LAPACK_LIBRARIES_IS_MKL)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  foreach(def ${LAPACK_COMPILE_DEFINITIONS})
    set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS};-D${def}")
  endforeach()
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  unset(MADNESS_CAN_INCLUDE_MKL_H CACHE)
  check_cxx_source_compiles(
    "
    #include <mkl.h>
    int main(int argc, char** argv) {
      return 0;
    }
    " MADNESS_CAN_INCLUDE_MKL_H)
  if(NOT MADNESS_CAN_INCLUDE_MKL_H)
    madness_report_lapack_error("LAPACK provides MKL but cannot include its headers; ensure that corresponding LAPACK_INCLUDE_DIRS, LAPACK_COMPILE_DEFINITIONS, or LAPACK_COMPILE_OPTIONS were provided.")
  endif()
  cmake_pop_check_state()
endif(USER_LAPACK_LIBRARIES_IS_MKL)

if(USER_LAPACK_LIBRARIES_IS_ARMPL AND LAPACK_INCLUDE_DIRS)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  unset(MADNESS_CAN_INCLUDE_ARMPL_H CACHE)
  check_cxx_source_compiles(
    "
    #include <armpl.h>
    int main(int argc, char** argv) {
      return 0;
    }
    " MADNESS_CAN_INCLUDE_ARMPL_H)
  cmake_pop_check_state()
endif()

if(USER_LAPACK_LIBRARIES_IS_BLIS AND LAPACK_INCLUDE_DIRS)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  set(CMAKE_REQUIRED_LIBRARIES ${PROCESSED_LAPACK_LIBRARIES} Threads::Threads)
  unset(MADNESS_CAN_INCLUDE_BLIS_H CACHE)
  check_cxx_source_compiles(
    "
    #include <blis.h>
    int main(int argc, char** argv) {
      return 0;
    }
    " MADNESS_CAN_INCLUDE_BLIS_H)
  if(NOT MADNESS_CAN_INCLUDE_BLIS_H)
    madness_report_lapack_error("LAPACK provides BLIS but cannot include its headers; ensure that corresponding LAPACK_INCLUDE_DIRS, LAPACK_COMPILE_DEFINITIONS, or LAPACK_COMPILE_OPTIONS were provided.")
  endif()
  if(NOT CMAKE_CROSSCOMPILING)
    unset(MADNESS_BLIS_IS_SERIAL CACHE)
    check_cxx_source_runs(
      "
      #include <blis.h>
      int main(int argc, char** argv) {
        // Return 0 if single-threaded (serial), 1 if multithreaded
        return (bli_info_get_enable_threading() == 0) ? 0 : 1;
      }
      " MADNESS_BLIS_IS_SERIAL)
    if(NOT MADNESS_BLIS_IS_SERIAL)
      madness_report_lapack_error("The detected BLIS library is multithreaded (OpenMP/pthreads). MADNESS strictly requires thread-safe sequential (single-thread) BLIS.")
    endif()
  endif()
  cmake_pop_check_state()
endif()

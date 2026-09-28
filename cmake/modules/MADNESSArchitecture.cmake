include_guard(GLOBAL)

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

set(MADNESS_TARGET_ARCH "default" CACHE STRING
    "Target CPU architecture: 'default' (the toolchain's baseline, portable), 'performance' (x86-64-v3 on x86: AVX2/FMA/BMI2, Haswell/Excavator and newer), 'native' (the build host), 'none'/OFF, or a custom compiler arch/cpu string")

set(MADNESS_TUNE_ARCH "" CACHE STRING
    "Optional CPU scheduling tuning (maps to -mtune on x86, e.g., 'native' or 'zen4')")

set(MADNESS_RELAXED_MATH "strict" CACHE STRING
    "Floating-point optimization mode: 'strict'/OFF (IEEE-754, the default), 'relaxed' (associative + FMA contraction, preserves NaN/Inf), or 'fast' (-ffast-math)")

# Check if target architecture flags are disabled
if (MADNESS_TARGET_ARCH MATCHES "^(none|NONE|OFF|off|False|FALSE|0)$")
  set(_madness_apply_arch OFF)
else()
  set(_madness_apply_arch ON)
endif()

# 'default' adds no flags, so the binaries run wherever the toolchain's own output runs: distro and cluster builds
# are routinely compiled on a newer host than the one they run on, and a configure-time check cannot see that.
# An architecture the user already chose through the compiler flags (by hand or via a toolchain file) counts as
# a deliberate choice, so no warning about the portable baseline is due then.
if (_madness_apply_arch AND MADNESS_TARGET_ARCH STREQUAL "default")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _madness_build_type)
  set(_madness_user_flags "${CMAKE_C_FLAGS} ${CMAKE_CXX_FLAGS}")
  if (_madness_build_type)
    string(APPEND _madness_user_flags " ${CMAKE_C_FLAGS_${_madness_build_type}} ${CMAKE_CXX_FLAGS_${_madness_build_type}}")
  endif()
  if (_madness_user_flags MATCHES "(^| )-(march|mcpu)=")
    set(_madness_apply_arch OFF)
    set(_madness_arch_from_user_flags ON)
    message(STATUS "MADNESS: CMAKE_<LANG>_FLAGS select the target architecture")
  endif()
endif()

if (_madness_apply_arch)
  set(_arch_flags_to_apply "")

  if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64|i.86)$")
    if (MADNESS_TARGET_ARCH STREQUAL "default")
      message(WARNING "MADNESS: building for the portable x86-64 baseline (no AVX2/FMA), which runs MADNESS "
                      "roughly 15-20% slower than an AVX2 build. If every machine that will run this build has AVX2 "
                      "(Intel Haswell / AMD Excavator, 2013, or newer), configure with "
                      "-DMADNESS_TARGET_ARCH=performance, or with -DMADNESS_TARGET_ARCH=native to target this host only.")
    elseif (MADNESS_TARGET_ARCH STREQUAL "performance")
      check_cxx_compiler_flag("-march=x86-64-v3" _HAS_X86_64_V3)
      if (_HAS_X86_64_V3)
        list(APPEND _arch_flags_to_apply "-march=x86-64-v3")
      else()
        # compilers older than GCC 11 / Clang 12 do not know the x86-64-v3 level; spell out its main extensions
        list(APPEND _arch_flags_to_apply "-mavx2" "-mfma" "-mbmi2")
      endif()
    elseif (MADNESS_TARGET_ARCH STREQUAL "native")
      list(APPEND _arch_flags_to_apply "-march=native")
    else()
      # Custom user-supplied architecture name
      list(APPEND _arch_flags_to_apply "-march=${MADNESS_TARGET_ARCH}")
    endif()

    if (MADNESS_TUNE_ARCH)
      list(APPEND _arch_flags_to_apply "-mtune=${MADNESS_TUNE_ARCH}")
    endif()

  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    if (APPLE)
      # On macOS Apple Silicon, Apple Clang natively targets Darwin ARM64, whose baseline
      # includes Apple Silicon capabilities (ARMv8.5-A+). Forcing a generic Linux-style
      # -march=armv8-a flag is unnecessary and inappropriately downgrades the baseline.
      if (MADNESS_TARGET_ARCH STREQUAL "native")
        check_cxx_compiler_flag("-mcpu=native" _HAS_MCPU_NATIVE)
        if (_HAS_MCPU_NATIVE)
          list(APPEND _arch_flags_to_apply "-mcpu=native")
        endif()
      elseif (NOT (MADNESS_TARGET_ARCH STREQUAL "default" OR MADNESS_TARGET_ARCH STREQUAL "performance"))
        if (MADNESS_TARGET_ARCH MATCHES "^armv")
          list(APPEND _arch_flags_to_apply "-march=${MADNESS_TARGET_ARCH}")
        else()
          list(APPEND _arch_flags_to_apply "-mcpu=${MADNESS_TARGET_ARCH}")
        endif()
      endif()
    else()
      # 'default' keeps the compiler's own baseline, which is never below armv8-a and may be higher
      # (e.g. a toolchain configured --with-cpu=neoverse-v2). There is no portable ARM counterpart of x86-64-v3
      # worth selecting by default, so 'performance' does the same; use 'native' or an explicit CPU instead.
      if (MADNESS_TARGET_ARCH STREQUAL "default" OR MADNESS_TARGET_ARCH STREQUAL "performance")
      elseif (MADNESS_TARGET_ARCH STREQUAL "native")
        check_cxx_compiler_flag("-mcpu=native" _HAS_MCPU_NATIVE)
        if (_HAS_MCPU_NATIVE)
          list(APPEND _arch_flags_to_apply "-mcpu=native")
        else()
          list(APPEND _arch_flags_to_apply "-march=native")
        endif()
      else()
        if (MADNESS_TARGET_ARCH MATCHES "^armv")
          list(APPEND _arch_flags_to_apply "-march=${MADNESS_TARGET_ARCH}")
        else()
          list(APPEND _arch_flags_to_apply "-mcpu=${MADNESS_TARGET_ARCH}")
        endif()
      endif()
    endif()

    if (MADNESS_TUNE_ARCH)
      list(APPEND _arch_flags_to_apply "-mtune=${MADNESS_TUNE_ARCH}")
    endif()
  endif()

  foreach(_flag IN LISTS _arch_flags_to_apply)
    string(REGEX REPLACE "[^a-zA-Z0-9]" "_" _flag_var "MADNESS_SUPPORTS_${_flag}")
    check_cxx_compiler_flag("${_flag}" ${_flag_var})
    if (${_flag_var})
      add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:${_flag}>")
      message(STATUS "MADNESS: Enabled target architecture flag: ${_flag}")
    else()
      message(WARNING "MADNESS: Compiler does not support target architecture flag '${_flag}'")
    endif()
  endforeach()

elseif (NOT _madness_arch_from_user_flags)
  message(STATUS "MADNESS: Target architecture flags disabled (MADNESS_TARGET_ARCH=${MADNESS_TARGET_ARCH})")
endif()

# Floating-point relaxation mode. Opt-in: relaxed modes change results between builds (the reduction order depends
# on the vectorization the compiler chooses), which moves numerical references and regression baselines.
if (MADNESS_RELAXED_MATH STREQUAL "relaxed")
  set(_math_flags_to_check
      "-fassociative-math"
      "-fno-signed-zeros"
      "-fno-trapping-math"
      "-ffp-contract=fast"
  )
elseif (MADNESS_RELAXED_MATH STREQUAL "fast")
  set(_math_flags_to_check "-ffast-math")
elseif (MADNESS_RELAXED_MATH MATCHES "^(strict|STRICT|none|NONE|OFF|off|False|FALSE|0)$")
  set(_math_flags_to_check "")
else()
  message(WARNING "MADNESS: Unknown MADNESS_RELAXED_MATH mode '${MADNESS_RELAXED_MATH}', ignoring")
  set(_math_flags_to_check "")
endif()

foreach(_flag IN LISTS _math_flags_to_check)
  string(REGEX REPLACE "[^a-zA-Z0-9]" "_" _flag_var "MADNESS_SUPPORTS_${_flag}")
  check_cxx_compiler_flag("${_flag}" ${_flag_var})
  if (${_flag_var})
    add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:${_flag}>")
    message(STATUS "MADNESS: Enabled floating-point flag: ${_flag}")
  else()
    message(WARNING "MADNESS: Compiler does not support floating-point flag '${_flag}'")
  endif()
endforeach()

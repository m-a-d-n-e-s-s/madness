include_guard(GLOBAL)

# MADNESS internal warning-policy target.
#
# Owns an INTERFACE library `madness_internal_warnings` that carries the
# warning flags applied to MADNESS's own translation units. The target is
# linked PRIVATE-ly to MAD*-obj object libraries (by add_mad_library) and
# to in-tree executables (by add_mad_executable). PRIVATE scope is load-
# bearing: it keeps these flags out of INTERFACE_COMPILE_OPTIONS on the
# installed/exported MAD* targets, so downstream consumers (TiledArray,
# MPQC, …) do not inherit -Werror through find_package(madness).
#
# It is also NOT installed/exported, which keeps it out of the package.

add_library(madness_internal_warnings INTERFACE)

if (MADNESS_WERROR)
  if (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|IntelLLVM)$")
    target_compile_options(madness_internal_warnings INTERFACE -Werror)
  else()
    message(WARNING "MADNESS_WERROR=ON but compiler '${CMAKE_CXX_COMPILER_ID}' is not in the supported set; ignoring.")
  endif()
endif()

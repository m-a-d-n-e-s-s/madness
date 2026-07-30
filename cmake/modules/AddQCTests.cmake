# Register a quantum-chemistry example deck as a ctest case (label `qctest`).
#
# A qctest case is a directory holding an input deck, a one-liner run script, a
# declarative list of result keys to check, and a reference run:
#
#   <case>/input                      the deck, verbatim (what a user copies)
#   <case>/run.sh                     one-liner; honours ${MADQC}, ${MPIEXEC}, ${NP}
#   <case>/check.json                 which calc_info.json keys matter, and to what tolerance
#   <case>/reference/*.calc_info.json compared numerically by bin/run_qctest.py
#   <case>/reference/*.out            checked in for humans/agents; NOT compared
#
# Usage
#   add_qctest(scf_h2o_hf madqc "qctest;short")
#
# _labels must contain `qctest` plus a cost tier (short|medium|long|verylong).
# short/medium cases are picked up by check-short-madness for free, since that
# target is `ctest -L "short|medium"`.
#
# The case runs in a per-case scratch directory under the build tree, so the
# source directory stays read-only and cases cannot collide with each other.
# Every case may self-skip with exit code 77 (missing MPI, too few threads, ...).

macro(add_qctest _case _binary _labels)

  if (NOT TARGET qctests-madness)
    add_custom_target_subproject(madness qctests)
  endif()

  # build the driver binary on demand, mirroring add_scripted_tests
  if (NOT TEST madness/test/qc/build)
    add_test(NAME madness/test/qc/build
        COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${_binary})
  endif()
  set_property(TEST madness/test/qc/build PROPERTY LABELS "${_labels}" APPEND)
  add_dependencies(qctests-madness ${_binary})

  # scratch directory: the case is copied here and run, leaving the source tree clean
  set(_qctest_workdir "${CMAKE_CURRENT_BINARY_DIR}/${_case}")
  file(MAKE_DIRECTORY "${_qctest_workdir}")

  add_test(NAME madness/test/qc/${_case}/run
      COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/bin/run_qctest.py"
              --case "${CMAKE_CURRENT_SOURCE_DIR}/${_case}"
      WORKING_DIRECTORY "${_qctest_workdir}")

  set_tests_properties(madness/test/qc/${_case}/run PROPERTIES
      LABELS "${_labels}"
      DEPENDS madness/test/qc/build
      ENVIRONMENT "MADQC=$<TARGET_FILE:${_binary}>"
      SKIP_RETURN_CODE 77)

  # Cost tier -> timeout. Tiers are budgets on measured wall time (short < 10 s,
  # medium < 30 s, long < 2 min, verylong beyond); the timeout is a hang detector,
  # so it sits roughly an order of magnitude above the tier ceiling to absorb slow
  # or loaded machines without turning a slow run into a spurious failure.
  # Anything not tagged gets the ctest default.
  if ("${_labels}" MATCHES "verylong")
    set_tests_properties(madness/test/qc/${_case}/run PROPERTIES TIMEOUT 7200)
  elseif ("${_labels}" MATCHES "long")
    set_tests_properties(madness/test/qc/${_case}/run PROPERTIES TIMEOUT 1200)
  elseif ("${_labels}" MATCHES "medium")
    set_tests_properties(madness/test/qc/${_case}/run PROPERTIES TIMEOUT 300)
  elseif ("${_labels}" MATCHES "short")
    set_tests_properties(madness/test/qc/${_case}/run PROPERTIES TIMEOUT 120)
  endif()

endmacro()

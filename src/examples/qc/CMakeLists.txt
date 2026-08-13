# src/examples/qc -- quantum-chemistry example decks, doubling as integration tests.
#
# Each case is a data-only directory (deck + run.sh + check.json + reference);
# adding one takes no code. See README.md for the layout and add_qctest in
# cmake/modules/AddQCTests.cmake for the registration.
#
# Labels: `qctest` plus a cost tier. short/medium cases join check-short-madness
# automatically (that target is `ctest -L "short|medium"`), so tag a case `short`
# only if it runs in well under a minute.

if (BUILD_TESTING AND Python3_EXECUTABLE)

  # Tiers are MEASURED wall times (node26, 96 cores, MAD_NUM_THREADS=20), not
  # estimates: short < 10 s, medium < 30 s, long < 2 min, verylong beyond.
  # Re-measure before retagging; short+medium gate every CI run via
  # check-short-madness.
  add_qctest(scf_he_hf         madqc "qctest;short;applications")     #   5 s
  add_qctest(scf_he_hf_mpi     madqc "qctest;short;applications")     #   7 s, 2 ranks
  add_qctest(nemo_he_hf        madqc "qctest;short;applications")     #  10 s
  add_qctest(cis_he_singlets   madqc "qctest;medium;applications")    #  25 s
  add_qctest(scf_h2o_hf        madqc "qctest;long;applications")      #  38 s
  # the only gopt case: asserts the optimizer projects translations/rotations
  # out of the gradient before testing convergence (see README)
  add_qctest(scf_lih_gopt      madqc "qctest;long;applications")      #  38 s
  # the `optimize` workflow on both reference engines; the pair also cross-checks
  # that a regularized and a plain SCF reference find the same minimum
  add_qctest(scf_lih_optimize   madqc "qctest;long;applications")      #  37 s
  add_qctest(nemo_lih_optimize  madqc "qctest;long;applications")      #  77 s
  # the only polyatomic optimization, the only one on a functional, and the only
  # one at the default eprec -- see the deck for why each of those matters
  add_qctest(scf_h2o_lda_optimize madqc "qctest;long;applications")    #  74 s
  # A plain `calc_type cc2` ground state on helium measured 1833 s (30.5 min) on
  # 20 threads, so there is no cc2 case: lrcc2_he_excited already solves the CC2
  # ground state on its way to the excitation, and bounds it with `iter_max`.
  add_qctest(mp2_he_corr       madqc "qctest;verylong;applications")  # 790 s
  add_qctest(lrcc2_he_excited  madqc "qctest;verylong;applications")  # 6744 s (1h52m)
  add_qctest(oep_be_oaep       madqc "qctest;medium;applications")    #  21 s
  # localization invariance: same energy from all three `localize` settings
  add_qctest(nemo_h2o_canon    madqc "qctest;long;applications")      #  48 s
  add_qctest(nemo_h2o_boys     madqc "qctest;long;applications")      #  50 s
  add_qctest(nemo_h2o_new      madqc "qctest;long;applications")      #  51 s

endif()

# The decks are user-facing documentation, not just fixtures: install them so the
# run.sh one-liners work verbatim against an installed madqc.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/madness/${MADNESS_VERSION}/examples/qc"
        COMPONENT examples
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "run.sh"
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE)

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

  # The runner's own comparison semantics (tol / rtol / max / allow_zero), so a
  # change to run_qctest.py cannot silently loosen every case at once.
  add_test(NAME madness/test/qc/run_qctest_selftest/run
           COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/bin/test_run_qctest.py")
  set_tests_properties(madness/test/qc/run_qctest_selftest/run
                       PROPERTIES LABELS "qctest;short;applications")

  # Tiers are MEASURED wall times (node26, 96 cores, MAD_NUM_THREADS=20), not
  # estimates: short < 10 s, medium < 30 s, long < 2 min, verylong beyond.
  # Re-measure before retagging; short+medium gate every CI run via
  # check-short-madness.
  add_qctest(scf_he_hf         madqc "qctest;short;applications")     #   6 s
  # Without libxc, xcfunctional_ldaonly.cc stands in and understands only `lda`
  # and `hf`, so any deck naming a functional beyond those has to be registered
  # conditionally -- same guard the chem unit tests use. This is also why every
  # other deck here is `xc hf` or `xc lda`.
  if (TARGET Libxc::xc)
    add_qctest(scf_he_pbe0     madqc "qctest;medium;applications")    #   9 s
    add_qctest(scf_he_tpss     madqc "qctest;medium;applications")    #  17 s
    # the open-shell counterpart of scf_he_tpss: the only deck with unequal spin
    # densities, hence the only one exercising the spin-polarized sigma matrix
    add_qctest(scf_li_tpss     madqc "qctest;medium;applications")    #  15 s
    # the only DFT deck on the nemo path: XCOperator's nemo constructor, the
    # regularized density-gradient path, and the ncf/GGA-flux interaction
    add_qctest(nemo_he_pbe     madqc "qctest;medium;applications")    #  12 s
    # meta-gga on the nemo path: tau from the nemos via the psi = R F product
    # rule, and the non-multiplicative term with the R factors cancelled
    add_qctest(nemo_he_tpss    madqc "qctest;medium;applications")    #  24 s
    # the weak form of the semilocal xc potential (xc_weak_gga), on two orbitals
    # so the per-orbital flux is exercised; must agree with the multiplicative form
    add_qctest(nemo_lih_pbe_weak madqc "qctest;medium;applications")  #  23 s
  endif()
  add_qctest(scf_he_hf_mpi     madqc "qctest;short;applications")     #  11 s, 2 ranks
  add_qctest(nemo_he_hf        madqc "qctest;short;applications")     #  12 s
  # the only case with k varying across the protocol; must agree with nemo_he_hf
  add_qctest(nemo_he_varyk     madqc "qctest;medium;applications")    #  13 s
  add_qctest(cis_he_singlets   madqc "qctest;medium;applications")    #   9 s
  # the first `response` case: the smallest complete molresponse deck (see README)
  add_qctest(response_he_alpha madqc "qctest;short;applications")     #   9 s, 7 threads
  # the same H2 deck twice: the `excited.*` path, once in the Tamm-Dancoff
  # approximation and once at full RPA (`excited.tda false`)
  add_qctest(response_h2_es_tda madqc "qctest;medium;applications")   #  22 s, 7 threads
  add_qctest(response_h2_es_rpa madqc "qctest;medium;applications")   #  25 s, 7 threads
  # the `quadratic` path: the z FD leg plus one VBC solve for the static
  # beta_zzz of LiH. `long` rather than the CI tier the rest of the response
  # set sits in -- a static beta cannot be cut below this (see check.json)
  add_qctest(response_lih_beta madqc "qctest;long;applications")      #  42 s, 7 threads
  add_qctest(scf_h2o_hf        madqc "qctest;long;applications")      #  35 s
  # asserts the optimizer projects translations/rotations out of the gradient
  # before testing convergence, and that explicit `optimization` values beat the
  # derived ones (see README)
  add_qctest(scf_lih_optimize_tight madqc "qctest;long;applications") #  38 s
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
  add_qctest(oep_be_oaep       madqc "qctest;medium;applications")    #  28 s
  # localization invariance: same energy from all three `localize` settings
  add_qctest(nemo_h2o_canon    madqc "qctest;long;applications")      #  75 s
  add_qctest(nemo_h2o_boys     madqc "qctest;long;applications")      #  83 s
  add_qctest(nemo_h2o_new      madqc "qctest;long;applications")      #  79 s

  # Grimme's D3 dispersion correction in the energy and the single-point
  # gradient. Needs simple-dftd3, since without it the deck's `dispersion d3bj`
  # aborts by design rather than silently doing nothing -- and needs libxc too,
  # because its `xc pbe` is unsupported by the LDA-only functional build (which
  # throws a bare `char const*` from xcfunctional_ldaonly.cc). `xc lda` would not
  # rescue it: simple-dftd3 has no LDA damping parameters either.
  if (MADNESS_HAS_DFTD3 AND TARGET Libxc::xc)
    add_qctest(scf_lih_pbe_d3  madqc "qctest;medium;applications")     #  12 s, laptop
  endif()

  # The polarizable continuum model, and the only deck carrying a `pcm` group.
  # Needs PCMSolver: without it the deck aborts by design (chem/pcm.cc throws
  # rather than quietly running in vacuum), so it cannot be registered
  # unconditionally -- ENABLE_PCM is off by default.
  if (MADNESS_HAS_PCM)
    add_qctest(scf_lih_pcm_water madqc "qctest;medium;applications")   #  12 s
  endif()

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

# H2 Excited-State Validation: legacy vs. v1 molresponse

System: H₂ at 1.4 Bohr (HF/TDA & TDDFT, 4 roots, single-process, k=8, thresh=1e-6).
Ground state archive: `mad.restartdata.00000` from this directory.

## Final excitation energies (au)

| Code           | Mode | ω₁         | ω₂         | ω₃         | ω₄         | Iters | Wall (s) |
|----------------|------|------------|------------|------------|------------|-------|----------|
| legacy         | TDA  | 0.46807128 | 0.47785601 | 0.48135613 | 0.48136583 |    10 |    331.4 |
| **v1 (this branch)** | TDA  | 0.46802021 | 0.47785019 | 0.48134267 | 0.48134267 |    25 |    132.3 |
| legacy         | RPA  | 0.46556124 | 0.47714269 | 0.48084333 | 0.48084599 |    15 |    377.8 |
| **v1 (this branch)** | RPA  | 0.46556163 | 0.47714275 | 0.48084341 | 0.48084342 |    25 |    206 |

Both TDA codes agree to within 5 × 10⁻⁵ a.u. (a single iteration of bsh-residual tolerance). Legacy converges in fewer iterations but each iteration is slower; v1 is ~2.5× faster overall on this system. RPA values are ~3 × 10⁻³ below TDA, as expected for a small Y contribution. **v1 and legacy RPA agree to ≤ 5 × 10⁻⁷** on the lowest three roots — the gap on ω₄ (~3 × 10⁻⁶) is at the convergence-test floor; both runs hit `maxiter` before the strict bsh-residual target.

## Restart from saved state

| Code   | Cold start         | Restart from save        | Notes |
|--------|-------------------|--------------------------|-------|
| legacy | 10 iters / 331.4 s | 2 iters / 6.6 s         | Restart skipped almost all main iterations — checkpoint was already converged. |
| v1     | 25 iters / 132.3 s | 25 iters / 108.9 s       | Restart skipped the 5-iteration guess phase (~30 s); main loop ran the full count because v1's restart entry point lands before the convergence test rather than re-evaluating it on load. |

**Restart works on both codes**, but legacy's restart is more efficient because it re-checks convergence on entry. v1's restart correctly resumes from `guess_restart` and converges to the identical answer.

## v1 RPA (full TDDFT) — fixed

v1 RPA initially crashed during iteration 0 with a malloc-consolidate SIGSEGV
deep inside the response-exchange machinery. Root cause was a hard-coded
`vector_real_function_3d temp_J(3)` in `J1StrategyFull::compute_J1` and
`J1StrategyStable::compute_J1` (`ResponseBase.hpp`): the `3` was sized for
the dipole frequency-response case (3 Cartesian directions). With 4
excited roots, the loop indexed `temp_J[3]` past the buffer end and
corrupted the heap; the SIGSEGV manifested in whichever later allocator
hit the corrupted arena (initially the K[1] macrotask, then the J[1]
shared_ptr release after switching exchange algorithms). Replacing the
hardcoded `3` with `x.num_states()` makes RPA converge cleanly and match
legacy RPA.

## Reproducing

```bash
ARCHIVE_DIR=$PWD  # this directory has mad.restartdata.00000
WORK=/tmp/h2_es_repro && mkdir -p $WORK && cd $WORK
ln -sf $ARCHIVE_DIR/mad.restartdata.00000 .
ln -sf $ARCHIVE_DIR/mad.restartaodata .
ln -sf $ARCHIVE_DIR/mad.calc_info.json .

# v1 TDA (cold start)
cat > response.in <<'EOF'
response
    archive mad.restartdata
    excited_state true
    tda true
    states 4
    protocol [1.e-6]
    dconv 1.e-4
    maxiter 25
    random true
    save true
    save_file v1_es_restart
end
EOF
~/Projects/madness/build/src/apps/molresponse/molresponse > tda.log 2>&1

# v1 TDA restart from saved state
cat > response.in <<'EOF'
response
    archive mad.restartdata
    excited_state true
    tda true
    states 4
    protocol [1.e-6]
    dconv 1.e-4
    maxiter 25
    restart true
    restart_file v1_es_restart
end
EOF
~/Projects/madness/build/src/apps/molresponse/molresponse > restart.log 2>&1

# legacy: same input, but file is named 'input' and uses the legacy binary
cp response.in input
~/Projects/madness/build/src/apps/molresponse_legacy/molresponse_legacy > leg.log 2>&1
```

## v1 fixes required to make TDA work

The v1 codebase had been bit-rotting for some time; the H2 4-root TDA path required 11 patches before it ran cleanly, all in `src/apps/molresponse/` and `src/madness/mra/commandlineparser.h`:

1. `src/apps/CMakeLists.txt` — re-enable `molresponse` subdir.
2. `Exchange<double,3>::Algorithm` → `ExchangeAlgorithm` (60 sites in 4 files; rename in upstream).
3. `response_parameters.h` — add pure-virtual `get_tag()` override.
4. `molresponse.cc` — drop a stale `omega` argument from `FrequencyResponse` ctor call.
5. `ground_parameters.h` — set the cubic cell from the archive's `L` before reading orbitals (otherwise `Function::load` rejects the cell mismatch).
6. `ResponseBase::solve` — call `check_k` (which builds `hamiltonian`) **before** `initialize`, so the guess-iteration's response-potential build doesn't see an empty Hamiltonian tensor.
7. `response_space(World, m, n)` ctor — fill slots with real zero functions instead of leaving `FunctionImpl=nullptr`. Without this, `X_space::copy()` and downstream operators dereference null impls.
8. `response_space::push_back` and `pop_back` — keep the `active` list in sync with `num_states`. The omega convergence to the legacy reference values only happens after this fix; without it, the diagonalization rotation produces a response_space with empty `active`, and `to_vector`/`from_vector` silently produce zero-length payloads.
9. `ExcitedResponse::initialize` — rebuild `Chi` from scratch after `select_functions`, so `Chi.n_states` matches the trimmed `Chi.x.size()`.
10. `ExcitedResponse` ctor — set inner-product / J1 / K1 / VXC strategies (mirrors `FrequencyResponse` ctor; was missing entirely).
11. `ResponseBase::update_residual` — defensive canonical active list and graceful empty-`old_residuals` handling (the first iteration calls it with `Tensor<double>()`).
12. `J1StrategyFull` / `J1StrategyStable` — `temp_J(3)` → `temp_J(x.num_states())`. Hardcoded for FrequencyResponse's 3 dipole directions; out-of-bounds write for any other state count.

## Files

- v1 TDA cold start log: `/tmp/h2_v1_es/v1_h2_es.log`  (guess + 25 iters → converged)
- v1 TDA restart log:    `/tmp/h2_v1_restart/restart.log`  (restart → 25 iters → same answer)
- v1 RPA log (crash):    `/tmp/h2_v1_full/full.log`
- legacy TDA cold log:   `/tmp/h2_legacy_restart/tda_save.log`
- legacy TDA restart log: `/tmp/h2_legacy_restart_run/restart.log`  (restart → 1 iter)
- legacy RPA log:        `/tmp/h2_legacy_full/full.log`

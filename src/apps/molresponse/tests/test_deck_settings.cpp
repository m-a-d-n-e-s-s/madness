// =========================================================================
// apply_deck_es_knobs — the deck's excited.* iteration keys reach the executor.
//
// review/findings C2: excited.maxiter, excited.guess_max_iter and
// excited.maxsub were parsed, printed and then dropped; the ES solves ran on
// response.maxiter and the ExecutorSettings defaults whatever the deck said.
// A key the deck does not set keeps the ExecutorSettings default (the
// ResponseParameters defaults are not the solver's).
//
// Pure parameter mapping; no World needed.
// =========================================================================

#include "../deck_settings.hpp"

#include <cstdio>

namespace {

int failed = 0;

void expect(bool cond, const char *label) {
  std::printf("  [%s]  %s\n", cond ? "PASS" : "FAIL", label);
  if (!cond) ++failed;
}

} // namespace

int main() {
  using molresponse_v3::ExecutorSettings;
  using molresponse_v3::apply_deck_es_knobs;

  std::printf("=== deck sets all three keys ===\n");
  {
    ResponseParameters rp;
    rp.set_user_defined_value<size_t>("excited.maxiter", 17);
    rp.set_user_defined_value<size_t>("excited.guess_max_iter", 3);
    rp.set_user_defined_value<size_t>("excited.maxsub", 4);
    ExecutorSettings s;
    apply_deck_es_knobs(rp, s);
    expect(s.es_max_iters == 17, "excited.maxiter -> es_max_iters");
    expect(s.es_tda_warmup_iters == 3, "excited.guess_max_iter -> es_tda_warmup_iters");
    expect(s.es_kain_maxsub == 4, "excited.maxsub -> es_kain_maxsub");
  }

  std::printf("=== deck sets none ===\n");
  {
    ResponseParameters rp;
    const ExecutorSettings def;
    ExecutorSettings s;
    apply_deck_es_knobs(rp, s);
    expect(s.es_max_iters == def.es_max_iters && s.es_max_iters == 0,
           "unset excited.maxiter inherits response.maxiter (es_max_iters 0)");
    expect(s.es_tda_warmup_iters == def.es_tda_warmup_iters,
           "unset excited.guess_max_iter keeps the executor default");
    expect(s.es_kain_maxsub == def.es_kain_maxsub,
           "unset excited.maxsub keeps the executor default");
  }

  std::printf("=== ES budget ===\n");
  {
    ExecutorSettings s;
    s.max_iters = 30;
    expect(s.es_iter_budget() == 30, "es_max_iters 0 -> max_iters");
    s.es_max_iters = 12;
    expect(s.es_iter_budget() == 12, "es_max_iters > 0 wins");
  }

  // review/findings C3: deck `kain` / `maxrotn` reached only FD; both ES
  // executors overwrote them before building the main ES policy.
  std::printf("=== deck kain / maxrotn reach ES ===\n");
  {
    ResponseParameters rp;
    rp.set_user_defined_value<bool>("kain", false);
    rp.set_user_defined_value<double>("maxrotn", 0.3);
    ExecutorSettings s;
    apply_deck_es_knobs(rp, s);
    expect(!s.es_kain, "deck kain false -> es_kain false");
    expect(s.es_maxrotn == 0.3, "deck maxrotn -> es_maxrotn");

    ResponseParameters unset;
    const ExecutorSettings def;
    ExecutorSettings u;
    apply_deck_es_knobs(unset, u);
    expect(u.es_kain == def.es_kain && u.es_kain,
           "unset kain keeps ES KAIN on (the deck default kain=false is not applied)");
    expect(u.es_maxrotn == def.es_maxrotn, "unset maxrotn keeps es_maxrotn");
  }

  std::printf("=== ES warmup / main policies ===\n");
  {
    using molresponse_v3::detail_exec::es_iteration_policies;
    ExecutorSettings s;
    s.policy.kain_min_residual = 0.02;
    s.es_kain = false;
    s.es_maxrotn = 0.3;
    s.es_kain_maxsub = 4;
    s.es_tda_warmup_iters = 3;
    s.es_main_kain_delay = 2;
    const auto p = es_iteration_policies(s);
    expect(!p.main.kain, "es_kain false -> main ES solve runs without KAIN");
    expect(p.warm.tda_warmup_iters == 3,
           "warmup is a KAIN-free window of es_tda_warmup_iters");
    expect(p.main.maxrotn == 0.3 && p.warm.maxrotn == 0.3, "es_maxrotn -> both");
    expect(p.main.kain_maxsub == 4 && p.warm.kain_maxsub == 4, "es_kain_maxsub -> both");
    expect(p.main.tda_warmup_iters == 2, "main KAIN delay = es_main_kain_delay");
    expect(p.main.kain_min_residual == 0.02, "kain.min_residual reaches the main solve");
    expect(p.main.lock_converged == s.es_lock_converged && !p.warm.lock_converged,
           "only the main solve locks converged roots");

    const auto d = es_iteration_policies(ExecutorSettings{});
    expect(d.main.kain && d.warm.kain, "defaults: KAIN on");
  }

  std::printf("\n%s (%d failure(s))\n", failed ? "FAILED" : "PASSED", failed);
  return failed ? 1 : 0;
}

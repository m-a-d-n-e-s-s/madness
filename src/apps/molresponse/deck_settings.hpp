#pragma once

// Deck (ResponseParameters) -> ExecutorSettings mapping for the excited.*
// iteration keys. Separate from madqc_adapter.hpp so it can be unit-tested
// without an SCF.

#include "calc/calc_executor.hpp"

#include <madness/chem/ResponseParameters.hpp>

namespace molresponse_v3 {

/// Only keys the deck SETS are applied: the ResponseParameters defaults
/// (guess_max_iter 5, maxiter 20) are not the executor's (warmup 10, ES budget
/// inherited from response.maxiter), and the executor defaults are the
/// validated ones.
inline void apply_deck_es_knobs(const ResponseParameters &rp,
                                ExecutorSettings &s) {
  if (rp.is_user_defined("excited.maxiter"))
    s.es_max_iters = static_cast<int>(rp.excited_maxiter());
  if (rp.is_user_defined("excited.guess_max_iter"))
    s.es_tda_warmup_iters = static_cast<int>(rp.excited_guess_max_iter());
  if (rp.is_user_defined("excited.maxsub"))
    s.es_kain_maxsub = static_cast<int>(rp.excited_maxsub());
}

} // namespace molresponse_v3

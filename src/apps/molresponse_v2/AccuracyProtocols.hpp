#ifndef ACCURACYPROTOCOLS_HPP
#define ACCURACYPROTOCOLS_HPP

#include "ResponseManager.hpp"
#include "ResponseState.hpp"
#inclulde <
#include <iostream>
#include <vector>

class AccuracyProtocols {

public:
  virtual void refine(ResponseState &state, ResponseManager &manager) = 0;
};

class IncrementalAccuracy : public AccuracyProtocols {
public:
  IncrementalAccuracy(std::vector<double> thresholds)
      : thresholds_(thresholds) {}

  void refine(ResponseState &state, ResponseManager &manager) override {
    for (auto threshold : thresholds_) {
      state.convergence_threshold = threshold;
      if (!state.is_converged) {
        manager.computeState(state);
      }
    }
  }

private:
  std::vector<double> thresholds_;
};

#endif // ACCURACYPROTOCOLS_HPP

#ifndef RESPONSE_STATE_HPP
#define RESPONSE_STATE_HPP
#include <iostream>

struct ResponseState {
  enum class Perturbation { X, Y, Z };

  Perturbation perturbation;

  double frequency;
  double convergence_threshold;
  bool is_converged;

  [[nodiscard]] std::string getIdentifer() const {

    return perturbationToString(perturbation) + "_" +
           std::to_string(frequency) + "_" +
           std::to_string(convergence_threshold);
  }

private:
  static std::string perturbationToString(Perturbation p) {
    switch (p) {
    case Perturbation::X:
      return "X";
    case Perturbation::Y:
      return "Y";
    case Perturbation::Z:
      return "Z";
    }
    return "Unknown";
  }
};

#endif // RESPONSE_STATE_HPP


// Tasks.hpp defines the available tasks in MADQC which are
// energy, optimization, property
//
//
// We will implement energy and optimization tasks using a strategy design
// pattern
//

#ifndef MADQC_TASKS_HPP
#define MADQC_TASKS_HPP

#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>

namespace madness {

struct TaskParameters : public QCCalculationParametersBase {
  TaskParameters(const TaskParameters& other) = default;

  TaskParameters(World& world, const commandlineparser& parser)
      : TaskParameters() {
    read_input_and_commandline_options(world, parser, "task");
  }
  TaskParameters() {
    initialize<std::string>("driver", "energy",
                            "energy, optimization, property,'oniom");
    initialize<std::string>("method", "dft", "cc2, tddft, mp2");
  }

  using QCCalculationParametersBase::read_input_and_commandline_options;

  [[nodiscard]] std::string get_driver() const {
    return get<std::string>("driver");
  }
  [[nodiscard]] std::string get_method() const {
    return get<std::string>("method");
  }
};

class EnergyTaskStrategy {
 public:
  virtual double compute_energy() = 0;
  virtual ~EnergyTaskStrategy() = default;
};

class DFTEnergyTask : public EnergyTaskStrategy {

 public:
  double compute_energy() override { 



    return 0.0; }
};

class MP2EnergyTask : public EnergyTaskStrategy {
  double compute_energy() override { return 0.0; }
};

class CC2EnergyTask : public EnergyTaskStrategy {
  double compute_energy() override { return 0.0; }
};

// context class
class EnergyCalculator {
 private:
  EnergyTaskStrategy* strategy;

 public:
  explicit EnergyCalculator(EnergyTaskStrategy* strategy) : strategy(strategy) {}

  double compute_energy() { return strategy->compute_energy(); }

  void set_strategy(EnergyTaskStrategy* strategy) { this->strategy = strategy; }
};

}  // namespace madness
#endif  // MADQC_TASKS_HPP

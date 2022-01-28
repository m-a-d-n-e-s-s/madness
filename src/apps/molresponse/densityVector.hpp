//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_DENSITYVECTOR_HPP
#define MADNESS_DENSITYVECTOR_HPP
#include "molresponse/global_functions.h"
#include "property.h"

#include "mra.h"
#include "ResponseBase.hpp"

class DensityVector {
 public:
  DensityVector() = default;
  virtual void compute(World& world) = 0;
  virtual ~DensityVector() = default;
};

class FrequencyVector;

class FrequencySolver : public ResponseBase {
 public:
  virtual ~FrequencySolver() = default;
  void compute(World&world,FrequencyVector& density);
 private:
  PropertyBase p;                   // for frequency calculations
  X_space PQ;
  Tensor<double> omega;
  response_space PropertyRHS(World& world, PropertyBase& p) const;
  void iterate(World& world);
  void compute_and_print_polarizability(World& world, std::string message);
  void frequency_to_json(json& j_mol_in,
                         size_t iter,
                         const Tensor<double>& res_X,
                         const Tensor<double>& res_Y,
                         const Tensor<double>& density_res,
                         const Tensor<double>& omega);
  void update_x_space_response(World& world,
                               X_space& Chi,
                               X_space& res,
                               XCOperator<double, 3>& xc,
                               vector<poperatorT>& bsh_x_ops,
                               vector<poperatorT>& bsh_y_ops,
                               QProjector<double, 3>& projector,
                               double& x_shifts,
                               double& omega_n,
                               NonLinearXsolver& kain_x_space,
                               vector<X_vector> Xvector,
                               vector<X_vector> Xresidual,
                               Tensor<double>& bsh_residualsX,
                               Tensor<double>& bsh_residualsY,
                               size_t iteration,
                               Tensor<double>& maxrotn);
  X_space bsh_update_response(World& world,
                              X_space& theta_X,
                              vector<poperatorT>& bsh_x_ops,
                              vector<poperatorT>& bsh_y_ops,
                              QProjector<double, 3>& projector,
                              double& x_shifts);
};


class FrequencyVector : public DensityVector {
 public:
  FrequencyVector(World& world, std::unique_ptr<FrequencySolver> calc) : computing(std::move(calc)) {}

 protected:
  std::unique_ptr<FrequencySolver> computing;
 private:
  X_space PQ;
  Tensor<double> omega;
};






#endif  // MADNESS_DENSITYVECTOR_HPP

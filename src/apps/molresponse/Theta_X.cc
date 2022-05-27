

#include <stdio.h>


#include <memory>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/xcfunctional.h"
#include "molresponse/TDDFT.h"

#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"

#include "molresponse/x_space.h"
using json = nlohmann::json;


X_space TDDFT::Compute_Theta_X(World& world, X_space& Chi, XCOperator<double, 3> xc, std::string calc_type) {
  bool compute_Y = calc_type.compare("full") == 0;
  X_space Theta_X = X_space(world, Chi.num_states(), Chi.num_orbitals());
  // compute
  X_space V0X = compute_V0X(world, Chi, xc, compute_Y);

  V0X.truncate();
  if (r_params.print_level() >= 20) {
    print("---------------Theta ----------------");
    print("<X|V0|X>");
    print(inner(Chi, V0X));
  }

  X_space E0X(world, Chi.num_states(), Chi.num_orbitals());
  if (r_params.localize().compare("canon") == 0) {
    E0X = Chi.copy();
    E0X.truncate();
    E0X.X = E0X.X * ham_no_diag;
    if (compute_Y) {
      E0X.Y = E0X.Y * ham_no_diag;
    }

    E0X.truncate();
  }

  if (r_params.print_level() >= 20) {
    print("<X|(E0-diag(E0)|X>");
    print(inner(Chi, E0X));
  }

  X_space gamma;
  // compute
  if (calc_type.compare("full") == 0) {
    gamma = compute_gamma_full(world, Chi, xc);
  } else if (calc_type.compare("static") == 0) {
    gamma = compute_gamma_static(world, Chi, xc);
  } else {
    gamma = compute_gamma_tda(world, Chi, xc);
  }

  Theta_X = (V0X - E0X) + gamma;
  Theta_X.truncate();

  if (r_params.print_level() >= 20) {
    print("<X|Theta|X>");
    print(inner(Chi, Theta_X));
  }

  return Theta_X;
}
void TDDFT::output_json() const {
  std::ofstream ofs("response.json");
  ofs << j_molresponse;
}
void TDDFT::excited_to_json(json& j_mol_in,
                            size_t iter,
                            const Tensor<double>& res_X,
                            const Tensor<double>& res_Y,
                            const Tensor<double>& density_res,
                            const Tensor<double>& omega) {
  json j = {};
  j["iter"] = iter;
  j["res_X"] = tensor_to_json(res_X);
  j["res_Y"] = tensor_to_json(res_Y);
  j["density_residuals"] = tensor_to_json(density_res);
  j["omega"] = tensor_to_json(omega);
  auto index = j_mol_in["protocol_data"].size() - 1;
  j_mol_in["protocol_data"][index]["iter_data"].push_back(j);
}
void TDDFT::frequency_to_json(json& j_mol_in,
                              size_t iter,
                              const Tensor<double>& res_X,
                              const Tensor<double>& res_Y,
                              const Tensor<double>& density_res,
                              const Tensor<double>& omega) {
  json j = {};
  j["iter"] = iter;
  j["res_X"] = tensor_to_json(res_X);
  j["res_Y"] = tensor_to_json(res_Y);
  j["density_residuals"] = tensor_to_json(density_res);
  j["polar"] = tensor_to_json(omega);
  auto index = j_mol_in["protocol_data"].size() - 1;
  j_mol_in["protocol_data"][index]["iter_data"].push_back(j);
}

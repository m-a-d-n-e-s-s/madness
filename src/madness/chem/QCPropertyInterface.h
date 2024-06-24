#ifndef MADNESS_QCPROPERTYINTERFACE_H
#define MADNESS_QCPROPERTYINTERFACE_H

#include <madness/mra/mra.h>
#include<madness/chem/molecule.h>
//class NuclearCorrelationFactor;

namespace madness {

/// class implementing properties of QC models
class QCPropertyInterface {

  virtual std::string name() const = 0;

  virtual bool selftest() = 0;

  virtual Tensor<double> nuclear_derivative(
      const real_function_3d &density, const Molecule &molecule,
      const std::shared_ptr<NuclearCorrelationFactor> ncf = 0) const {
    print("nuclear_derivative not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
    return Tensor<double>();
  }

  virtual real_function_3d density() const {
    print("density not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
  }

  virtual real_function_3d no_cusp_density() const {
    print("no_cusp_density not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
  }

  virtual real_function_3d spindensity(const int spin) const {
    MADNESS_ASSERT(spin == 0 or spin == 1); // alpha or beta
    print("spindensity not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
  }

  virtual real_function_3d no_cusp_spindensity(const int spin) const {
    MADNESS_ASSERT(spin == 0 or spin == 1); // alpha or beta
    print("no_cusp_spindensity not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
  }

  virtual std::vector<double> multipole_moment(
      const real_function_3d &density, const int l, const Molecule &molecule,
      const std::shared_ptr<NuclearCorrelationFactor> ncf = 0) const {
    print("dipole moment not implemented in ", name());
    MADNESS_EXCEPTION("feature not implemented", 1);
  }
};

}  // namespace madness

#endif //MADNESS_QCPROPERTYINTERFACE_H

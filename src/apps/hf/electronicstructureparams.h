#ifndef ELECTRONICSTRUCTUREPARAMS_H_
#define ELECTRONICSTRUCTUREPARAMS_H_

#include <mra/mra.h>
#include <misc/ran.h>
#include <misc/misc.h>
#include "mentity.h"
#include "molecularbasis.h"

namespace madness {

struct ElectronicStructureParams
{
  // Size of the cubic box (this needs to change)
  double L;
  // Amount of electronic charge
  int ncharge;
  // 1 - LDA; 2 - Hartree-Fock
  int functional;
  // Low value in the BSH / Coulomb fit
  double lo;
  // Spin-polarized
  bool spinpol;
  // Periodic system
  bool periodic;
  // Maximum number of interations
  int maxits;
  // Thresh
  double thresh;
  // Order of wavelets
  int waveorder;
  // Number of empty states
  int nempty;
  // Smearing parameter
  double smear;
  // Total number of bands
  int nbands;
  // Size of k-mesh (hardcoded for 3-d)
  int ngridk0, ngridk1, ngridk2;
  // Maximum occupation
  double maxocc;

  ElectronicStructureParams()
  {
    L = 10.0;
    ncharge = 1;
    functional = 1;
    lo = 1e-4;
    smear = 0.001;
    spinpol = false;
    periodic = false;
    maxits = 100;
    thresh = 1e-6;
    waveorder = 8;
    nempty = 2;
    nbands = ncharge + nempty;
    ngridk0 = 1; ngridk1 = 1; ngridk2 = 1;
    maxocc = 2.0;
  }

  template <typename Archive>
  void serialize(Archive& ar) {
      ar & L & ncharge & functional & lo & smear & spinpol & periodic &
      maxits & thresh & waveorder & nempty & nbands &
      ngridk0 & ngridk1 & ngridk2 & maxocc;
  }

  void read_file(const std::string& filename)
  {
    std::ifstream f(filename.c_str());
    position_stream(f, "dft");
    string s;
    while (f >> s)
    {
      if (s == "end")
      {
        break;
      }
      else if (s == "L")
      {
        f >> L;
      }
      else if (s == "functional")
      {
        f >> functional;
      }
      else if (s == "lo")
      {
        f >> lo;
      }
      else if (s == "smear")
      {
        f >> smear;
      }
      else if (s == "spinpol")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          spinpol = true;
        }
        else if (tempstr == "false")
        {
          spinpol = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- spinpol", 0);
        }
      }
      else if (s == "periodic")
      {
        std::string tempstr;
        f >> tempstr;
        if (tempstr == "true")
        {
          periodic = true;
        }
        else if (tempstr == "false")
        {
          periodic = false;
        }
        else
        {
          MADNESS_EXCEPTION("input error -- periodic", 0);
        }
      }
      else if (s == "maxits")
      {
        f >> maxits;
      }
      else if (s == "thresh")
      {
        f >> thresh;
      }
      else if (s == "waveorder")
      {
        f >> waveorder;
      }
      else if (s == "nempty")
      {
        f >> nempty;
      }
      else if (s == "ngridk")
      {
        f >> ngridk0; f >> ngridk1; f >> ngridk2;
      }
      else
      {
        std::cout << "moldft: unrecognized input keyword " << s << std::endl;
        MADNESS_EXCEPTION("input error", 0);
      }
    }
    // compute total number of bands
    nbands = ncharge + nempty;
    // maximum occupation
    maxocc = (spinpol) ? 1.0 : 2.0;
  }

  void set_molecular_info(const MolecularEntity& mentity, const AtomicBasisSet& aobasis) {
      lo = mentity.smallest_length_scale();
  }
};

}
#endif /* ELECTRONICSTRUCTUREPARAMS_H_ */

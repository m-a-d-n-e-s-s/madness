
#ifndef MADNESS_MADCHEM_H
#define MADNESS_MADCHEM_H

// numerical library
#include<madness.h>

// convenience
#include<madness/world/timing_utilities.h>
#include<madness/misc/info.h>

// the KAIN nonlinear solver
#include<madness/mra/nonlinsol.h>

// quantum chemistry headers
#include<madness/chem/AC.h>
#include<madness/chem/BSHApply.h>
#include<madness/chem/CC2.h>
#include<madness/chem/CCPotentials.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/CalculationParameters.h>
#include<madness/chem/ParameterManager.hpp>
#include<madness/chem/ESInterface.h>
#include<madness/chem/GuessFactory.h>
#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/NWChem.h>
#include<madness/chem/PNO.h>
#include<madness/chem/PNOF12Potentials.h>
#include<madness/chem/PNOGuessFunctions.h>
#include<madness/chem/PNOParameters.h>
#include<madness/chem/PNOStructures.h>
#include<madness/chem/PNOTensors.h>
#include"madness/mra/QCCalculationParametersBase.h"
#include<madness/chem/QCPropertyInterface.h>
#include<madness/chem/SCF.h>
#include<madness/chem/SCFOperators.h>
#include<madness/chem/SCFProtocol.h>
#include<madness/chem/TDHF.h>
#include<madness/chem/atomutil.h>
#include<madness/chem/basis.h>
#include<madness/chem/ccpairfunction.h>
#include"madness/mra/commandlineparser.h"
#include<madness/chem/corepotential.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/diamagneticpotentialfactor.h>
#include<madness/chem/distpm.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/exchangeoperator.h>
#include<madness/chem/gaussian.h>
#include<madness/chem/gth_pseudopotential.h>
#include<madness/chem/localizer.h>
#include<madness/chem/masks_and_boxes.h>
#include<madness/chem/molecular_functors.h>
#include<madness/chem/molecular_optimizer.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/molecule.h>
#include<madness/chem/molopt.h>
#include<madness/chem/mp2.h>
#include<madness/chem/nemo.h>
#include<madness/chem/oep.h>
#include<madness/chem/pcm.h>
#include<madness/chem/pointgroupoperator.h>
#include<madness/chem/pointgroupsymmetry.h>
#include<madness/chem/polynomial.h>
#include<madness/chem/potentialmanager.h>
#include<madness/chem/projector.h>
#include<madness/chem/vibanal.h>
#include<madness/chem/write_test_input.h>
#include<madness/chem/xcfunctional.h>
#include<madness/chem/zcis.h>
#include<madness/chem/znemo.h>

#endif //MADNESS_MADCHEM_H

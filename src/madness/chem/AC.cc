/*
 * AC.cc
 *
 *  Created on: Dec 15, 2016
 *      Author: msahre
 */
#include "AC.h"

namespace madness{

double slater_radius(int atomic_number){
	const double ang2b = 0.52917721092;
	if(atomic_number == 1){return 0.35/ang2b;}
	else if (atomic_number == 2) {return 0.35/ang2b;}
	else if (atomic_number == 3) {return 1.45/ang2b;}
	else if (atomic_number == 4) {return 1.05/ang2b;}
	else if (atomic_number == 5) {return 0.85/ang2b;}
	else if (atomic_number == 6) {return 0.70/ang2b;}
	else if (atomic_number == 7) {return 0.65/ang2b;}
	else if (atomic_number == 8) {return 0.60/ang2b;}
	else if (atomic_number == 9) {return 0.50/ang2b;}
	else if (atomic_number == 10) {return 1.60/ang2b;}
	else if (atomic_number == 11) {return 1.80/ang2b;}
	else if (atomic_number == 12) {return 1.50/ang2b;}
	else if (atomic_number == 13) {return 1.25/ang2b;}
	else if (atomic_number == 14) {return 1.10/ang2b;}
	else if (atomic_number == 15) {return 1.00/ang2b;}
	else if (atomic_number == 16) {return 1.00/ang2b;}
	else if (atomic_number == 17) {return 1.00/ang2b;}
	else if (atomic_number == 18) {return 0.71/ang2b;}
	else if (atomic_number == 19) {return 2.20/ang2b;}
	else if (atomic_number == 20) {return 1.80/ang2b;}
	else{
		MADNESS_EXCEPTION("Slater radius for element does not exist!",1);
	}
}

std::vector<atom_information<3> >make_atom_vec(const Molecule& molecule, double R1_, double R2_) {
	std::vector< atom_information<3> > atom_vec;
	for(const auto& x:molecule.get_atoms()){
		atom_information<3> tmp;
		tmp.coord = x.get_coords();
		tmp.charge = x.get_atomic_number();
		tmp.R1 = R1_*slater_radius(x.get_atomic_number());
		tmp.R2 = R2_*slater_radius(x.get_atomic_number());

		atom_vec.push_back(tmp);
	}
	return atom_vec;
}

} // end namespace madness

/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

#include<madness/chem/correlationfactor.h>

namespace madness{

	/// create and return a new nuclear correlation factor

	/// note there is also an Ad-hoc nuclear correlation factor, which can only be created directly
	/// @param[in]	world	the world
	/// @param[in]	calc	the calculation as read from the input file
	/// @return 	a nuclear correlation factor
	std::shared_ptr<NuclearCorrelationFactor>
	create_nuclear_correlation_factor(World& world,
			const Molecule& molecule,
			const std::shared_ptr<PotentialManager> potentialmanager,
			const std::string inputline) {

		std::stringstream ss(lowercase(inputline));
		std::string corrfac, factor;
		ss >> corrfac >> factor;

		// read the length scale factor if there is one
		double a=0.0;
		if (factor.size()>0) {
			std::stringstream fss(factor);
			if (not (fss >> a)) {
				if (world.rank()==0) print("could not read the length scale parameter a: ",a);
				MADNESS_EXCEPTION("input error in the nuclear correlation factor",1);
			}
		}

		typedef std::shared_ptr<NuclearCorrelationFactor> ncf_ptr;

		if (corrfac == "gaussslater") {
			return ncf_ptr(new GaussSlater(world, molecule));
		} else if (corrfac == "linearslater") {
			return ncf_ptr(new LinearSlater(world, molecule, a));
        } else if ((corrfac == "gradientalgaussslater") or (corrfac == "ggs")) {
            return ncf_ptr(new GradientalGaussSlater(world, molecule,a));
        } else if (corrfac == "slater") {
			return ncf_ptr(new Slater(world, molecule, a));
        } else if (corrfac == "poly4erfc") {
            return ncf_ptr(new poly4erfc(world, molecule, a));
		} else if (corrfac == "polynomial4") {
			return ncf_ptr(new Polynomial<4>(world, molecule, a ));
		} else if (corrfac == "polynomial5") {
			return ncf_ptr(new Polynomial<5>(world, molecule, a));
		} else if (corrfac == "polynomial6") {
			return ncf_ptr(new Polynomial<6>(world, molecule, a));
		} else if (corrfac == "polynomial7") {
			return ncf_ptr(new Polynomial<7>(world, molecule, a));
		} else if (corrfac == "polynomial8") {
			return ncf_ptr(new Polynomial<8>(world, molecule, a));
		} else if (corrfac == "polynomial9") {
			return ncf_ptr(new Polynomial<9>(world, molecule, a));
		} else if (corrfac == "polynomial10") {
			return ncf_ptr(new Polynomial<10>(world, molecule, a));
		} else if ((corrfac == "none") or (corrfac == "one")) {
			return ncf_ptr(new PseudoNuclearCorrelationFactor(world,
					molecule,potentialmanager,1.0));
		} else if (corrfac == "two") {
			return ncf_ptr(new PseudoNuclearCorrelationFactor(world,
					molecule,potentialmanager,2.0));
		} else if (corrfac == "linear") {
			return ncf_ptr(new PseudoNuclearCorrelationFactor(world,
					molecule,potentialmanager, a));
		} else {
			if (world.rank()==0) print(inputline);
			MADNESS_EXCEPTION("unknown nuclear correlation factor", 1);
			return ncf_ptr();
		}
	}

	std::shared_ptr<NuclearCorrelationFactor>
	create_nuclear_correlation_factor(World& world,
			const Molecule& molecule,
			const std::shared_ptr<PotentialManager> pm,
			const std::pair<std::string,double>& ncf) {
		std::stringstream ss;
		ss << ncf.first << " " << ncf.second;
		return create_nuclear_correlation_factor(world,molecule,pm,ss.str());
	}


}

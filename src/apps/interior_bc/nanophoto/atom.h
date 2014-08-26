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

  $Id: test_problems.h 1856 2010-04-06 14:03:52Z mgr522 $
*/

/** \file atom.h
    \brief Provides basis set information for computing the molecular density.

*/

#ifndef MADNESS_INTERIOR_BC_ATOM_H__INCLUDED
#define MADNESS_INTERIOR_BC_ATOM_H__INCLUDED

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include "basisfunction.h"
#include <string>

using namespace madness;

typedef std::shared_ptr<GaussianBF> BasisFunc;

/** \brief Abstract Atom class. */
class Atom {
    protected:
        Atom() {}

        std::vector<BasisFunc> basis;
        Vector<double, 3> center;

    public:
        /// \brief Sets up the basis functions for the atom.
        Atom(const Vector<double, 3> &c)
            : basis(0) {

            center[0] = c[0];
            center[1] = c[1];
            center[2] = c[2];
        }

        Atom(const Vector<double, 3> &c, int n)
            : basis(n) {

            center[0] = c[0];
            center[1] = c[1];
            center[2] = c[2];
        }

        virtual ~Atom() {
            for(std::vector<BasisFunc>::iterator iter = basis.begin();
                iter != basis.end();
                ++iter) {

                iter->reset();
            }
        }

        BasisFunc getBasisFunc(unsigned int n) {
            MADNESS_ASSERT(n >= 0 && n < basis.size());
            return basis[n];
        }

        const Vector<double, 3> &getCenter() const {
            return center;
        }

        virtual int dimBasis() const = 0;
};

/** \brief Hydrogen atom */
class Hydrogen : public Atom {
    private:
        Hydrogen() : Atom() {}

    public:
        Hydrogen(Vector<double, 3> &c)
            : Atom(c, 2) {

            std::vector<double> c1(3), c2(1);
            std::vector<double> e1(3), e2(1);

            c1[0] = 0.033494604338;
            c1[1] = 0.234726953484;
            c1[2] = 0.813757326146;
            e1[0] = 18.7311370;
            e1[1] = 2.8253944;
            e1[2] = 0.6401217;

            c2[0] = 1.0;
            e2[0] = 0.1612778;

            basis[0] = BasisFunc(new SBF(c1, e1, center));
            basis[1] = BasisFunc(new SBF(c2, e2, center));
        }

        virtual int dimBasis() const { return 2; }
};

#endif // MADNESS_INTERIOR_BC_ATOM_H__INCLUDED

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
    \brief

*/

#ifndef MADNESS_INTERIOR_BC_ATOM_H__INCLUDED
#define MADNESS_INTERIOR_BC_ATOM_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "basisfunction.h"
#include <string>

using namespace madness;

typedef SharedPtr<GaussianBF> BasisFunc;

/** \brief Abstract Atom class. */
class Atom {
    private:
        Atom() {}

    protected:
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

                delete *iter;
            }
        }

        BasisFunc getBasisFunc(unsigned int n) {
            assert(n >= 0 && n < basis.size());
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
        Hydrogen() : Atom(0) {}

    public:
        Hydrogen(Vector<double, 3> &c)
            : Atom(c, 2) {

            std::vector<double> c1(3), c2(1);
            std::vector<double> e1(3), e2(1);

            c1[0] = 0.033494604338;
            c1[1] = 0.234726953484;
            c1[2] = 0.813757326146;
            e1[0] = 13.0077340;
            e1[1] = 1.9620794;
            e1[2] = 0.4445290;

            c2[0] = 1.0;
            e2[0] = 0.1219492;

            basis[0] = BasisFunc(new SBF(c1, e1, center));
            basis[1] = BasisFunc(new SBF(c2, e2, center));
        }

        virtual int dimBasis() const { return 2; }
};

/** \brief Carbon atom */
class Carbon : public Atom {
    private:
        Carbon() : Atom(0) {}

    public:
        Carbon(Vector<double, 3> &c)
            : Atom(c, 8) {

            std::vector<double> c1s(3), c1p(3), e1(3);
            std::vector<double> c2s(1), c2p(1), e2(1);

            c1s[0] = -0.201519462736;
            c1s[1] = 0.111217608663;
            c1s[2] = 0.976798591283;
            c1p[0] = 0.129228777258;
            c1p[1] = 0.415605896012;
            c1p[2] = 0.607426909344;
            e1[0] = 4.2860000;
            e1[1] = 1.0460000;
            e1[2] = 0.3447000;

            c2s[0] = 1.0;
            c2p[0] = 1.0;
            e2[0] = 0.1128000;

            basis[0] = BasisFunc(new SBF(c1s, e1, center));
            basis[1] = BasisFunc(new PBF(c1p, e1, center, PBF::X));
            basis[2] = BasisFunc(new PBF(c1p, e1, center, PBF::Y));
            basis[3] = BasisFunc(new PBF(c1p, e1, center, PBF::Z));
            basis[4] = BasisFunc(new SBF(c2s, e2, center));
            basis[5] = BasisFunc(new PBF(c2p, e2, center, PBF::X));
            basis[6] = BasisFunc(new PBF(c2p, e2, center, PBF::Y));
            basis[7] = BasisFunc(new PBF(c2p, e2, center, PBF::Z));
        }

        virtual int dimBasis() const { return 8; }
};

/** \brief Silicon atom */
class Silicon : public Atom {
    private:
        Silicon() : Atom(0) {}

    public:
        Silicon(Vector<double, 3> &c)
            : Atom(c, 8) {

            std::vector<double> c1s(3), c1p(3), e1(3);
            std::vector<double> c2s(1), c2p(1), e2(1);

            c1s[0] = -0.422951019072;
            c1s[1] = 0.240668175467;
            c1s[2] = 1.014688250150;
            c1p[0] = -0.118957567934;
            c1p[1] = 0.334854995370;
            c1p[2] = 0.795847246205;
            e1[0] = 1.1670000;
            e1[1] = 0.5268000;
            e1[2] = 0.1807000;

            c2s[0] = 1.0;
            c2p[0] = 1.0;
            e2[0] = 0.0648000;

            basis[0] = BasisFunc(new SBF(c1s, e1, center));
            basis[1] = BasisFunc(new PBF(c1p, e1, center, PBF::X));
            basis[2] = BasisFunc(new PBF(c1p, e1, center, PBF::Y));
            basis[3] = BasisFunc(new PBF(c1p, e1, center, PBF::Z));
            basis[4] = BasisFunc(new SBF(c2s, e2, center));
            basis[5] = BasisFunc(new PBF(c2p, e2, center, PBF::X));
            basis[6] = BasisFunc(new PBF(c2p, e2, center, PBF::Y));
            basis[7] = BasisFunc(new PBF(c2p, e2, center, PBF::Z));
        }

        virtual int dimBasis() const { return 8; }
};

#endif // MADNESS_INTERIOR_BC_ATOM_H__INCLUDED

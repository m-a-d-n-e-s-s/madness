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

/// \file molsys.h
/// \brief represents molecular system
/// \defgroup moldft The molecular density funcitonal and Hartree-Fock code

#ifndef MADNESS_MOLSYS_H
#define MADNESS_MOLSYS_H

#include <madness/mra/mra.h>
#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>

struct MolecularSystem {
    // typedefs
    typedef madness::Function<double,3> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef madness::Tensor<double> tensorT;

    // fields
    bool spin_restricted;           // true if restricted, and information of beta orbitals are not keeped
    vecfuncT amo, bmo;              // Molecular orbitals (alpha/beta)
    std::vector<int> aset, bset;    // degenerate (or same principal quantum number?) set identifier (alpha/beta)
    tensorT aocc, bocc;             // occupation number of electron of orbitals (alpha/beta)
    tensorT aeps, beps;             // eigenvalues (alpha/beta)
    int nio;                        // Number of Archive I/O


    // member functions

    // constructor
    MolecularSystem () {};

    // save/load with restartdata
    void save (World & world, const char * filename);
    void load (World & world, const char * filename, unsigned int naplha, unsigned int nvalpha, unsigned int nbeta, unsigned int nvbeta, unsigned int n_core);
};

#endif


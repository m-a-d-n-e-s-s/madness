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


#ifndef MADNESS_RESPONSE_H
#define MADNESS_RESPONSE_H

/// \file response.h
/// \brief Coupled-Purturbed HF/KS

#include <mra/mra.h>
#include <vector>
#include <tkato/scf.h>

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr<operatorT> poperatorT;


struct PurturbationOperator
{
    World & world;
    functionT & op;
    int axis;
    double freq;

    PurturbationOperator (World & world, functionT op, int axis, double freq)
        : world(world), op(op), axis(axis), freq(freq) {}
};

struct CoupledPurturbation
{
    World & world;          ///< World object
    SCF & calc;     ///< SCF object including solved MOs with HF/DFT
    double freq;            ///< frequency
    vector< vector<vecfuncT> > rmo;     ///< virtual orbitals
    vector< vector<vecfuncT> > Vrmo;    ///< potential applied virtual orbitals
    int nXY;                ///< Tamm-Dancoff 1: ON 2: OFF

    CoupledPurturbation (World & world, Calculation & calc, double freq, bool enableTD = false)
        : world(world), calc(calc), freq(freq), rmo(3), Vrmo(3), nXY(enableTD?1:2) {}

    void guess_excite (int axis);
    void project_diag_space (int axis);
    vecfuncT calcBSH (vecfuncT & rmo, double eval);
    double max_norm (int axis);
    void make_vpsi (int axis);
    void solve (int axis, PurturbationOperator & purturb);
    void dipole_response ();
};

#endif


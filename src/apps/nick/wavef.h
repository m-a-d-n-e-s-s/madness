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
#ifndef WAVEF_H
#define WAVEF_H
//\file wavef.h
//\brief The hydrogenic bound and continuum states
/************************************************************************
 * Here is a madness representation of the hydrogenic wave functions.
 * The bound states come from the Gnu Scientific Library. The unbound
 * states are generated with the confluent hypergeometric function which
 * uses gmp and mpfr for extended precision
 * 
 * Using: Gnu Scientific Library          http://www.gnu.org/software/gsl/
 *        GNU Multiple Precision library  http://gmplib.org/
 *        mpfr                            http://www.mpfr.org/
 * By:    Nick Vence
 ************************************************************************/
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include "interp.h"
#include <complex>
#include <iostream>
#include <stdio.h>

const int NDIM  = 3;
typedef std::complex<double> complexd;
typedef madness::SharedPtr< madness::FunctionFunctorInterface<complexd,NDIM> > functorT;
typedef madness::Vector<double,NDIM> vector3D;
const complexd I(0,1);
const double PI = M_PI;
const complexd one(1,0);
void     testFact(madness::World&);
complexd gamma(complexd AA);
complexd gamma(double re, double im);
void     testGamma(madness::World&);


class baseWF : public madness::FunctionFunctorInterface<complexd,NDIM> {
public:
    typedef madness::Vector<double,NDIM> vector3D;
    complexd operator()(const vector3D& x) const = 0;
};

/******************************************
 * Scattering WaveFunction
 ******************************************/
class ScatteringWF : public baseWF { 
public:
    ScatteringWF(const double Z, const vector3D& kVec );
    ScatteringWF(const double Z, const vector3D& kVec, const double cutoff );
    ScatteringWF(madness::World& world, const double Z, const vector3D& kVec, const double cutoff );
    complexd operator()(const vector3D& x) const;
    complexd aForm3(complexd ZZ) const;
    complexd f11(double x) const;
    double   Z;
    vector3D kVec;
    double   k;
    double   domain;
    double   ra;
    double   cutoff;
private:
    CubicInterpolationTable<complexd > fit1F1;
    complexd expmPIZ_k;
    complexd expPIZ_2k;
    complexd expPIZ_2kXgamma1pIZ_k;
    complexd gamma1pIZ_k;
    complexd gammamIZ_k;
    complexd one;
    double   dx;
    int n;
    struct MemberFuncPtr {
        ScatteringWF* obj;
        MemberFuncPtr(ScatteringWF* obj) : obj(obj) {}
        complexd operator()(double x) {return obj->f11(x);}
    };
};

class phikl : public ScatteringWF {
public:
    complexd operator()(const vector3D& x) const;
};    

class phiK : public ScatteringWF {
public:
    complexd operator()(const vector3D& x) const;
};    

/******************************************
 * Bound WaveFunction
 ******************************************/
class BoundWF : public madness::FunctionFunctorInterface<complexd,NDIM> {
public:
    typedef madness::Vector<double,NDIM> vector3D;
    BoundWF(double Z, int nn, int ll, int mm );
    complexd operator()(const vector3D& x) const;
private:
    double Z;
    int n;
    int l;
    int m;
};

/******************************************
 *Exp[ I*(k.r) ]
 ******************************************/
class Expikr : public madness::FunctionFunctorInterface<complexd,NDIM>
{
public:
    typedef madness::Vector<double,NDIM> vector3D;
    Expikr(const vector3D& kVec);
    complexd operator()(const vector3D& r) const;
private:
    vector3D kVec;
    double k;
    double costhK;
};

/******************************************
 *Exp[ -I*(kr + k.r) ]
 ******************************************/
class Expikr2 : public madness::FunctionFunctorInterface<complexd,NDIM>
{
public:
    typedef madness::Vector<double,NDIM> vector3D;
    Expikr2(const vector3D& kVec);
    complexd operator()(const vector3D& r) const;
private:
    vector3D kVec;
    double k;
    double costhK;
};
#endif

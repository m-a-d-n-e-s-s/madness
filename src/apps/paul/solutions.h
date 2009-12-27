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
// various solutions to the incompressible navier-stokes equations

#include<mra/mra.h> 

namespace madness {

typedef Vector<double, 3> coordT3d;
typedef Vector<double, 3> coordT3d;
typedef Vector<double, 1> coordT1d;
typedef Function<double, 3> functionT;
typedef std::vector<functionT> functT;

const double mu = 1; // Effective Viscosity

double mytime = 0.0; // Global variable for the current time

const double Lx = 2. * WST_PI;
const double Ly = 1. * WST_PI;
const double Lz = 1. * WST_PI;

// This should be passed in thru the class or app context

//*****************************************************************************
static double init_zero(const coordT3d& r) {
	return 0.0;
}
//*****************************************************************************
//


//*****************************************************************************
//*****************************************************************************
//  Example # 1 from the Wuhan Paper - Analytical Cosines example
//*****************************************************************************
//*****************************************************************************
/*
static double uxexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(x) * sin(x) * (sin(2. * y) * sin(z) * sin(z) - sin(y)
			* sin(y) * sin(2. * z));
}
static double uyexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(y) * sin(y) * (sin(2. * z) * sin(x) * sin(x) - sin(z)
			* sin(z) * sin(2. * x));
}
static double uzexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(z) * sin(z) * (sin(2. * x) * sin(y) * sin(y) - sin(x)
			* sin(x) * sin(2. * y));
}
static double fxexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

 double value =     -sin(t) * pow(sin(x), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) - cos(t) * sin(x) * sin(y) * cos(z) + 0.2e1 * pow(cos(t), 0.2e1) * pow(sin(x), 0.3e1) * pow(sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z), 0.2e1) * cos(x) + pow(cos(t), 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x)) * pow(sin(x), 0.2e1) * (0.2e1 * cos(0.2e1 * y) * pow(sin(z), 0.2e1) - 0.2e1 * sin(y) * sin(0.2e1 * z) * cos(y)) + pow(cos(t), 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y)) * pow(sin(x), 0.2e1) * (0.2e1 * sin(0.2e1 * y) * sin(z) * cos(z) - 0.2e1 * pow(sin(y), 0.2e1) * cos(0.2e1 * z)) - mu * (0.2e1 * cos(t) * pow(cos(x), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) - 0.2e1 * cos(t) * pow(sin(x), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) + cos(t) * pow(sin(x), 0.2e1) * (-0.4e1 * sin(0.2e1 * y) * pow(sin(z), 0.2e1) - 0.2e1 * pow(cos(y), 0.2e1) * sin(0.2e1 * z) + 0.2e1 * pow(sin(y), 0.2e1) * sin(0.2e1 * z)) + cos(t) * pow(sin(x), 0.2e1) * (0.2e1 * sin(0.2e1 * y) * pow(cos(z), 0.2e1) - 0.2e1 * sin(0.2e1 * y) * pow(sin(z), 0.2e1) + 0.4e1 * pow(sin(y), 0.2e1) * sin(0.2e1 * z)));

	/*
	2*sin(t)*pow(sin(x),2)*sin(y)*sin(y - z)*sin(z) -
   cos(t)*(2*pow(cos(z),2)*pow(sin(x),2)*sin(2*y) +
      cos(z)*sin(x)*(sin(y) + 2*(3 - 5*cos(2*y))*sin(x)*sin(z)) -
      4*sin(z)*(pow(cos(x),2)*sin(y)*sin(y - z) + 2*pow(sin(x),2)*sin(2*y)*sin(z))) ;
*/
/*
  return  value ;
}
static double fyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double value =      -sin(t) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x)) + cos(t) * cos(x) * cos(y) * cos(z) + pow(cos(t), 0.2e1) * pow(sin(x), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) * pow(sin(y), 0.2e1) * (0.2e1 * sin(0.2e1 * z) * sin(x) * cos(x) - 0.2e1 * pow(sin(z), 0.2e1) * cos(0.2e1 * x)) + 0.2e1 * pow(cos(t), 0.2e1) * pow(sin(y), 0.3e1) * pow(sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x), 0.2e1) * cos(y) + pow(cos(t), 0.2e1) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y)) * pow(sin(y), 0.2e1) * (0.2e1 * cos(0.2e1 * z) * pow(sin(x), 0.2e1) - 0.2e1 * sin(z) * sin(0.2e1 * x) * cos(z)) - mu * (cos(t) * pow(sin(y), 0.2e1) * (0.2e1 * sin(0.2e1 * z) * pow(cos(x), 0.2e1) - 0.2e1 * sin(0.2e1 * z) * pow(sin(x), 0.2e1) + 0.4e1 * pow(sin(z), 0.2e1) * sin(0.2e1 * x)) + 0.2e1 * cos(t) * pow(cos(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x)) - 0.2e1 * cos(t) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x)) + cos(t) * pow(sin(y), 0.2e1) * (-0.4e1 * sin(0.2e1 * z) * pow(sin(x), 0.2e1) - 0.2e1 * pow(cos(z), 0.2e1) * sin(0.2e1 * x) + 0.2e1 * pow(sin(z), 0.2e1) * sin(0.2e1 * x)));

	
	/*-2*sin(t)*sin(x)*pow(sin(y),2)*sin(x - z)*sin(z) +
   cos(t)*(cos(x)*cos(y)*cos(z) + 2*pow(cos(z),2)*sin(2*x)*pow(sin(y),2) +
      2*(3 - 5*cos(2*y))*sin(x)*sin(x - z)*sin(z) -
      2*pow(cos(x),2)*pow(sin(y),2)*sin(2*z)) ;
*/
/*
  return  value ;
}
static double fzexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double value =        -sin(t) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y)) - cos(t) * cos(x) * sin(y) * sin(z) + pow(cos(t), 0.2e1) * pow(sin(x), 0.2e1) * (sin(0.2e1 * y) * pow(sin(z), 0.2e1) - pow(sin(y), 0.2e1) * sin(0.2e1 * z)) * pow(sin(z), 0.2e1) * (0.2e1 * cos(0.2e1 * x) * pow(sin(y), 0.2e1) - 0.2e1 * sin(x) * sin(0.2e1 * y) * cos(x)) + pow(cos(t), 0.2e1) * pow(sin(y), 0.2e1) * (sin(0.2e1 * z) * pow(sin(x), 0.2e1) - pow(sin(z), 0.2e1) * sin(0.2e1 * x)) * pow(sin(z), 0.2e1) * (0.2e1 * sin(0.2e1 * x) * sin(y) * cos(y) - 0.2e1 * pow(sin(x), 0.2e1) * cos(0.2e1 * y)) + 0.2e1 * pow(cos(t), 0.2e1) * pow(sin(z), 0.3e1) * pow(sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y), 0.2e1) * cos(z) - mu * (cos(t) * pow(sin(z), 0.2e1) * (-0.4e1 * sin(0.2e1 * x) * pow(sin(y), 0.2e1) - 0.2e1 * pow(cos(x), 0.2e1) * sin(0.2e1 * y) + 0.2e1 * pow(sin(x), 0.2e1) * sin(0.2e1 * y)) + cos(t) * pow(sin(z), 0.2e1) * (0.2e1 * sin(0.2e1 * x) * pow(cos(y), 0.2e1) - 0.2e1 * sin(0.2e1 * x) * pow(sin(y), 0.2e1) + 0.4e1 * pow(sin(x), 0.2e1) * sin(0.2e1 * y)) + 0.2e1 * cos(t) * pow(cos(z), 0.2e1) * (sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y)) - 0.2e1 * cos(t) * pow(sin(z), 0.2e1) * (sin(0.2e1 * x) * pow(sin(y), 0.2e1) - pow(sin(x), 0.2e1) * sin(0.2e1 * y)));

	/*
	2*sin(t)*sin(x)*sin(x - y)*sin(y)*pow(sin(z),2) +
   cos(t)*(4*pow(cos(z),2)*sin(x)*sin(x - y)*sin(y) +
      sin(z)*(-(cos(x)*sin(y)) + (3*sin(2*x) - 5*sin(2*(x - y)) - 3*sin(2*y))*sin(z))) ;
*/
/*
  return  value ;
}
static double pexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * cos(x) * sin(y) * cos(z);
	//return cos(2.*WST_PI*y)*cos(2.*WST_PI*z) ;
}
//*****************************************************************************


//*****************************************************************************
//*****************************************************************************
//  Example # 2 from the Wuhan Paper
//*****************************************************************************
/*
static double pexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value = cos(t)*cos(Xx)*cos(Zz)*sin(Yy)  ;
     return value ;
}
static double uxexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value = cos(t)*pow(sin(Xx),2)*(sin(2*Yy)*pow(sin(Zz),2) - 
                    pow(sin(Yy),2)*sin(2*Zz)) ;
     return value ;
}
static double uyexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value = cos(t)*pow(sin(Yy),2)*(-(sin(2*Xx)*pow(sin(Zz),2)) + 
                    pow(sin(Xx),2)*sin(2*Zz)) ;
     return value ;
}
static double uzexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value = cos(t)*(sin(2*Xx)*pow(sin(Yy),2) - pow(sin(Xx),2)*sin(2*Yy))*
                    pow(sin(Zz),2) ;
     return value ;
}
static double fxexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value =-(cos(t)*cos(Zz)*sin(Xx)*sin(Yy)) + 
   pow(cos(t),2)*pow(sin(Xx),2)*
    (sin(2*Xx)*pow(sin(Yy),2) - pow(sin(Xx),2)*sin(2*Yy))*
    pow(sin(Zz),2)*(-2*cos(2*Zz)*pow(sin(Yy),2) + 
      2*cos(Zz)*sin(2*Yy)*sin(Zz)) + 
   pow(cos(t),2)*pow(sin(Xx),2)*pow(sin(Yy),2)*
    (-(sin(2*Xx)*pow(sin(Zz),2)) + pow(sin(Xx),2)*sin(2*Zz))*
    (2*cos(2*Yy)*pow(sin(Zz),2) - 2*cos(Yy)*sin(Yy)*sin(2*Zz)) - 
   2*cos(t)*pow(cos(Xx),2)*(sin(2*Yy)*pow(sin(Zz),2) - 
      pow(sin(Yy),2)*sin(2*Zz)) + 
   2*cos(t)*pow(sin(Xx),2)*(sin(2*Yy)*pow(sin(Zz),2) - 
      pow(sin(Yy),2)*sin(2*Zz)) - 
   sin(t)*pow(sin(Xx),2)*(sin(2*Yy)*pow(sin(Zz),2) - 
      pow(sin(Yy),2)*sin(2*Zz)) + 
   2*pow(cos(t),2)*cos(Xx)*pow(sin(Xx),3)*
    pow(sin(2*Yy)*pow(sin(Zz),2) - pow(sin(Yy),2)*sin(2*Zz),2) - 
   cos(t)*pow(sin(Xx),2)*(-4*sin(2*Yy)*pow(sin(Zz),2) - 
      2*pow(cos(Yy),2)*sin(2*Zz) + 2*pow(sin(Yy),2)*sin(2*Zz)) - 
   cos(t)*pow(sin(Xx),2)*(2*pow(cos(Zz),2)*sin(2*Yy) - 
      2*sin(2*Yy)*pow(sin(Zz),2) + 4*pow(sin(Yy),2)*sin(2*Zz)) ; 

     return value ;
}
static double fyexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value =  cos(t)*cos(Xx)*cos(Yy)*cos(Zz) + 
   pow(cos(t),2)*pow(sin(Yy),2)*
    (sin(2*Xx)*pow(sin(Yy),2) - pow(sin(Xx),2)*sin(2*Yy))*
    pow(sin(Zz),2)*(2*cos(2*Zz)*pow(sin(Xx),2) - 
      2*cos(Zz)*sin(2*Xx)*sin(Zz)) - 
   cos(t)*pow(sin(Yy),2)*(-2*pow(cos(Zz),2)*sin(2*Xx) + 
      2*sin(2*Xx)*pow(sin(Zz),2) - 4*pow(sin(Xx),2)*sin(2*Zz)) - 
   cos(t)*pow(sin(Yy),2)*(4*sin(2*Xx)*pow(sin(Zz),2) + 
      2*pow(cos(Xx),2)*sin(2*Zz) - 2*pow(sin(Xx),2)*sin(2*Zz)) - 
   2*cos(t)*pow(cos(Yy),2)*(-(sin(2*Xx)*pow(sin(Zz),2)) + 
      pow(sin(Xx),2)*sin(2*Zz)) + 
   2*cos(t)*pow(sin(Yy),2)*(-(sin(2*Xx)*pow(sin(Zz),2)) + 
      pow(sin(Xx),2)*sin(2*Zz)) - 
   sin(t)*pow(sin(Yy),2)*(-(sin(2*Xx)*pow(sin(Zz),2)) + 
      pow(sin(Xx),2)*sin(2*Zz)) + 
   2*pow(cos(t),2)*cos(Yy)*pow(sin(Yy),3)*
    pow(-(sin(2*Xx)*pow(sin(Zz),2)) + pow(sin(Xx),2)*sin(2*Zz),2)\
    + pow(cos(t),2)*pow(sin(Xx),2)*pow(sin(Yy),2)*
    (-2*cos(2*Xx)*pow(sin(Zz),2) + 2*cos(Xx)*sin(Xx)*sin(2*Zz))*
    (sin(2*Yy)*pow(sin(Zz),2) - pow(sin(Yy),2)*sin(2*Zz)) ; 

     return value ;
}
static double fzexact(const coordT3d& r) {
	const double Xx = r[0], Yy = r[1], Zz = r[2];
	double t = mytime;

     double value =  -2*cos(t)*pow(cos(Zz),2)*(sin(2*Xx)*pow(sin(Yy),2) - 
      pow(sin(Xx),2)*sin(2*Yy)) - cos(t)*cos(Xx)*sin(Yy)*sin(Zz) + 
   2*cos(t)*(sin(2*Xx)*pow(sin(Yy),2) - pow(sin(Xx),2)*sin(2*Yy))*
    pow(sin(Zz),2) - sin(t)*
    (sin(2*Xx)*pow(sin(Yy),2) - pow(sin(Xx),2)*sin(2*Yy))*
    pow(sin(Zz),2) - cos(t)*
    (-4*sin(2*Xx)*pow(sin(Yy),2) - 2*pow(cos(Xx),2)*sin(2*Yy) + 
      2*pow(sin(Xx),2)*sin(2*Yy))*pow(sin(Zz),2) - 
   cos(t)*(2*pow(cos(Yy),2)*sin(2*Xx) - 2*sin(2*Xx)*pow(sin(Yy),2) + 
      4*pow(sin(Xx),2)*sin(2*Yy))*pow(sin(Zz),2) + 
   2*pow(cos(t),2)*cos(Zz)*pow(sin(2*Xx)*pow(sin(Yy),2) - 
      pow(sin(Xx),2)*sin(2*Yy),2)*pow(sin(Zz),3) + 
   pow(cos(t),2)*pow(sin(Yy),2)*
    (-2*cos(2*Yy)*pow(sin(Xx),2) + 2*cos(Yy)*sin(2*Xx)*sin(Yy))*
    pow(sin(Zz),2)*(-(sin(2*Xx)*pow(sin(Zz),2)) + 
      pow(sin(Xx),2)*sin(2*Zz)) + 
   pow(cos(t),2)*pow(sin(Xx),2)*
    (2*cos(2*Xx)*pow(sin(Yy),2) - 2*cos(Xx)*sin(Xx)*sin(2*Yy))*
    pow(sin(Zz),2)*(sin(2*Yy)*pow(sin(Zz),2) - 
      pow(sin(Yy),2)*sin(2*Zz)) ;  

     return value ;
}
*/


// two interacting vortices

static double uxexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(x) * sin(x) * (sin(2. * y) * sin(z) * sin(z) - sin(y)
			* sin(y) * sin(2. * z));
	//return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double uyexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(y) * sin(y) * (sin(2. * z) * sin(x) * sin(x) - sin(z)
			* sin(z) * sin(2. * x));
	//return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double uzexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

	return cos(t) * sin(z) * sin(z) * (sin(2. * x) * sin(y) * sin(y) - sin(x)
			* sin(x) * sin(2. * y));
	//return  t*cos(WST_PI*y)*cos(WST_PI*z) ;
}
static double fxexact(const coordT3d& r)
{
  return 0. ;
}
static double fyexact(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
        double t = mytime ;

  double x0 = Lx / 2. ;
  double y0 = Ly / 2. ;
  double Lparam = 0.25 ;

  // Creates two interacting vortices
     double myval =  18.* exp(- ( (x-x0)*(x-x0) + (y-y0)*(y-y0) ) / Lparam)  ;
  return myval ;
}
static double fzexact(const coordT3d& r)
{
  return 0. ;
}

static double pexact(const coordT3d& r) {
	const double x = r[0], y = r[1], z = r[2];
	double t = mytime;

        return cos(t) * cos(x) * sin(y) * cos(z);
}


//*****************************************************************************
//*****************************************************************************

inline functionT div(const functT& uint) {
	return diff(uint[0], 0) + diff(uint[1], 1) + diff(uint[2], 2);
}

inline functionT lap(const functionT& uint) {
	return diff(diff(uint, 0),0) + diff(diff(uint,1), 1) + diff(diff(uint,2), 2);
}

}


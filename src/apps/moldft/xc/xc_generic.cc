/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
  $Id: funcimpl.h 465 2008-01-30 22:29:17Z wsttiger $
*/


/// \file xc_generic.cc
/// \brief wrapper for DL XC functional library

#include "./dlxc.h"

/*
  Here is some info on the DL XC interface design extracted from 

  http://www.cse.scitech.ac.uk/ccg/dft/design.html

  Please refer to that page for additional details of the XC interface

    * The level of derivatives required (ideriv valid values are 0,1, and 2)
    * The number of grid points (npt)
    * The alpha density (rhoa)
    * The beta density (rhob)
    * The dot product of the alpha density gradient with itself (sigmaaa)
    * The dot product of the beta density gradient with itself (sigmabb)
    * The dot product of the alpha density gradient with the beta density gradient (sigmaab) 

and return

    * The value of the functional (zk)
    * The derivative of the functional with respect to the alpha density (vrhoa)
    * The derivative of the functional with respect to the beta density (vrhob)
    * The derivative of the functional with respect to gamma-alpha-alpha (vsigmaaa)
    * The derivative of the functional with respect to gamma-beta-beta (vsigmabb)
    * The derivative of the functional with respect to gamma-alpha-beta (vsigmaab)
    * The 2nd derivative of the functional with respect to the alpha density (v2rhoa2)
    * The 2nd derivative of the functional with respect to the beta density (v2rhob2)
    * The 2nd derivative of the functional with respect to the alpha and beta density (v2rhoab)
    * The 2nd derivative of the functional with respect to rho-alpha and gamma-alpha-alpha (v2rhoasigmaaa)
    * The 2nd derivative of the functional with respect to rho-alpha and gamma-alpha-beta (v2rhoasigmaab)
    * The 2nd derivative of the functional with respect to rho-alpha and gamma-beta-beta (v2rhoasigmabb)
    * The 2nd derivative of the functional with respect to rho-beta and gamma-beta-beta (v2rhobsigmabb)
    * The 2nd derivative of the functional with respect to rho-beta and gamma-alpha-beta (v2rhobsigmaab)
    * The 2nd derivative of the functional with respect to rho-beta and gamma-alpha-alpha (v2rhobsigmaaa)
    * The 2nd derivative of the functional with respect to gamma-alpha-alpha (v2sigmaaa2)
    * The 2nd derivative of the functional with respect to gamma-alpha-alpha and gamma-alpha-beta (v2sigmaaaab)
    * The 2nd derivative of the functional with respect to gamma-alpha-alpha and gamma-beta-beta (v2sigmaaabb)
    * The 2nd derivative of the functional with respect to gamma-alpha-beta (v2sigmaab2)
    * The 2nd derivative of the functional with respect to gamma-alpha-beta and gamma-beta-beta (v2sigmaabbb)
    * The 2nd derivative of the functional with respect to gamma-beta-beta (v2sigmabb2) 

where the name in brackets is the variable to be used in the codes to
represent the quantities.

In the closed shell case however less quantities need to be
computed. Also some of the quantities one would expect can be combined
at the functional level. The closed shell case routines only return
these compound quantities which are defined as closed shell quantity
equivalent in terms of open shell quantities

vsigmaaa	:=	2*vsigmaaa+vsigmaab

v2rhoa2	:=	v2rhoa2+v2rhoab

v2rhoasigmaaa	:=	v2rhoasigmaaa+v2rhoasigmaab+v2rhoasigmabb

v2sigmaaa2	:=	2*v2sigmaaa2+4*v2sigmaaaab+2*v2sigmaaabb+v2sigmaab2

 */


typedef int (*DL_XC_RKS) (integer *ideriv, integer *npt, doublereal *rhoa1, 
                          doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
                          doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
                          doublereal *v2sigmaaa2);


static const double THRESH_RHO = 1e-12;
static const double THRESH_GRHO = 1e-20;

static void munge_grho(int npoint, double *rho, double *grho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
        if ((rho[i] <=THRESH_RHO) || 
            (grho[i] < THRESH_GRHO)) grho[i] = THRESH_GRHO;            
    }
}

static void munge_rho(int npoint, double *rho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
    }
}


/// Computes and combines the list of closed shell functionals and potentials

/// Also munges the density and its derivative to avoid numerical grief.
///
/// If you are doing just LDA you should call xc_rks_generic_lda instead
/// since that doees not require gamma.
///
/// The result arrays are assume initialized before entry and this routine 
/// accumulates into them.
///
/// THIS ROUTINE IS NOT YET TESTED!
void xc_rks_generic(const std::vector<DL_XC_RKS>& functionals, ///< Functionals to invoke
                    const std::vector<double>& c,              ///< Coefficients for functionals
                    const Tensor<double>& rho_alpha,           ///< Alpha-spin density at each grid point
                    const Tensor<double>& gamma_alpha,         ///< |del rho_alpha|^2
                    Tensor<double>& f,                         ///< Value of functional at each grid point
                    Tensor<double>& df_drho,                   ///< Derivative of functional w.r.t. rho_alpha
                    Tensor<double>& df_dgamma,                 ///< Derivative of functional w.r.t. gamma_alpha
                    Tensor<double>& d2f_drho2,                 ///< 2nd derivative of functional 
                    Tensor<double>& d2f_drhodgamma,            ///< 2nd derivative of functional 
                    Tensor<double>& d2f_dgamma2)               ///< 2nd derivative of functional 
{
    integer ideriv = 2;
    integer npt = rho_alpha.dim[0];
    
    for (unsigned int i=0; i<functionals.size(); i++) {
        Tensor<double> tf(npt);
        Tensor<double> tdf_drho(npt);
        Tensor<double> tdf_dgamma(npt);
        Tensor<double> td2f_drho2(npt);
        Tensor<double> td2f_drhodgamma(npt);
        Tensor<double> td2f_dgamma2(npt);
        
        functionals[i](&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(), 
                       tf.ptr(), 
                       tdf_drho.ptr(), tdf_dgamma.ptr(), 
                       td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());

        f.gaxpy(1.0, tf, c[i]);
        df_drho.gaxpy(1.0, tdf_drho, c[i]);
        df_dgamma.gaxpy(1.0, tdf_dgamma, c[i]);
        d2f_drho2.gaxpy(1.0, td2f_drho2, c[i]);
        d2f_drhodgamma.gaxpy(1.0, td2f_drhodgamma, c[i]);
        d2f_dgamma2.gaxpy(1.0, td2f_dgamma2, c[i]);
    }
}

void xc_rks_generic_lda(const Tensor<double>& rho_alpha,           ///< Alpha-spin density at each grid point
                        Tensor<double>& f,                         ///< Value of functional at each grid point
                        Tensor<double>& df_drho)                   ///< Derivative of functional w.r.t. rho_alpha
{
    integer ideriv = 1;
    integer npt = rho_alpha.dim[0];
    
    Tensor<double> gamma_alpha(npt);
    Tensor<double> tf(npt);
    Tensor<double> tdf_drho(npt);
    Tensor<double> tdf_dgamma(npt);
    Tensor<double> td2f_drho2(npt);
    Tensor<double> td2f_drhodgamma(npt);
    Tensor<double> td2f_dgamma2(npt);
        
    xc_x_lda(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(), 
             tf.ptr(), 
             tdf_drho.ptr(), tdf_dgamma.ptr(), 
             td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());
    
    f.gaxpy(1.0, tf, 1.0);
    df_drho.gaxpy(1.0, tdf_drho, 1.0);

    tf.fill(0.0);
    tdf_drho.fill(0.0);

    xc_c_vwn5(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(), 
              tf.ptr(), 
              tdf_drho.ptr(), tdf_dgamma.ptr(), 
              td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());
    
    f.gaxpy(1.0, tf, 1.0);
    df_drho.gaxpy(1.0, tdf_drho, 1.0);
}


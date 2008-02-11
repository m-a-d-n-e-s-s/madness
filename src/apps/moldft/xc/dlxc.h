#ifndef MADNESS_DLXC_H
#define MADNESS_DLXC_H

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

#include "./f2c.h"

typedef int (*DL_XC_UKS)(integer *ideriv, integer *npt, 
                         doublereal *rhoa1, doublereal *rhob1, doublereal *sigmaaa1, doublereal *sigmabb1, 
                         doublereal *sigmaab1, doublereal *zk, doublereal *vrhoa, doublereal *
                         vrhob, doublereal *vsigmaaa, doublereal *vsigmabb, doublereal *
                         vsigmaab, doublereal *v2rhoa2, doublereal *v2rhob2, doublereal *
                         v2rhoab, doublereal *v2rhoasigmaaa, doublereal *v2rhoasigmaab, 
                         doublereal *v2rhoasigmabb, doublereal *v2rhobsigmabb, doublereal *
                         v2rhobsigmaab, doublereal *v2rhobsigmaaa, doublereal *v2sigmaaa2, 
                         doublereal *v2sigmaaaab, doublereal *v2sigmaaabb, doublereal *
                         v2sigmaab2, doublereal *v2sigmaabbb, doublereal *v2sigmabb2);

typedef int (*DL_XC_RKS) (integer *ideriv, integer *npt, doublereal *rhoa1, 
                          doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
                          doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
                          doublereal *v2sigmaaa2);

extern "C" {

    
    DL_XC_UKS uks_c_lyp__;

    DL_XC_UKS uks_c_lyp__;
    
    DL_XC_RKS_c_lyp__;
    
    DL_XC_UKS uks_c_p86__;
    
    DL_XC_RKS_c_p86__;
    
    DL_XC_UKS uks_c_pbe__;
    
    DL_XC_RKS_c_pbe__;
    
    DL_XC_UKS uks_c_pw91__;
    
    DL_XC_RKS_c_pw91__;
    
    DL_XC_UKS uks_c_pw92__;
    
    DL_XC_RKS_c_pw92__;
    
    DL_XC_UKS uks_c_pz81__;
    
    DL_XC_RKS_c_pz81__;
    
    DL_XC_UKS uks_c_vwn5__;
    
    DL_XC_RKS_c_vwn5__;
    
    DL_XC_UKS uks_c_vwn5rpa__;
    
    DL_XC_RKS_c_vwn5rpa__;
    
    DL_XC_UKS uks_x_b3__;
    
    DL_XC_RKS_x_b3__;
    
    DL_XC_UKS uks_x_b88__;
    
    DL_XC_RKS_x_b88__;
    
    DL_XC_UKS uks_xc_b3lyp__;
    
    DL_XC_RKS_xc_b3lyp__;
    
    DL_XC_UKS uks_xc_b97_1__;
    
    DL_XC_RKS_xc_b97_1__;
    
    DL_XC_UKS uks_xc_b97_2__;
    
    DL_XC_RKS_xc_b97_2__;
    
    DL_XC_UKS uks_xc_b97__;
    
    DL_XC_RKS_xc_b97__;
    
    DL_XC_UKS uks_xc_edf1__;
    
    DL_XC_RKS_xc_edf1__;
    
    DL_XC_UKS uks_xc_hcth120__;
    
    DL_XC_RKS_xc_hcth120__;
    
    DL_XC_UKS uks_xc_hcth147__;
    
    DL_XC_RKS_xc_hcth147__;
    
    DL_XC_UKS uks_xc_hcth407__;
    
    DL_XC_RKS_xc_hcth407__;
    
    DL_XC_UKS uks_xc_hcth__;
    
    DL_XC_RKS_xc_hcth__;
    
    DL_XC_UKS uks_xc_pw91__;
    
    DL_XC_RKS_xc_pw91__;
    
    DL_XC_UKS uks_x_lda__;
    
    DL_XC_RKS_x_lda__;
    
    DL_XC_UKS uks_x_pbe__;
    
    DL_XC_RKS_x_pbe__;
    
    DL_XC_UKS uks_x_pw91__;
    
    DL_XC_RKS_x_pw91__;

 
}
#endif

/*
 * libxc.h
 *
 *  Created on: Nov 23, 2008
 *      Author: wsttiger
 */

#ifndef LIBXC_H_
#define LIBXC_H_

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <world/world.h>
#include "xc.h"

using namespace madness;

const double THRESH_RHO = 1e-8;
const double THRESH_GRHO = 1e-20;

//***************************************************************************
inline void wst_munge_grho(int npoint, double *rho, double *grho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
        if ((rho[i] <=THRESH_RHO) ||
            (grho[i] < THRESH_GRHO)) grho[i] = THRESH_GRHO;
    }
}
//***************************************************************************

//***************************************************************************
inline void wst_munge_rho(int npoint, double *rho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
    }
}
//***************************************************************************

//***************************************************************************
inline void xc_generic_lda(Tensor<double> rho_alpha,           ///< Alpha-spin density at each grid point
                          Tensor<double> f,                         ///< Value of functional at each grid point
                          Tensor<double> df_drho,                   ///< Derivative of functional w.r.t. rho_alpha
                          bool spinpol)
    {
    MADNESS_ASSERT(rho_alpha.iscontiguous());
    MADNESS_ASSERT(f.iscontiguous());
    MADNESS_ASSERT(df_drho.iscontiguous());

    rho_alpha = rho_alpha.flat();
    f = f.flat();
    df_drho = df_drho.flat();

    XC(lda_type) xc_c_func;

    int npt = rho_alpha.dim[0];

    Tensor<double> tf(npt);
    Tensor<double> tdf_drho(npt);
    double* rhoptr = rho_alpha.ptr();
    double* tfptr = tf.ptr();
    double* tdf_drhoptr = tdf_drho.ptr();

    tf.fill(0.0);
    tdf_drho.fill(0.0);
    f.fill(0.0);
    df_drho.fill(0.0);

    wst_munge_rho(npt, rhoptr);

    xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
    for (int i = 0; i < npt; i++)
    {
      xc_lda_vxc(&xc_c_func, &rhoptr[i], &tfptr[i], &tdf_drhoptr[i]);
    }

    f.gaxpy(1.0, tf, 1.0);
    df_drho.gaxpy(1.0, tdf_drho, 1.0);

    tf.fill(0.0);
    tdf_drho.fill(0.0);

    xc_lda_x_init(&xc_c_func, XC_UNPOLARIZED, 3, 0);
    for (int i = 0; i < npt; i++)
    {
      xc_lda_vxc(&xc_c_func, &rhoptr[i], &tfptr[i], &tdf_drhoptr[i]);
    }

    f.gaxpy(1.0, tf, 1.0);
    df_drho.gaxpy(1.0, tdf_drho, 1.0);
}
  //***************************************************************************

  //***************************************************************************
  template <int NDIM>
 inline void xc_lda_V(const Key<NDIM>& key, Tensor<double>& t)
  {
    Tensor<double> enefunc = copy(t);
    Tensor<double> V = copy(t);
    ::xc_generic_lda(t, enefunc, V, false);
    t(___) = V(___);
  }
  //***************************************************************************

  //***************************************************************************
  template <int NDIM>
 inline void xc_lda_ene(const Key<NDIM>& key, Tensor<double>& t)
  {
    Tensor<double> V = copy(t);
    Tensor<double> enefunc = copy(t);
    ::xc_generic_lda(t, enefunc, V, false);
    t(___) = enefunc(___);
  }
  //***************************************************************************



#endif /* LIBXC_H_ */

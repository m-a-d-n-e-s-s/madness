#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "dft.h"
#include "util.h"
#include <moldft/xc/f2c.h>
#include <vector>
#include "poperator.h"

typedef madness::Vector<double,3> coordT;

////***************************************************************************
//int rks_x_lda__ (integer *ideriv, integer *npt, doublereal *rhoa1,
//                  doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
//                  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
//                  doublereal *v2sigmaaa2);
////***************************************************************************
//
////***************************************************************************
//int rks_c_vwn5__ (integer *ideriv, integer *npt, doublereal *rhoa1,
//                  doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
//                  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
//                  doublereal *v2sigmaaa2);
////***************************************************************************

/* Subroutine */ int rks_c_vwn5__(integer *ideriv, integer *npt, doublereal *
  rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
  doublereal *v2sigmaaa2)
{
  // WSTHORNTON
  static doublereal c_b2 = .16666666666666666;
  static doublereal c_b3 = .33333333333333331;
   /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), atan(
      doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, t1, t2, t4, t6, t7, t10, t20, t11, t22, t13, t23,
      t16, t17, t25, t19, t26, t28, t29, t32, t33, t37, t38, t40, t41,
      t43, t51, t52, t27, t34, t39, t46, t47, t48, t49, t53, t55, t56,
      t58, t60, t63, t66, t67, t68, t69, t77, t79, t80, t81, t88, t92,
      t94, t102, t103, t105, t107, t125, t134, t138, rho;


/*     S.H. Vosko, L. Wilk, and M. Nusair */
/*     Accurate spin-dependent electron liquid correlation energies for */
/*     local spin density calculations: a critical analysis */
/*     Can. J. Phys. 58 (1980) 1200-1211 */


/*     CITATION: */

/*     Functionals were obtained from the Density Functional Repository */
/*     as developed and distributed by the Quantum Chemistry Group, */
/*     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD */
/*     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or */
/*     Paul Sherwood for further information. */

/*     COPYRIGHT: */

/*     Users may incorporate the source code into software packages and */
/*     redistribute the source code provided the source code is not */
/*     changed in anyway and is properly cited in any documentation or */
/*     publication related to its use. */

/*     ACKNOWLEDGEMENT: */

/*     The source code was generated using Maple 8 through a modified */
/*     version of the dfauto script published in: */

/*        R. Strange, F.R. Manby, P.J. Knowles */
/*        Automatic code generation in density functional theory */
/*        Comp. Phys. Comm. 136 (2001) 310-318. */

    /* Parameter adjustments */
    --v2sigmaaa2;
    --v2rhoasigmaaa;
    --v2rhoa2;
    --vsigmaaa;
    --vrhoa;
    --zk;
    --sigmaaa1;
    --rhoa1;

    /* Function Body */
    if (*ideriv == 0) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t7 = 1 / (t2 * .6203504908994 + t4 * 2.935818660072219 +
      12.9352);
    t10 = log(t2 * .6203504908994 * t7);
    t16 = atan(6.15199081975908 / (t4 * 1.575246635799487 +
      3.72744));
/* Computing 2nd power */
    d__1 = t4 * .7876233178997433 + .10498;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    zk[i__] = rho * (t10 * .0310907 + t16 * .03878329487811301 +
      t22 * 9.690227711544374e-4);
      } else {
/* rho */
    zk[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 1) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t6 = t2 * .6203504908994 + t4 * 2.935818660072219 + 12.9352;
    t7 = 1 / t6;
    t10 = log(t2 * .6203504908994 * t7);
    t11 = t10 * .0310907;
    t13 = t4 * 1.575246635799487 + 3.72744;
    t16 = atan(6.15199081975908 / t13);
    t17 = t16 * .03878329487811301;
    t19 = t4 * .7876233178997433 + .10498;
/* Computing 2nd power */
    d__1 = t19;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    t23 = t22 * 9.690227711544374e-4;
    zk[i__] = rho * (t11 + t17 + t23);
/* Computing 2nd power */
    d__1 = t2;
    t25 = d__1 * d__1;
    t26 = 1 / t25;
/* Computing 2nd power */
    d__1 = rho;
    t28 = d__1 * d__1;
    t29 = 1 / t28;
/* Computing 2nd power */
    d__1 = t6;
    t32 = d__1 * d__1;
    t33 = 1 / t32;
/* Computing 2nd power */
    d__1 = t4;
    t37 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t37;
    t38 = d__1 * d__1;
    t40 = 1 / t38 / t4;
    t41 = t40 * t29;
    t43 = t26 * -.2067834969664667 * t29 - t41 *
      .4893031100120365;
/* Computing 2nd power */
    d__1 = t13;
    t51 = d__1 * d__1;
    t52 = 1 / t51;
    vrhoa[i__] = t11 + t17 + t23 + rho * ((t26 *
      -.2067834969664667 * t7 * t29 - t2 * .6203504908994 *
      t33 * t43) * .05011795824473985 / t2 * t6 + t52 *
      .0626408570946439 * t40 * t29 / (t52 * 37.8469910464
      + 1.) + (t19 * -.2625411059665811 * t7 * t41 - t20 *
      1. * t33 * t43) * 9.690227711544374e-4 / t20 * t6);
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 2) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t6 = t2 * .6203504908994 + t4 * 2.935818660072219 + 12.9352;
    t7 = 1 / t6;
    t10 = log(t2 * .6203504908994 * t7);
    t11 = t10 * .0310907;
    t13 = t4 * 1.575246635799487 + 3.72744;
    t16 = atan(6.15199081975908 / t13);
    t17 = t16 * .03878329487811301;
    t19 = t4 * .7876233178997433 + .10498;
/* Computing 2nd power */
    d__1 = t19;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    t23 = t22 * 9.690227711544374e-4;
    zk[i__] = rho * (t11 + t17 + t23);
/* Computing 2nd power */
    d__1 = t2;
    t25 = d__1 * d__1;
    t26 = 1 / t25;
    t27 = t26 * t7;
/* Computing 2nd power */
    d__1 = rho;
    t28 = d__1 * d__1;
    t29 = 1 / t28;
/* Computing 2nd power */
    d__1 = t6;
    t32 = d__1 * d__1;
    t33 = 1 / t32;
    t34 = t2 * t33;
/* Computing 2nd power */
    d__1 = t4;
    t37 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t37;
    t38 = d__1 * d__1;
    t39 = t38 * t4;
    t40 = 1 / t39;
    t41 = t40 * t29;
    t43 = t26 * -.2067834969664667 * t29 - t41 *
      .4893031100120365;
    t46 = t27 * -.2067834969664667 * t29 - t34 * .6203504908994 *
      t43;
    t47 = 1 / t2;
    t48 = t46 * t47;
    t49 = t48 * t6;
/* Computing 2nd power */
    d__1 = t13;
    t51 = d__1 * d__1;
    t52 = 1 / t51;
    t53 = t52 * t40;
    t55 = t52 * 37.8469910464 + 1.;
    t56 = 1 / t55;
    t58 = t53 * t29 * t56;
    t60 = t19 * t7;
    t63 = t20 * t33;
    t66 = t60 * -.2625411059665811 * t41 - t63 * 1. * t43;
    t67 = 1 / t20;
    t68 = t66 * t67;
    t69 = t68 * t6;
    vrhoa[i__] = t11 + t17 + t23 + rho * (t49 *
      .05011795824473985 + t58 * .0626408570946439 + t69 *
      9.690227711544374e-4);
    t77 = 1 / t25 / t1;
/* Computing 2nd power */
    d__1 = t28;
    t79 = d__1 * d__1;
    t80 = 1 / t79;
    t81 = t77 * t7 * t80;
    t88 = 1 / t28 / rho;
    t92 = 1 / t32 / t6;
/* Computing 2nd power */
    d__1 = t43;
    t94 = d__1 * d__1;
    t102 = 1 / t39 / t1;
    t103 = t102 * t80;
    t105 = t40 * t88;
    t107 = t77 * -.1378556646443111 * t80 + t26 *
      .4135669939329333 * t88 - t103 * .4077525916766971 +
      t105 * .978606220024073;
    t125 = t80 * t56;
/* Computing 2nd power */
    d__1 = t51;
    t134 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t55;
    t138 = d__1 * d__1;
    s1 = t49 * .2004718329789594 + t58 * .2505634283785756;
    v2rhoa2[i__] = s1 + t69 * .00387609108461775 + rho * 2. * ((
      t81 * -.1378556646443111 + t26 * .4135669939329333 *
      t33 * t29 * t43 + t27 * .4135669939329333 * t88 + t2 *
       1.2407009817988 * t92 * t94 - t34 * .6203504908994 *
      t107) * .05011795824473985 * t47 * t6 + t46 *
      .01670598608157995 / t2 / t1 * t6 * t29 + t48 *
      .05011795824473985 * t43 + .03289159980064473 / t51 /
      t13 * t77 * t125 + t52 * .05220071424553658 * t102 *
      t125 - t53 * .1252817141892878 * t88 * t56 -
      1.244848083156773 / t134 / t13 * t77 * t80 / t138 + (
      t81 * .03446391616107778 + t19 * .5250822119331622 *
      t33 * t41 * t43 - t60 * .2187842549721509 * t103 +
      t60 * .5250822119331622 * t105 + t20 * 2. * t92 * t94
      - t63 * 1. * t107) * 9.690227711544374e-4 * t67 * t6
      + t66 * 2.544083100456872e-4 / t20 / t19 * t6 * t40 *
      t29 + t68 * 9.690227711544374e-4 * t43);
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
    v2rhoa2[i__] = 0.;
      }
/* rho */
  }
    }
/* ideriv */
    return 0;
} /* rks_c_vwn5__ */

/* Subroutine */ int rks_x_lda__(integer *ideriv, integer *npt, doublereal *
  rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
  doublereal *v2sigmaaa2)
{
  // WSTHORNTON
  static doublereal c_b2 = .33333333333333331;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal t1, t5, rho;


/*     P.A.M. Dirac */
/*     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376 */


/*     CITATION: */

/*     Functionals were obtained from the Density Functional Repository */
/*     as developed and distributed by the Quantum Chemistry Group, */
/*     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD */
/*     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or */
/*     Paul Sherwood for further information. */

/*     COPYRIGHT: */

/*     Users may incorporate the source code into software packages and */
/*     redistribute the source code provided the source code is not */
/*     changed in anyway and is properly cited in any documentation or */
/*     publication related to its use. */

/*     ACKNOWLEDGEMENT: */

/*     The source code was generated using Maple 8 through a modified */
/*     version of the dfauto script published in: */

/*        R. Strange, F.R. Manby, P.J. Knowles */
/*        Automatic code generation in density functional theory */
/*        Comp. Phys. Comm. 136 (2001) 310-318. */

    /* Parameter adjustments */
    --v2sigmaaa2;
    --v2rhoasigmaaa;
    --v2rhoa2;
    --vsigmaaa;
    --vrhoa;
    --zk;
    --sigmaaa1;
    --rhoa1;

    /* Function Body */
    if (*ideriv == 0) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
      } else {
/* rho */
    zk[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 1) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
    vrhoa[i__] = t1 * -.9847450218426965;
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 2) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
    vrhoa[i__] = t1 * -.9847450218426965;
/* Computing 2nd power */
    d__1 = t1;
    t5 = d__1 * d__1;
    v2rhoa2[i__] = -.6564966812284644 / t5;
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
    v2rhoa2[i__] = 0.;
      }
/* rho */
  }
    }
/* ideriv */
    return 0;
} /* rks_x_lda__ */

namespace madness
{

const double THRESH_RHO = 1e-12;
const double THRESH_GRHO = 1e-20;

void wst_munge_grho(int npoint, double *rho, double *grho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
        if ((rho[i] <=THRESH_RHO) ||
            (grho[i] < THRESH_GRHO)) grho[i] = THRESH_GRHO;
    }
}

void wst_munge_rho(int npoint, double *rho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
    }
}

//***************************************************************************
  void xc_rks_generic_lda(Tensor<double> rho_alpha,           ///< Alpha-spin density at each grid point
                          Tensor<double> f,                         ///< Value of functional at each grid point
                          Tensor<double> df_drho)                   ///< Derivative of functional w.r.t. rho_alpha
  {
    MADNESS_ASSERT(rho_alpha.iscontiguous());
    MADNESS_ASSERT(f.iscontiguous());
    MADNESS_ASSERT(df_drho.iscontiguous());

    rho_alpha = rho_alpha.flat();
    f = f.flat();
    df_drho = df_drho.flat();

      integer ideriv = 2;
      integer npt = rho_alpha.dim[0];

      Tensor<double> gamma_alpha(npt);
      Tensor<double> tf(npt);
      Tensor<double> tdf_drho(npt);
      Tensor<double> tdf_dgamma(npt);
      Tensor<double> td2f_drho2(npt);
      Tensor<double> td2f_drhodgamma(npt);
      Tensor<double> td2f_dgamma2(npt);

      wst_munge_rho(npt, rho_alpha.ptr());

      f.fill(0.0);
      df_drho.fill(0.0);

      double* rhoptr = rho_alpha.ptr();
      for (int i = 0; i < npt; i++)
      {
        rhoptr[i] *= 2.0;
      }

      int returnvalue = ::rks_x_lda__(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(),
               tf.ptr(),
               tdf_drho.ptr(), tdf_dgamma.ptr(),
               td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());

      f.gaxpy(1.0, tf, 1.0);
      df_drho.gaxpy(1.0, tdf_drho, 1.0);

      tf.fill(0.0);
      tdf_drho.fill(0.0);

      returnvalue = ::rks_c_vwn5__(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(),
                tf.ptr(),
                tdf_drho.ptr(), tdf_dgamma.ptr(),
                td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());

      f.gaxpy(1.0, tf, 1.0);
      df_drho.gaxpy(1.0, tdf_drho, 1.0);

  }
  //***************************************************************************

  //***************************************************************************
  template <int NDIM>
  void dft_xc_lda_V(const Key<NDIM>& key, Tensor<double>& t)
  {
    Tensor<double> enefunc = copy(t);
    Tensor<double> V = copy(t);
    xc_rks_generic_lda(t, enefunc, V);
    t(___) = V(___);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearChargeDensityOp<T,NDIM>::DFTNuclearChargeDensityOp(World& world, funcT rhon,
      double coeff, double thresh, bool periodic) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearChargeDensityOp");
    _rhon = rhon;
    SeparatedConvolution<T,NDIM>* cop;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
          1e-4, thresh, L);
    }
    else
    {
      cop =
        CoulombOperatorPtr<T,NDIM>(world,
            FunctionDefaults<NDIM>::get_k(), 1e-8, thresh);
    }
    // Apply operator to get potential
    _Vnuc = apply(*cop, rhon);
    delete cop;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearPotentialOp<T,NDIM>::DFTNuclearPotentialOp(World& world, funcT V,
      double coeff, double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearPotentialOp");
    _V = V;
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombOp<T,NDIM>::DFTCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("CoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    _cop = CoulombOperatorPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-4, thresh);
    // Initialize potential
    _Vc = FunctionFactory<T,NDIM>(world);
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombPeriodicOp<T,NDIM>::DFTCoulombPeriodicOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("CoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
    _cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-4, thresh, L);
    // Initialize potential
    _Vc = FunctionFactory<T,NDIM>(world);
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombPeriodicOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearChargeDensityOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _Vnuc * psi;
    return rfunc;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearPotentialOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _V * psi;
    return rfunc;
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    // If not spin-polarized, multiply by 2.0
    double factor = (_spinpol) ? 1.0 : 2.0;
//    if (_world.rank() == 0) printf("Applying Coulomb operator ...\n\n");
    printf("Applying Coulomb operator ...\n\n");
    funcT rfunc = factor * _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombPeriodicOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    // If not spin-polarized, multiply by 2.0
    double factor = (_spinpol) ? 1.0 : 2.0;
//    if (_world.rank() == 0) printf("Applying Coulomb operator ...\n\n");
    printf("Applying Coulomb operator ...\n\n");
    funcT rfunc = factor * _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  XCFunctionalLDA<T,NDIM>::XCFunctionalLDA(World& world, double coeff, double thresh)
    : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("XCFunctionalLDA");
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> XCFunctionalLDA<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT V_rho = copy(rho);
    V_rho.reconstruct();
    V_rho.unaryop(&dft_xc_lda_V<NDIM>);
    funcT rfunc = V_rho * psi;
    // WSTHORNTON DEBUG
    double pnorm = psi.norm2();
    double rfuncnorm = rfunc.norm2();
    printf("pnorm = %.5f\trfuncnorm = %.5f\n\n", pnorm, rfuncnorm);
    return rfunc;
  }
  //***************************************************************************
}

namespace madness
{
//***************************************************************************
template <int NDIM>
void dft_xc_lda_ene(const Key<NDIM>& key, Tensor<double>& t)
{
  Tensor<double> V = copy(t);
  Tensor<double> enefunc = copy(t);
  xc_rks_generic_lda(t, enefunc, V);
  t(___) = enefunc(___);
}
//***************************************************************************

//  //***************************************************************************
//  template <typename T, int NDIM>
//  DFT<T,NDIM>::DFT(World& world, funcT V, std::vector<funcT> phis,
//      std::vector<double> eigs, double thresh, bool periodic)
//  : _world(world), _V(V), _thresh(thresh)
//  {
//
//    if (world.rank() == 0 && !periodic) printf("DFT constructor (non-peridic) ...\n\n");
//    if (world.rank() == 0 && periodic) printf("DFT constructor (periodic) ...\n\n");
//
//    // Create ops list
//    std::vector<EigSolverOp<T,NDIM>*> ops;
//    // Add nuclear potential to ops list
//    ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
//    if (periodic)
//      ops.push_back(new DFTCoulombPeriodicOp<T,NDIM>(world, 1.0, thresh));
//    else
//      ops.push_back(new DFTCoulombOp<T,NDIM>(world, 1.0, thresh));
//    _xcfunc = new XCFunctionalLDA<T,NDIM>(world, 1.0, thresh);
//    ops.push_back(_xcfunc);
//
//    // Create solver
//    if (periodic)
//    {
//      std::vector<kvecT> kpoints;
//      kvecT gammap(0.0);
//      kpoints.push_back(gammap);
//      _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, kpoints, thresh);
//    }
//    else
//    {
//      _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, thresh);
//    }
//    _solver->addObserver(this);
//
//  }
//  //***************************************************************************
//
  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::DFT(World& world, funcT rhon, std::vector<funcT> phis,
      std::vector<double> eigs, double thresh, bool periodic)
  : _world(world), _rhon(rhon), _thresh(thresh)
  {

    if (world.rank() == 0 && !periodic) printf("DFT constructor (non-peridic) ...\n\n");
    if (world.rank() == 0 && periodic) printf("DFT constructor (periodic) ...\n\n");

    // Create ops list
    std::vector<EigSolverOp<T,NDIM>*> ops;
    // Add nuclear potential to ops list
//    ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
    ops.push_back(new DFTNuclearChargeDensityOp<T,NDIM>(world, rhon, 1.0, thresh, periodic));
    if (periodic)
    {
      ops.push_back(new DFTCoulombPeriodicOp<T,NDIM>(world, 1.0, thresh));
    }
    else
    {
      ops.push_back(new DFTCoulombOp<T,NDIM>(world, 1.0, thresh));
    }
    _xcfunc = new XCFunctionalLDA<T,NDIM>(world, 1.0, thresh);
    ops.push_back(_xcfunc);

    // Create solver
    if (periodic)
    {
      std::vector<kvecT> kpoints;
      kvecT gammap(0.0);
      kpoints.push_back(gammap);
      _solver = new EigSolver<T,NDIM>(world, rhon, phis, eigs, ops, kpoints, thresh);
    }
    else
    {
      _solver = new EigSolver<T,NDIM>(world, rhon, phis, eigs, ops, thresh);
    }
    _solver->addObserver(this);

  }
  //***************************************************************************

//  //***************************************************************************
//  template <typename T, int NDIM>
//  DFT<T,NDIM>::DFT(World& world, funcT V, std::vector<funcT> phis,
//      std::vector<double> eigs, std::vector<kvecT> kpoints, double thresh)
//  : _world(world), _V(V), _thresh(thresh)
//  {
//    // Create ops list
//    std::vector<EigSolverOp<T,NDIM>*> ops;
//    // Add nuclear potential to ops list
//    ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
//    ops.push_back(new DFTCoulombPeriodicOp<T,NDIM>(world, 1.0, thresh));
//    _xcfunc = new XCFunctionalLDA<T,NDIM>(world, 1.0, thresh);
//    ops.push_back(_xcfunc);
//
//    // Create solver
//    _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, kpoints, thresh);
//    _solver->addObserver(this);
//
//  }
//  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::~DFT()
  {
    delete _solver;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::solve(int maxits)
  {
    _solver->solve(maxits);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_ke_sp(funcT psi, bool periodic)
  {
    // Do calculation
    double kenergy = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi = diff(psi, axis);
      kenergy += 0.5 * inner(dpsi, dpsi);
    }
    return kenergy;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_ke_sp(const std::vector<funcT>& phis, bool spinpol,
      bool periodic)
  {
    double tot_ke = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate kinetic energy contribution from psi
      tot_ke += calculate_ke_sp(psi);
    }
    if (!spinpol) tot_ke *= 2.0;
    return tot_ke;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_pe_sp(const World& world, const Function<double,NDIM>& rho,
      const Function<double,NDIM>& rhon, bool spinpol, const double thresh, bool periodic)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM>* op;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      op = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh, L);
    }
    else
    {
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    // Apply Coulomb operator and trace with the density
    funcT Vnuc = apply(*op, rhon);
    double tot_pe = Vnuc.inner(rho);
    if (!spinpol) tot_pe *= 2.0;
    delete op;
    return tot_pe;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_coulomb_energy(const World& world, const Function<double, NDIM>& rho,
      bool spinpol, const double thresh, bool periodic)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM>* op;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      op = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh, L);
    }
    else
    {
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    // Apply Coulomb operator and trace with the density
    funcT Vc = apply(*op, rho);
    double tot_ce = Vc.inner(rho);
    if (!spinpol) tot_ce *= 2.0;
    delete op;
    return tot_ce;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_xc_energy(const Function<double, NDIM>& rho)
  {
    funcT enefunc = copy(rho);
    enefunc.reconstruct();
    enefunc.unaryop(&dft_xc_lda_ene<NDIM>);
    return enefunc.trace();
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs, const Function<double, NDIM>& rho,
      const int& iter, bool periodic)
  {
    if (iter%3 == 0)
    {
      if (world().rank() == 0) printf("Calculating energies ...\n");
      if (world().rank() == 0) printf("Calculating KE ...\n");
      double ke = DFT::calculate_tot_ke_sp(phis, false, periodic);
      if (world().rank() == 0) printf("Calculating PE ...\n");
      double pe = DFT::calculate_tot_pe_sp(_world, rho, _rhon, false, _thresh, periodic);
      if (world().rank() == 0) printf("Calculating CE ...\n");
      double ce = DFT::calculate_tot_coulomb_energy(_world, rho, false, _thresh, periodic);
      if (world().rank() == 0) printf("Calculating EE ...\n");
      double xce = DFT::calculate_tot_xc_energy(rho);
      if (world().rank() == 0) printf("Calculating NE ...\n");
      double ne = 0.0;
      if (world().rank() == 0) printf("Kinetic energy:\t\t\t %.8f\n", ke);
      if (world().rank() == 0) printf("Potential energy:\t\t %.8f\n", pe);
      if (world().rank() == 0) printf("Coulomb energy:\t\t\t %.8f\n", ce);
      if (world().rank() == 0) printf("XC energy:\t\t\t %.8f\n", xce);
      if (world().rank() == 0) printf("Total energy:\t\t\t %.8f\n", ke + pe + ce + xce + ne);
      if (world().rank() == 0) printf("gs ene = %.4f\n", eigs[0]);
      if (world().rank() == 0) printf("1st es ene = %.4f\n", eigs[1]);
      T mtxe = matrix_element(phis[0], phis[0]);
      if (world().rank() == 0) printf("\nKS matrix element:\t\t\t%.8f\n\n", mtxe);
      print_matrix_elements(phis[0], phis[0]);
    }
  }
  //***************************************************************************

  //***************************************************************************
  template class DFT<double, 1>;
  template class DFT<double, 2>;
  template class DFT<double, 3>;


//  template class DFT< std::complex<double> >;
  //***************************************************************************
}

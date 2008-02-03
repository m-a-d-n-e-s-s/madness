/* xc_x_pw91.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = .33333333333333331;

/* :X_PW91subrstart */
/*    Generated: Sun Oct 24 15:07:52 BST 2004 */
/* Subroutine */ int uks_x_pw91__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *rhob1, doublereal *sigmaaa1, doublereal *sigmabb1, 
	doublereal *sigmaab1, doublereal *zk, doublereal *vrhoa, doublereal *
	vrhob, doublereal *vsigmaaa, doublereal *vsigmabb, doublereal *
	vsigmaab, doublereal *v2rhoa2, doublereal *v2rhob2, doublereal *
	v2rhoab, doublereal *v2rhoasigmaaa, doublereal *v2rhoasigmaab, 
	doublereal *v2rhoasigmabb, doublereal *v2rhobsigmabb, doublereal *
	v2rhobsigmaab, doublereal *v2rhobsigmaaa, doublereal *v2sigmaaa2, 
	doublereal *v2sigmaaaab, doublereal *v2sigmaaabb, doublereal *
	v2sigmaab2, doublereal *v2sigmaabbb, doublereal *v2sigmabb2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t6, t7, t8, t10, t11, t12, t13, t14, 
	    t40, t25, t17, t26, t16, t19, t27, t28, t38, t39, t42, t44, t46, 
	    t47, t48, t50, t53, t61, t62, t15, t20, t23, t24, t29, t32, t33, 
	    t43, t51, t57, t66, t67, t76, t78, t21, t22, t31, t34, t35, t41, 
	    t55, t56, t59, t60, t65, t68, t69, t79, t82, t86, t87, t89, t93, 
	    t102, t103, t116, t119, t123, t124, t126, t130, t139, t140, t149, 
	    t151, t170, t172, t36, t49, t63, t70, t85, t90, t91, t107, t109, 
	    t112, t118, t120, t153, t157, t159, t171, t72, t80, t92, t98, t99,
	     t105, t106, t113, t117, t122, t129, t135, t136, t142, t143, t146,
	     t152, t158, t163, t164, t167, t173, t178, t179, t184, t185, t201,
	     t203, rho, t206, t210, t212, t214, t218, t233, t234, t256, t258, 
	    t261, t265, t267, t269, t273, t288, t289, t305, t307, t308, t310, 
	    t345, t347, t348, t350, t382, t386, t388, t400, t413, t417, t419, 
	    t431, rhoa, rhob, sigma, sigmaaa, sigmaab, sigmabb;


/*     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, */
/*     M.R. Pederson, D.J. Singh, C. Fiolhais */
/*     Atoms, molecules, solids and surfaces: */
/*     Applications of the generalized gradient approximation */
/*     for exchange and correlation */
/*     Phys. Rev. B 46 (1992) 6671--6687 */


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
    --v2sigmabb2;
    --v2sigmaabbb;
    --v2sigmaab2;
    --v2sigmaaabb;
    --v2sigmaaaab;
    --v2sigmaaa2;
    --v2rhobsigmaaa;
    --v2rhobsigmaab;
    --v2rhobsigmabb;
    --v2rhoasigmabb;
    --v2rhoasigmaab;
    --v2rhoasigmaaa;
    --v2rhoab;
    --v2rhob2;
    --v2rhoa2;
    --vsigmaab;
    --vsigmabb;
    --vsigmaaa;
    --vrhob;
    --vrhoa;
    --zk;
    --sigmaab1;
    --sigmabb1;
    --sigmaaa1;
    --rhob1;
    --rhoa1;

    /* Function Body */
    if (*ideriv == 0) {
	i__1 = *npt;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = 0., d__2 = rhoa1[i__];
	    rhoa = max(d__1,d__2);
/* Computing MAX */
	    d__1 = 0., d__2 = rhob1[i__];
	    rhob = max(d__1,d__2);
	    rho = rhoa + rhob;
	    if (rho > 1e-20) {
		if (rhoa < 1e-20) {
		    rho = rhob;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmabb;
		    t2 = pow_dd(&rhob, &c_b2);
		    t3 = t2 * rhob;
		    t4 = sqrt(sigmabb);
		    t6 = t4 / t3;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhob;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t17 = exp(sigmabb * -1.645530784602056 * t14);
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    zk[i__] = t3 * -.9305257363491 * (t10 + 1. + (.2743 - t17 
			    * .1508) * .01645530784602056 * sigmabb * t14) / (
			    t10 + 1. + t25 * 1.083108625229223e-6 / t2 / t26 /
			     rhob);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t4 = sqrt(sigmaaa);
		    t6 = t4 / t3;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t17 = exp(sigmaaa * -1.645530784602056 * t14);
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    zk[i__] = t3 * -.9305257363491 * (t10 + 1. + (.2743 - t17 
			    * .1508) * .01645530784602056 * sigmaaa * t14) / (
			    t10 + 1. + t25 * 1.083108625229223e-6 / t2 / t26 /
			     rhoa);
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigmaab = sigmaab1[i__];
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmaaa + sigmabb + sigmaab * 2.;
		    t4 = pow_dd(&rhoa, &c_b2);
		    t5 = t4 * rhoa;
		    t6 = sqrt(sigmaaa);
		    t8 = t6 / t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t10 = log(t8 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t12 = t8 * .02520026100493014 * t10;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t13 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t14 = d__1 * d__1;
		    t16 = 1 / t14 / t13;
		    t19 = exp(sigmaaa * -1.645530784602056 * t16);
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t27 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t13;
		    t28 = d__1 * d__1;
		    t38 = pow_dd(&rhob, &c_b2);
		    t39 = t38 * rhob;
		    t40 = sqrt(sigmabb);
		    t42 = t40 / t39;
/* Computing 2nd power */
		    d__1 = t42;
		    t44 = log(t42 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t46 = t42 * .02520026100493014 * t44;
/* Computing 2nd power */
		    d__1 = rhob;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t38;
		    t48 = d__1 * d__1;
		    t50 = 1 / t48 / t47;
		    t53 = exp(sigmabb * -1.645530784602056 * t50);
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t62 = d__1 * d__1;
		    zk[i__] = t5 * -.9305257363491 * (t12 + 1. + (.2743 - t19 
			    * .1508) * .01645530784602056 * sigmaaa * t16) / (
			    t12 + 1. + t27 * 1.083108625229223e-6 / t4 / t28 /
			     rhoa) - t39 * .9305257363491 * (t46 + 1. + (
			    .2743 - t53 * .1508) * .01645530784602056 * 
			    sigmabb * t50) / (t46 + 1. + t61 * 
			    1.083108625229223e-6 / t38 / t62 / rhob);
		}
/* rhoa,rhob */
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
	    rhoa = max(d__1,d__2);
/* Computing MAX */
	    d__1 = 0., d__2 = rhob1[i__];
	    rhob = max(d__1,d__2);
	    rho = rhoa + rhob;
	    if (rho > 1e-20) {
		if (rhoa < 1e-20) {
		    rho = rhob;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmabb;
		    t2 = pow_dd(&rhob, &c_b2);
		    t3 = t2 * rhob;
		    t4 = sqrt(sigmabb);
		    t5 = 1 / t3;
		    t6 = t4 * t5;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhob;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t15 = sigmabb * t14;
		    t17 = exp(t15 * -1.645530784602056);
		    t19 = .2743 - t17 * .1508;
		    t20 = t19 * sigmabb;
		    t23 = t10 + 1. + t20 * .01645530784602056 * t14;
		    t24 = t3 * t23;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    t29 = 1 / t2 / t26 / rhob;
		    t32 = t10 + 1. + t25 * 1.083108625229223e-6 * t29;
		    t33 = 1 / t32;
		    zk[i__] = t24 * -.9305257363491 * t33;
		    vrhoa[i__] = 0.;
		    t43 = t4 * .03360034800657352 / t2 / t11 * t8;
		    t46 = 1 / t12 / t11 / rhob;
		    t50 = sqrt(t15 * 1.0000117555961 + 1.);
		    t51 = 1 / t50;
		    t53 = sigmabb * .03360054550205309 * t46 * t51;
		    t57 = t25 / t2 / t26 / t11;
/* Computing 2nd power */
		    d__1 = t32;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    vrhob[i__] = t2 * -1.2407009817988 * t23 * t33 - t3 * 
			    .9305257363491 * (-t43 - t53 - t57 * 
			    .01088885204563779 * t17 - t20 * 
			    .04388082092272149 * t46) * t33 + t24 * 
			    .9305257363491 * t67 * (-t43 - t53 - t57 * 
			    5.776579334555855e-6);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t76 = .01260013050246507 / t4 * t5 * t8;
		    t78 = t14 * .01260020456326991 * t51;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t76 + t78 + t29 * 
			    .00408331951711417 * t17 * sigmabb + t19 * 
			    .01645530784602056 * t14) * t33 + t24 * 
			    .9305257363491 * t67 * (t76 + t78 + sigmabb * 
			    2.166217250458446e-6 * t29);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t4 = sqrt(sigmaaa);
		    t5 = 1 / t3;
		    t6 = t4 * t5;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t15 = sigmaaa * t14;
		    t17 = exp(t15 * -1.645530784602056);
		    t19 = .2743 - t17 * .1508;
		    t20 = t19 * sigmaaa;
		    t23 = t10 + 1. + t20 * .01645530784602056 * t14;
		    t24 = t3 * t23;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    t29 = 1 / t2 / t26 / rhoa;
		    t32 = t10 + 1. + t25 * 1.083108625229223e-6 * t29;
		    t33 = 1 / t32;
		    zk[i__] = t24 * -.9305257363491 * t33;
		    t43 = t4 * .03360034800657352 / t2 / t11 * t8;
		    t46 = 1 / t12 / t11 / rhoa;
		    t50 = sqrt(t15 * 1.0000117555961 + 1.);
		    t51 = 1 / t50;
		    t53 = sigmaaa * .03360054550205309 * t46 * t51;
		    t57 = t25 / t2 / t26 / t11;
/* Computing 2nd power */
		    d__1 = t32;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    vrhoa[i__] = t2 * -1.2407009817988 * t23 * t33 - t3 * 
			    .9305257363491 * (-t43 - t53 - t57 * 
			    .01088885204563779 * t17 - t20 * 
			    .04388082092272149 * t46) * t33 + t24 * 
			    .9305257363491 * t67 * (-t43 - t53 - t57 * 
			    5.776579334555855e-6);
		    vrhob[i__] = 0.;
		    t76 = .01260013050246507 / t4 * t5 * t8;
		    t78 = t14 * .01260020456326991 * t51;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t76 + t78 + t29 * 
			    .00408331951711417 * t17 * sigmaaa + t19 * 
			    .01645530784602056 * t14) * t33 + t24 * 
			    .9305257363491 * t67 * (t76 + t78 + sigmaaa * 
			    2.166217250458446e-6 * t29);
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigmaab = sigmaab1[i__];
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmaaa + sigmabb + sigmaab * 2.;
		    t4 = pow_dd(&rhoa, &c_b2);
		    t5 = t4 * rhoa;
		    t6 = sqrt(sigmaaa);
		    t7 = 1 / t5;
		    t8 = t6 * t7;
/* Computing 2nd power */
		    d__1 = t8;
		    t10 = log(t8 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t12 = t8 * .02520026100493014 * t10;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t13 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t14 = d__1 * d__1;
		    t16 = 1 / t14 / t13;
		    t17 = sigmaaa * t16;
		    t19 = exp(t17 * -1.645530784602056);
		    t21 = .2743 - t19 * .1508;
		    t22 = t21 * sigmaaa;
		    t25 = t12 + 1. + t22 * .01645530784602056 * t16;
		    t26 = t5 * t25;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t27 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t13;
		    t28 = d__1 * d__1;
		    t31 = 1 / t4 / t28 / rhoa;
		    t34 = t12 + 1. + t27 * 1.083108625229223e-6 * t31;
		    t35 = 1 / t34;
		    t38 = pow_dd(&rhob, &c_b2);
		    t39 = t38 * rhob;
		    t40 = sqrt(sigmabb);
		    t41 = 1 / t39;
		    t42 = t40 * t41;
/* Computing 2nd power */
		    d__1 = t42;
		    t44 = log(t42 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t46 = t42 * .02520026100493014 * t44;
/* Computing 2nd power */
		    d__1 = rhob;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t38;
		    t48 = d__1 * d__1;
		    t50 = 1 / t48 / t47;
		    t51 = sigmabb * t50;
		    t53 = exp(t51 * -1.645530784602056);
		    t55 = .2743 - t53 * .1508;
		    t56 = t55 * sigmabb;
		    t59 = t46 + 1. + t56 * .01645530784602056 * t50;
		    t60 = t39 * t59;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t62 = d__1 * d__1;
		    t65 = 1 / t38 / t62 / rhob;
		    t68 = t46 + 1. + t61 * 1.083108625229223e-6 * t65;
		    t69 = 1 / t68;
		    zk[i__] = t26 * -.9305257363491 * t35 - t60 * 
			    .9305257363491 * t69;
		    t79 = t6 * .03360034800657352 / t4 / t13 * t10;
		    t82 = 1 / t14 / t13 / rhoa;
		    t86 = sqrt(t17 * 1.0000117555961 + 1.);
		    t87 = 1 / t86;
		    t89 = sigmaaa * .03360054550205309 * t82 * t87;
		    t93 = t27 / t4 / t28 / t13;
/* Computing 2nd power */
		    d__1 = t34;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
		    vrhoa[i__] = t4 * -1.2407009817988 * t25 * t35 - t5 * 
			    .9305257363491 * (-t79 - t89 - t93 * 
			    .01088885204563779 * t19 - t22 * 
			    .04388082092272149 * t82) * t35 + t26 * 
			    .9305257363491 * t103 * (-t79 - t89 - t93 * 
			    5.776579334555855e-6);
		    t116 = t40 * .03360034800657352 / t38 / t47 * t44;
		    t119 = 1 / t48 / t47 / rhob;
		    t123 = sqrt(t51 * 1.0000117555961 + 1.);
		    t124 = 1 / t123;
		    t126 = sigmabb * .03360054550205309 * t119 * t124;
		    t130 = t61 / t38 / t62 / t47;
/* Computing 2nd power */
		    d__1 = t68;
		    t139 = d__1 * d__1;
		    t140 = 1 / t139;
		    vrhob[i__] = t38 * -1.2407009817988 * t59 * t69 - t39 * 
			    .9305257363491 * (-t116 - t126 - t130 * 
			    .01088885204563779 * t53 - t56 * 
			    .04388082092272149 * t119) * t69 + t60 * 
			    .9305257363491 * t140 * (-t116 - t126 - t130 * 
			    5.776579334555855e-6);
		    t149 = .01260013050246507 / t6 * t7 * t10;
		    t151 = t16 * .01260020456326991 * t87;
		    vsigmaaa[i__] = t5 * -.9305257363491 * (t149 + t151 + t31 
			    * .00408331951711417 * t19 * sigmaaa + t21 * 
			    .01645530784602056 * t16) * t35 + t26 * 
			    .9305257363491 * t103 * (t149 + t151 + sigmaaa * 
			    2.166217250458446e-6 * t31);
		    vsigmaab[i__] = 0.;
		    t170 = .01260013050246507 / t40 * t41 * t44;
		    t172 = t50 * .01260020456326991 * t124;
		    vsigmabb[i__] = t39 * -.9305257363491 * (t170 + t172 + 
			    t65 * .00408331951711417 * t53 * sigmabb + t55 * 
			    .01645530784602056 * t50) * t69 + t60 * 
			    .9305257363491 * t140 * (t170 + t172 + sigmabb * 
			    2.166217250458446e-6 * t65);
		}
/* rhoa,rhob */
	    } else {
/* rho */
		zk[i__] = 0.;
		vrhoa[i__] = 0.;
		vrhob[i__] = 0.;
		vsigmaaa[i__] = 0.;
		vsigmaab[i__] = 0.;
		vsigmabb[i__] = 0.;
	    }
/* rho */
	}
    } else if (*ideriv == 2) {
	i__1 = *npt;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = 0., d__2 = rhoa1[i__];
	    rhoa = max(d__1,d__2);
/* Computing MAX */
	    d__1 = 0., d__2 = rhob1[i__];
	    rhob = max(d__1,d__2);
	    rho = rhoa + rhob;
	    if (rho > 1e-20) {
		if (rhoa < 1e-20) {
		    rho = rhob;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmabb;
		    t2 = pow_dd(&rhob, &c_b2);
		    t3 = t2 * rhob;
		    t4 = sqrt(sigmabb);
		    t5 = 1 / t3;
		    t6 = t4 * t5;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhob;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t15 = sigmabb * t14;
		    t17 = exp(t15 * -1.645530784602056);
		    t19 = .2743 - t17 * .1508;
		    t20 = t19 * sigmabb;
		    t23 = t10 + 1. + t20 * .01645530784602056 * t14;
		    t24 = t3 * t23;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    t29 = 1 / t2 / t26 / rhob;
		    t32 = t10 + 1. + t25 * 1.083108625229223e-6 * t29;
		    t33 = 1 / t32;
		    zk[i__] = t24 * -.9305257363491 * t33;
		    vrhoa[i__] = 0.;
		    t36 = t2 * t23;
		    t43 = t4 * .03360034800657352 / t2 / t11 * t8;
		    t44 = t11 * rhob;
		    t46 = 1 / t12 / t44;
		    t49 = t15 * 1.0000117555961 + 1.;
		    t50 = sqrt(t49);
		    t51 = 1 / t50;
		    t53 = sigmabb * .03360054550205309 * t46 * t51;
		    t57 = t25 / t2 / t26 / t11;
		    t62 = -t43 - t53 - t57 * .01088885204563779 * t17 - t20 * 
			    .04388082092272149 * t46;
		    t63 = t3 * t62;
/* Computing 2nd power */
		    d__1 = t32;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t69 = -t43 - t53 - t57 * 5.776579334555855e-6;
		    t70 = t67 * t69;
		    vrhob[i__] = t36 * -1.2407009817988 * t33 - t63 * 
			    .9305257363491 * t33 + t24 * .9305257363491 * t70;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t76 = .01260013050246507 / t4 * t5 * t8;
		    t78 = t14 * .01260020456326991 * t51;
		    t79 = t29 * t17;
		    t85 = t3 * (t76 + t78 + t79 * .00408331951711417 * 
			    sigmabb + t19 * .01645530784602056 * t14);
		    t90 = t76 + t78 + sigmabb * 2.166217250458446e-6 * t29;
		    t91 = t67 * t90;
		    vsigmabb[i__] = t85 * -.9305257363491 * t33 + t24 * 
			    .9305257363491 * t91;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t107 = t4 * .07840081201533821 / t2 / t44 * t8;
		    t109 = 1 / t12 / t26;
		    t112 = sigmabb * .1680027275102654 * t109 * t51;
		    t116 = t25 / t2 / t26 / t44;
		    t118 = 1 / t50 / t49;
		    t120 = t116 * .04480125399532632 * t118;
/* Computing 2nd power */
		    d__1 = t26;
		    t124 = d__1 * d__1;
		    t139 = 1 / t66 / t32;
/* Computing 2nd power */
		    d__1 = t69;
		    t140 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t12 * t23 * t33 - t2 *
			     2.4814019635976 * t62 * t33 + t36 * 
			    2.4814019635976 * t70 - t3 * .9305257363491 * (
			    t107 + t112 - t120 + t116 * .09799966841074009 * 
			    t17 - t25 * .04778117666686413 * sigmabb / t124 / 
			    t11 * t17 + t20 * .1608963433833121 * t109) * t33 
			    + t63 * 1.8610514726982 * t70 - t24 * 
			    1.8610514726982 * t139 * t140 + t24 * 
			    .9305257363491 * t67 * (t107 + t112 - t120 + t116 
			    * 3.658500245218708e-5);
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t153 = .006300065251232535 / t4 / sigmabb * t5 * t8;
		    t157 = .006300102281634954 / sigmabb * t14 * t51;
		    t159 = t29 * .006300176343092764 * t118;
/* Computing 2nd power */
		    d__1 = t90;
		    t171 = d__1 * d__1;
		    v2sigmabb2[i__] = t3 * -.9305257363491 * (-t153 + t157 - 
			    t159 - .006719227968777768 / t124 * t17 * sigmabb 
			    + t79 * .008166639034228341) * t33 + t85 * 
			    1.8610514726982 * t91 - t24 * 1.8610514726982 * 
			    t139 * t171 + t24 * .9305257363491 * t67 * (-t153 
			    + t157 - t159 + t29 * 2.166217250458446e-6);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t4 = sqrt(sigmaaa);
		    t5 = 1 / t3;
		    t6 = t4 * t5;
/* Computing 2nd power */
		    d__1 = t6;
		    t8 = log(t6 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t10 = t6 * .02520026100493014 * t8;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t11 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t12 = d__1 * d__1;
		    t14 = 1 / t12 / t11;
		    t15 = sigmaaa * t14;
		    t17 = exp(t15 * -1.645530784602056);
		    t19 = .2743 - t17 * .1508;
		    t20 = t19 * sigmaaa;
		    t23 = t10 + 1. + t20 * .01645530784602056 * t14;
		    t24 = t3 * t23;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t11;
		    t26 = d__1 * d__1;
		    t29 = 1 / t2 / t26 / rhoa;
		    t32 = t10 + 1. + t25 * 1.083108625229223e-6 * t29;
		    t33 = 1 / t32;
		    zk[i__] = t24 * -.9305257363491 * t33;
		    t36 = t2 * t23;
		    t43 = t4 * .03360034800657352 / t2 / t11 * t8;
		    t44 = t11 * rhoa;
		    t46 = 1 / t12 / t44;
		    t49 = t15 * 1.0000117555961 + 1.;
		    t50 = sqrt(t49);
		    t51 = 1 / t50;
		    t53 = sigmaaa * .03360054550205309 * t46 * t51;
		    t57 = t25 / t2 / t26 / t11;
		    t62 = -t43 - t53 - t57 * .01088885204563779 * t17 - t20 * 
			    .04388082092272149 * t46;
		    t63 = t3 * t62;
/* Computing 2nd power */
		    d__1 = t32;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t69 = -t43 - t53 - t57 * 5.776579334555855e-6;
		    t70 = t67 * t69;
		    vrhoa[i__] = t36 * -1.2407009817988 * t33 - t63 * 
			    .9305257363491 * t33 + t24 * .9305257363491 * t70;
		    vrhob[i__] = 0.;
		    t76 = .01260013050246507 / t4 * t5 * t8;
		    t78 = t14 * .01260020456326991 * t51;
		    t79 = t29 * t17;
		    t85 = t3 * (t76 + t78 + t79 * .00408331951711417 * 
			    sigmaaa + t19 * .01645530784602056 * t14);
		    t90 = t76 + t78 + sigmaaa * 2.166217250458446e-6 * t29;
		    t91 = t67 * t90;
		    vsigmaaa[i__] = t85 * -.9305257363491 * t33 + t24 * 
			    .9305257363491 * t91;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t107 = t4 * .07840081201533821 / t2 / t44 * t8;
		    t109 = 1 / t12 / t26;
		    t112 = sigmaaa * .1680027275102654 * t109 * t51;
		    t116 = t25 / t2 / t26 / t44;
		    t118 = 1 / t50 / t49;
		    t120 = t116 * .04480125399532632 * t118;
/* Computing 2nd power */
		    d__1 = t26;
		    t124 = d__1 * d__1;
		    t139 = 1 / t66 / t32;
/* Computing 2nd power */
		    d__1 = t69;
		    t140 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t12 * t23 * t33 - t2 *
			     2.4814019635976 * t62 * t33 + t36 * 
			    2.4814019635976 * t70 - t3 * .9305257363491 * (
			    t107 + t112 - t120 + t116 * .09799966841074009 * 
			    t17 - t25 * .04778117666686413 * sigmaaa / t124 / 
			    t11 * t17 + t20 * .1608963433833121 * t109) * t33 
			    + t63 * 1.8610514726982 * t70 - t24 * 
			    1.8610514726982 * t139 * t140 + t24 * 
			    .9305257363491 * t67 * (t107 + t112 - t120 + t116 
			    * 3.658500245218708e-5);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t153 = .006300065251232535 / t4 / sigmaaa * t5 * t8;
		    t157 = .006300102281634954 / sigmaaa * t14 * t51;
		    t159 = t29 * .006300176343092764 * t118;
/* Computing 2nd power */
		    d__1 = t90;
		    t171 = d__1 * d__1;
		    v2sigmaaa2[i__] = t3 * -.9305257363491 * (-t153 + t157 - 
			    t159 - .006719227968777768 / t124 * t17 * sigmaaa 
			    + t79 * .008166639034228341) * t33 + t85 * 
			    1.8610514726982 * t91 - t24 * 1.8610514726982 * 
			    t139 * t171 + t24 * .9305257363491 * t67 * (-t153 
			    + t157 - t159 + t29 * 2.166217250458446e-6);
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    v2sigmabb2[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigmaab = sigmaab1[i__];
/* Computing MAX */
		    d__1 = 0., d__2 = sigmabb1[i__];
		    sigmabb = max(d__1,d__2);
		    sigma = sigmaaa + sigmabb + sigmaab * 2.;
		    t4 = pow_dd(&rhoa, &c_b2);
		    t5 = t4 * rhoa;
		    t6 = sqrt(sigmaaa);
		    t7 = 1 / t5;
		    t8 = t6 * t7;
/* Computing 2nd power */
		    d__1 = t8;
		    t10 = log(t8 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t12 = t8 * .02520026100493014 * t10;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t13 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t14 = d__1 * d__1;
		    t16 = 1 / t14 / t13;
		    t17 = sigmaaa * t16;
		    t19 = exp(t17 * -1.645530784602056);
		    t21 = .2743 - t19 * .1508;
		    t22 = t21 * sigmaaa;
		    t25 = t12 + 1. + t22 * .01645530784602056 * t16;
		    t26 = t5 * t25;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t27 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t13;
		    t28 = d__1 * d__1;
		    t31 = 1 / t4 / t28 / rhoa;
		    t34 = t12 + 1. + t27 * 1.083108625229223e-6 * t31;
		    t35 = 1 / t34;
		    t38 = pow_dd(&rhob, &c_b2);
		    t39 = t38 * rhob;
		    t40 = sqrt(sigmabb);
		    t41 = 1 / t39;
		    t42 = t40 * t41;
/* Computing 2nd power */
		    d__1 = t42;
		    t44 = log(t42 * 1.000005877780776 + sqrt(d__1 * d__1 * 
			    1.0000117555961 + 1));
		    t46 = t42 * .02520026100493014 * t44;
/* Computing 2nd power */
		    d__1 = rhob;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t38;
		    t48 = d__1 * d__1;
		    t50 = 1 / t48 / t47;
		    t51 = sigmabb * t50;
		    t53 = exp(t51 * -1.645530784602056);
		    t55 = .2743 - t53 * .1508;
		    t56 = t55 * sigmabb;
		    t59 = t46 + 1. + t56 * .01645530784602056 * t50;
		    t60 = t39 * t59;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t62 = d__1 * d__1;
		    t65 = 1 / t38 / t62 / rhob;
		    t68 = t46 + 1. + t61 * 1.083108625229223e-6 * t65;
		    t69 = 1 / t68;
		    zk[i__] = t26 * -.9305257363491 * t35 - t60 * 
			    .9305257363491 * t69;
		    t72 = t4 * t25;
		    t76 = 1 / t4 / t13;
		    t79 = t6 * .03360034800657352 * t76 * t10;
		    t80 = t13 * rhoa;
		    t82 = 1 / t14 / t80;
		    t85 = t17 * 1.0000117555961 + 1.;
		    t86 = sqrt(t85);
		    t87 = 1 / t86;
		    t89 = sigmaaa * .03360054550205309 * t82 * t87;
		    t92 = 1 / t4 / t28 / t13;
		    t93 = t27 * t92;
		    t98 = -t79 - t89 - t93 * .01088885204563779 * t19 - t22 * 
			    .04388082092272149 * t82;
		    t99 = t5 * t98;
/* Computing 2nd power */
		    d__1 = t34;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
		    t105 = -t79 - t89 - t93 * 5.776579334555855e-6;
		    t106 = t103 * t105;
		    vrhoa[i__] = t72 * -1.2407009817988 * t35 - t99 * 
			    .9305257363491 * t35 + t26 * .9305257363491 * 
			    t106;
		    t109 = t38 * t59;
		    t113 = 1 / t38 / t47;
		    t116 = t40 * .03360034800657352 * t113 * t44;
		    t117 = t47 * rhob;
		    t119 = 1 / t48 / t117;
		    t122 = t51 * 1.0000117555961 + 1.;
		    t123 = sqrt(t122);
		    t124 = 1 / t123;
		    t126 = sigmabb * .03360054550205309 * t119 * t124;
		    t129 = 1 / t38 / t62 / t47;
		    t130 = t61 * t129;
		    t135 = -t116 - t126 - t130 * .01088885204563779 * t53 - 
			    t56 * .04388082092272149 * t119;
		    t136 = t39 * t135;
/* Computing 2nd power */
		    d__1 = t68;
		    t139 = d__1 * d__1;
		    t140 = 1 / t139;
		    t142 = -t116 - t126 - t130 * 5.776579334555855e-6;
		    t143 = t140 * t142;
		    vrhob[i__] = t109 * -1.2407009817988 * t69 - t136 * 
			    .9305257363491 * t69 + t60 * .9305257363491 * 
			    t143;
		    t146 = 1 / t6;
		    t149 = t146 * .01260013050246507 * t7 * t10;
		    t151 = t16 * .01260020456326991 * t87;
		    t152 = t31 * t19;
		    t157 = t149 + t151 + t152 * .00408331951711417 * sigmaaa 
			    + t21 * .01645530784602056 * t16;
		    t158 = t5 * t157;
		    t163 = t149 + t151 + sigmaaa * 2.166217250458446e-6 * t31;
		    t164 = t103 * t163;
		    vsigmaaa[i__] = t158 * -.9305257363491 * t35 + t26 * 
			    .9305257363491 * t164;
		    vsigmaab[i__] = 0.;
		    t167 = 1 / t40;
		    t170 = t167 * .01260013050246507 * t41 * t44;
		    t172 = t50 * .01260020456326991 * t124;
		    t173 = t65 * t53;
		    t178 = t170 + t172 + t173 * .00408331951711417 * sigmabb 
			    + t55 * .01645530784602056 * t50;
		    t179 = t39 * t178;
		    t184 = t170 + t172 + sigmabb * 2.166217250458446e-6 * t65;
		    t185 = t140 * t184;
		    vsigmabb[i__] = t179 * -.9305257363491 * t69 + t60 * 
			    .9305257363491 * t185;
		    t201 = t6 * .07840081201533821 / t4 / t80 * t10;
		    t203 = 1 / t14 / t28;
		    t206 = sigmaaa * .1680027275102654 * t203 * t87;
		    t210 = t27 / t4 / t28 / t80;
		    t212 = 1 / t86 / t85;
		    t214 = t210 * .04480125399532632 * t212;
/* Computing 2nd power */
		    d__1 = t28;
		    t218 = d__1 * d__1;
		    t233 = 1 / t102 / t34;
/* Computing 2nd power */
		    d__1 = t105;
		    t234 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t14 * t25 * t35 - t4 *
			     2.4814019635976 * t98 * t35 + t72 * 
			    2.4814019635976 * t106 - t5 * .9305257363491 * (
			    t201 + t206 - t214 + t210 * .09799966841074009 * 
			    t19 - t27 * .04778117666686413 * sigmaaa / t218 / 
			    t13 * t19 + t22 * .1608963433833121 * t203) * t35 
			    + t99 * 1.8610514726982 * t106 - t26 * 
			    1.8610514726982 * t233 * t234 + t26 * 
			    .9305257363491 * t103 * (t201 + t206 - t214 + 
			    t210 * 3.658500245218708e-5);
		    t256 = t40 * .07840081201533821 / t38 / t117 * t44;
		    t258 = 1 / t48 / t62;
		    t261 = sigmabb * .1680027275102654 * t258 * t124;
		    t265 = t61 / t38 / t62 / t117;
		    t267 = 1 / t123 / t122;
		    t269 = t265 * .04480125399532632 * t267;
/* Computing 2nd power */
		    d__1 = t62;
		    t273 = d__1 * d__1;
		    t288 = 1 / t139 / t68;
/* Computing 2nd power */
		    d__1 = t142;
		    t289 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t48 * t59 * t69 - t38 
			    * 2.4814019635976 * t135 * t69 + t109 * 
			    2.4814019635976 * t143 - t39 * .9305257363491 * (
			    t256 + t261 - t269 + t265 * .09799966841074009 * 
			    t53 - t61 * .04778117666686413 * sigmabb / t273 / 
			    t47 * t53 + t56 * .1608963433833121 * t258) * t69 
			    + t136 * 1.8610514726982 * t143 - t60 * 
			    1.8610514726982 * t288 * t289 + t60 * 
			    .9305257363491 * t140 * (t256 + t261 - t269 + 
			    t265 * 3.658500245218708e-5);
		    v2rhoab[i__] = 0.;
		    t305 = t146 * .01680017400328676 * t76 * t10;
		    t307 = t82 * .05040081825307963 * t87;
		    t308 = sigmaaa * t92;
		    t310 = t308 * .01680047024824737 * t212;
		    v2rhoasigmaaa[i__] = t4 * -1.2407009817988 * t157 * t35 + 
			    t72 * 1.2407009817988 * t164 - t5 * 
			    .9305257363491 * (-t305 - t307 + t310 - t92 * 
			    .03266655613691336 * t19 * sigmaaa + t27 * 
			    .01791794125007405 / t218 / rhoa * t19 - t21 * 
			    .04388082092272149 * t82) * t35 + t99 * 
			    .9305257363491 * t164 + t158 * .9305257363491 * 
			    t106 - t26 * 1.8610514726982 * t233 * t105 * t163 
			    + t26 * .9305257363491 * t103 * (-t305 - t307 + 
			    t310 - t308 * 1.155315866911171e-5);
		    v2rhoasigmaab[i__] = 0.;
		    v2rhoasigmabb[i__] = 0.;
		    v2rhobsigmaaa[i__] = 0.;
		    v2rhobsigmaab[i__] = 0.;
		    t345 = t167 * .01680017400328676 * t113 * t44;
		    t347 = t119 * .05040081825307963 * t124;
		    t348 = sigmabb * t129;
		    t350 = t348 * .01680047024824737 * t267;
		    v2rhobsigmabb[i__] = t38 * -1.2407009817988 * t178 * t69 
			    + t109 * 1.2407009817988 * t185 - t39 * 
			    .9305257363491 * (-t345 - t347 + t350 - t129 * 
			    .03266655613691336 * t53 * sigmabb + t61 * 
			    .01791794125007405 / t273 / rhob * t53 - t55 * 
			    .04388082092272149 * t119) * t69 + t136 * 
			    .9305257363491 * t185 + t179 * .9305257363491 * 
			    t143 - t60 * 1.8610514726982 * t288 * t142 * t184 
			    + t60 * .9305257363491 * t140 * (-t345 - t347 + 
			    t350 - t348 * 1.155315866911171e-5);
		    t382 = .006300065251232535 / t6 / sigmaaa * t7 * t10;
		    t386 = .006300102281634954 / sigmaaa * t16 * t87;
		    t388 = t31 * .006300176343092764 * t212;
/* Computing 2nd power */
		    d__1 = t163;
		    t400 = d__1 * d__1;
		    v2sigmaaa2[i__] = t5 * -.9305257363491 * (-t382 + t386 - 
			    t388 - .006719227968777768 / t218 * t19 * sigmaaa 
			    + t152 * .008166639034228341) * t35 + t158 * 
			    1.8610514726982 * t164 - t26 * 1.8610514726982 * 
			    t233 * t400 + t26 * .9305257363491 * t103 * (
			    -t382 + t386 - t388 + t31 * 2.166217250458446e-6);
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t413 = .006300065251232535 / t40 / sigmabb * t41 * t44;
		    t417 = .006300102281634954 / sigmabb * t50 * t124;
		    t419 = t65 * .006300176343092764 * t267;
/* Computing 2nd power */
		    d__1 = t184;
		    t431 = d__1 * d__1;
		    v2sigmabb2[i__] = t39 * -.9305257363491 * (-t413 + t417 - 
			    t419 - .006719227968777768 / t273 * t53 * sigmabb 
			    + t173 * .008166639034228341) * t69 + t179 * 
			    1.8610514726982 * t185 - t60 * 1.8610514726982 * 
			    t288 * t431 + t60 * .9305257363491 * t140 * (
			    -t413 + t417 - t419 + t65 * 2.166217250458446e-6);
		}
/* rhoa,rhob */
	    } else {
/* rho */
		zk[i__] = 0.;
		vrhoa[i__] = 0.;
		vrhob[i__] = 0.;
		v2rhoa2[i__] = 0.;
		v2rhob2[i__] = 0.;
		v2rhoab[i__] = 0.;
		vsigmaaa[i__] = 0.;
		vsigmaab[i__] = 0.;
		vsigmabb[i__] = 0.;
		v2rhoasigmaaa[i__] = 0.;
		v2rhoasigmaab[i__] = 0.;
		v2rhoasigmabb[i__] = 0.;
		v2rhobsigmaaa[i__] = 0.;
		v2rhobsigmaab[i__] = 0.;
		v2rhobsigmabb[i__] = 0.;
		v2sigmaaa2[i__] = 0.;
		v2sigmaab2[i__] = 0.;
		v2sigmabb2[i__] = 0.;
		v2sigmaaaab[i__] = 0.;
		v2sigmaaabb[i__] = 0.;
		v2sigmaabbb[i__] = 0.;
	    }
/* rho */
	}
    }
/* ideriv */
    return 0;
} /* uks_x_pw91__ */

/* Subroutine */ int rks_x_pw91__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t6, t8, t10, t11, t12, t20, t14, t15, 
	    t25, t17, t26, t19, t23, t24, t29, t32, t33, t43, t46, t50, t51, 
	    t53, t57, t66, t67, t76, t78, t36, t40, t44, t49, t56, t62, t63, 
	    t69, t70, t73, t79, t84, t85, t90, t91, t120, t112, t140, t211, 
	    t124, t107, t116, t109, t118, t156, t139, t158, t159, t161, t193, 
	    t197, t199, rho, sigma;


/*     J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, */
/*     M.R. Pederson, D.J. Singh, C. Fiolhais */
/*     Atoms, molecules, solids and surfaces: */
/*     Applications of the generalized gradient approximation */
/*     for exchange and correlation */
/*     Phys. Rev. B 46 (1992) 6671--6687 */


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
/* Computing MAX */
		d__1 = 0., d__2 = sigmaaa1[i__];
		sigma = max(d__1,d__2);
		t2 = pow_dd(&rho, &c_b2);
		t3 = t2 * rho;
		t4 = sqrt(sigma);
		t6 = t4 / t3;
/* Computing 2nd power */
		d__1 = t6;
		t8 = log(t6 * 1.259928455434599 + sqrt(d__1 * d__1 * 
			1.587419712813814 + 1));
		t10 = t6 * .03175033930295641 * t8;
/* Computing 2nd power */
		d__1 = rho;
		t11 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t12 = d__1 * d__1;
		t14 = 1 / t12 / t11;
		t17 = exp(sigma * -2.61211729852336 * t14);
/* Computing 2nd power */
		d__1 = sigma;
		t25 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t11;
		t26 = d__1 * d__1;
		zk[i__] = t3 * -.7385587663820224 * (t10 + 1. + (.2743 - t17 *
			 .1508) * .0261211729852336 * sigma * t14) / (t10 + 
			1. + t25 * 2.72926271249799e-6 / t2 / t26 / rho);
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
/* Computing MAX */
		d__1 = 0., d__2 = sigmaaa1[i__];
		sigma = max(d__1,d__2);
		t2 = pow_dd(&rho, &c_b2);
		t3 = t2 * rho;
		t4 = sqrt(sigma);
		t5 = 1 / t3;
		t6 = t4 * t5;
/* Computing 2nd power */
		d__1 = t6;
		t8 = log(t6 * 1.259928455434599 + sqrt(d__1 * d__1 * 
			1.587419712813814 + 1));
		t10 = t6 * .03175033930295641 * t8;
/* Computing 2nd power */
		d__1 = rho;
		t11 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t12 = d__1 * d__1;
		t14 = 1 / t12 / t11;
		t15 = sigma * t14;
		t17 = exp(t15 * -2.61211729852336);
		t19 = .2743 - t17 * .1508;
		t20 = t19 * sigma;
		t23 = t10 + 1. + t20 * .0261211729852336 * t14;
		t24 = t3 * t23;
/* Computing 2nd power */
		d__1 = sigma;
		t25 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t11;
		t26 = d__1 * d__1;
		t29 = 1 / t2 / t26 / rho;
		t32 = t10 + 1. + t25 * 2.72926271249799e-6 * t29;
		t33 = 1 / t32;
		zk[i__] = t24 * -.7385587663820224 * t33;
		t43 = t4 * .08466757147455043 / t2 / t11 * t8;
		t46 = 1 / t12 / t11 / rho;
		t50 = sqrt(t15 * 1.587419712813815 + 1.);
		t51 = 1 / t50;
		t53 = sigma * .1066750825533289 * t46 * t51;
		t57 = t25 / t2 / t26 / t11;
/* Computing 2nd power */
		d__1 = t32;
		t66 = d__1 * d__1;
		t67 = 1 / t66;
		vrhoa[i__] = t2 * -.9847450218426965 * t23 * t33 - t3 * 
			.3692793831910112 * (-t43 - t53 - t57 * 
			.05487637560595959 * t17 - t20 * .1393129225879125 * 
			t46) * t33 + t24 * .3692793831910112 * t67 * (-t43 - 
			t53 - t57 * 2.911213559997856e-5);
		t76 = .06350067860591282 / t4 * t5 * t8;
		t78 = t14 * .08000631191499664 * t51;
		vsigmaaa[i__] = t3 * -.7385587663820224 * (t76 + t78 + t29 * 
			.0411572817044697 * t17 * sigma + t19 * 
			.1044846919409344 * t14) * t33 + t24 * 
			.7385587663820224 * t67 * (t76 + t78 + sigma * 
			2.183410169998392e-5 * t29);
	    } else {
/* rho */
		zk[i__] = 0.;
		vrhoa[i__] = 0.;
		vsigmaaa[i__] = 0.;
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
/* Computing MAX */
		d__1 = 0., d__2 = sigmaaa1[i__];
		sigma = max(d__1,d__2);
		t2 = pow_dd(&rho, &c_b2);
		t3 = t2 * rho;
		t4 = sqrt(sigma);
		t5 = 1 / t3;
		t6 = t4 * t5;
/* Computing 2nd power */
		d__1 = t6;
		t8 = log(t6 * 1.259928455434599 + sqrt(d__1 * d__1 * 
			1.587419712813814 + 1));
		t10 = t6 * .03175033930295641 * t8;
/* Computing 2nd power */
		d__1 = rho;
		t11 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t12 = d__1 * d__1;
		t14 = 1 / t12 / t11;
		t15 = sigma * t14;
		t17 = exp(t15 * -2.61211729852336);
		t19 = .2743 - t17 * .1508;
		t20 = t19 * sigma;
		t23 = t10 + 1. + t20 * .0261211729852336 * t14;
		t24 = t3 * t23;
/* Computing 2nd power */
		d__1 = sigma;
		t25 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t11;
		t26 = d__1 * d__1;
		t29 = 1 / t2 / t26 / rho;
		t32 = t10 + 1. + t25 * 2.72926271249799e-6 * t29;
		t33 = 1 / t32;
		zk[i__] = t24 * -.7385587663820224 * t33;
		t36 = t2 * t23;
		t40 = 1 / t2 / t11;
		t43 = t4 * .08466757147455043 * t40 * t8;
		t44 = t11 * rho;
		t46 = 1 / t12 / t44;
		t49 = t15 * 1.587419712813815 + 1.;
		t50 = sqrt(t49);
		t51 = 1 / t50;
		t53 = sigma * .1066750825533289 * t46 * t51;
		t56 = 1 / t2 / t26 / t11;
		t57 = t25 * t56;
		t62 = -t43 - t53 - t57 * .05487637560595959 * t17 - t20 * 
			.1393129225879125 * t46;
		t63 = t3 * t62;
/* Computing 2nd power */
		d__1 = t32;
		t66 = d__1 * d__1;
		t67 = 1 / t66;
		t69 = -t43 - t53 - t57 * 2.911213559997856e-5;
		t70 = t67 * t69;
		vrhoa[i__] = t36 * -.9847450218426965 * t33 - t63 * 
			.3692793831910112 * t33 + t24 * .3692793831910112 * 
			t70;
		t73 = 1 / t4;
		t76 = t73 * .06350067860591282 * t5 * t8;
		t78 = t14 * .08000631191499664 * t51;
		t79 = t29 * t17;
		t84 = t76 + t78 + t79 * .0411572817044697 * sigma + t19 * 
			.1044846919409344 * t14;
		t85 = t3 * t84;
		t90 = t76 + t78 + sigma * 2.183410169998392e-5 * t29;
		t91 = t67 * t90;
		vsigmaaa[i__] = t85 * -.7385587663820224 * t33 + t24 * 
			.7385587663820224 * t91;
		t107 = t4 * .395115333547902 / t2 / t44 * t8;
		t109 = 1 / t12 / t26;
		t112 = sigma * 1.066750825533289 * t109 * t51;
		t116 = t25 / t2 / t26 / t44;
		t118 = 1 / t50 / t49;
		t120 = t116 * .4515683437631874 * t118;
/* Computing 2nd power */
		d__1 = t26;
		t124 = d__1 * d__1;
		t139 = 1 / t66 / t32;
/* Computing 2nd power */
		d__1 = t69;
		t140 = d__1 * d__1;
		v2rhoa2[i__] = -.6564966812284644 / t12 * t23 * t33 - t2 * 
			1.969490043685393 * t62 * t33 + t36 * 
			1.969490043685393 * t70 - t3 * .3692793831910112 * (
			t107 + t112 - t120 + t116 * .9877747609072727 * t17 - 
			t25 * .764498826669826 * sigma / t124 / t11 * t17 + 
			t20 * 1.021628098978025 * t109) * t33 + t63 * 
			.7385587663820224 * t70 - t24 * .7385587663820224 * 
			t139 * t140 + t24 * .3692793831910112 * t67 * (t107 + 
			t112 - t120 + t116 * 3.687537175997285e-4);
		t156 = t73 * .1693351429491009 * t40 * t8;
		t158 = t46 * .6400504953199731 * t51;
		t159 = sigma * t56;
		t161 = t159 * .3386762578223905 * t118;
		v2rhoasigmaaa[i__] = t2 * -.9847450218426965 * t84 * t33 + 
			t36 * .9847450218426965 * t91 - t3 * 
			.3692793831910112 * (-t156 - t158 + t161 - t56 * 
			.6585165072715151 * t17 * sigma + t25 * 
			.5733741200023695 / t124 / rho * t17 - t19 * 
			.5572516903516501 * t46) * t33 + t63 * 
			.3692793831910112 * t91 + t85 * .3692793831910112 * 
			t70 - t24 * .7385587663820224 * t139 * t69 * t90 + 
			t24 * .3692793831910112 * t67 * (-t156 - t158 + t161 
			- t159 * 2.328970847998285e-4);
		t193 = .1270013572118256 / t4 / sigma * t5 * t8;
		t197 = .1600126238299933 / sigma * t14 * t51;
		t199 = t29 * .2540071933667929 * t118;
/* Computing 2nd power */
		d__1 = t90;
		t211 = d__1 * d__1;
		v2sigmaaa2[i__] = t3 * -.7385587663820224 * (-t193 + t197 - 
			t199 - .4300305900017772 / t124 * t17 * sigma + t79 * 
			.3292582536357576) * t33 + t85 * 1.477117532764045 * 
			t91 - t24 * 1.477117532764045 * t139 * t211 + t24 * 
			.7385587663820224 * t67 * (-t193 + t197 - t199 + t29 *
			 8.733640679993569e-5);
	    } else {
/* rho */
		zk[i__] = 0.;
		vrhoa[i__] = 0.;
		v2rhoa2[i__] = 0.;
		vsigmaaa[i__] = 0.;
		v2rhoasigmaaa[i__] = 0.;
		v2sigmaaa2[i__] = 0.;
	    }
/* rho */
	}
    }
/* ideriv */
    return 0;
} /* rks_x_pw91__ */


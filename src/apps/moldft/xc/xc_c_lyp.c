/* xc_c_lyp.f -- translated by f2c (version 20050501).
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

/* :C_LYPsubrstart */
/*    Generated: Wed Mar 10 10:16:01 GMT 2004 */
/* Subroutine */ int uks_c_lyp__(integer *ideriv, integer *npt, doublereal *
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
    double pow_dd(doublereal *, doublereal *), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t4, t5, t7, t8, t9, t10, t11, t30, t21, t14, t15, t23, 
	    t17, t24, t19, t25, t28, t29, t34, t56, t16, t22, t36, t41, t45, 
	    t49, t52, t57, t60, t63, t73, t84, t85, t91, t92, t93, t95, t96, 
	    t98, t18, t27, t32, t37, t43, t51, t67, t71, t76, t80, t86, t88, 
	    t97, t99, t100, t102, t111, t113, t103, t106, t115, t170, t154, 
	    t119, t171, t108, t185, t186, t121, t132, t140, t143, t151, t158, 
	    t162, t166, t172, t175, t177, t183, t188, t191, t193, t198, t200, 
	    t203, t206, t209, t212, t215, t222, t237, t239, t243, t245, t247, 
	    t248, t271, t285, t286, t287, t290, t293, t298, t300, t302, t305, 
	    t307, t308, t310, t313, t316, t320, t323, t325, t335, t337, t340, 
	    t343, t346, t357, rho, t377, t380, t383, t390, t393, t396, t397, 
	    t398, t399, t407, t411, t413, t414, t416, t417, t419, t420, t422, 
	    t423, t428, t429, t431, t432, t433, t434, t443, t446, t449, t460, 
	    t469, t470, rhoa, rhob, sigma, sigmaaa, sigmaab, sigmabb;


/*     C. Lee, W. Yang, and R.G. Parr */
/*     Development of the Colle-Salvetti correlation-energy formula into */
/*     a functional of the electron density */
/*     Phys. Rev. B37 (1988) 785-789 */


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
		    zk[i__] = 0.;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    zk[i__] = 0.;
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
		    t4 = pow_dd(&rho, &c_b2);
		    t5 = 1 / t4;
		    t8 = 1 / (t5 * .349 + 1.);
		    t10 = 1 / rho;
		    t11 = rhob * t10;
		    t14 = t5 * .2533;
		    t15 = exp(-t14);
/* Computing 2nd power */
		    d__1 = rho;
		    t17 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t19 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t23 = d__1 * d__1;
		    t24 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t34 = t5 * t8;
		    t56 = t17 * .6666666666666667;
		    zk[i__] = t8 * -.19672 * rhoa * t11 - t15 * .00649176 * 
			    t8 / t19 / t17 / rho * (rhoa * rhob * (t25 * 
			    36.46239897876478 * t23 + t30 * 36.46239897876478 
			    * t28 + (2.611111111111111 - t5 * 
			    .09850555555555556 - t34 * .1357222222222222) * 
			    sigma - (2.5 - t5 * .01407222222222222 - t34 * 
			    .01938888888888889) * 1. * (sigmaaa + sigmabb) - (
			    t14 + t34 * .349 - 11.) * .1111111111111111 * (
			    rhoa * t10 * sigmaaa + t11 * sigmabb)) - t17 * 
			    .6666666666666667 * sigma + (t56 - t23 * 1.) * 
			    sigmabb + (t56 - t28 * 1.) * sigmaaa);
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
		    zk[i__] = 0.;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    zk[i__] = 0.;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = 0.;
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
		    t4 = pow_dd(&rho, &c_b2);
		    t5 = 1 / t4;
		    t7 = t5 * .349 + 1.;
		    t8 = 1 / t7;
		    t9 = t8 * rhoa;
		    t10 = 1 / rho;
		    t11 = rhob * t10;
		    t14 = t5 * .2533;
		    t15 = exp(-t14);
		    t16 = t15 * t8;
/* Computing 2nd power */
		    d__1 = rho;
		    t17 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t19 = d__1 * d__1;
		    t21 = 1 / t19 / t17 / rho;
		    t22 = rhoa * rhob;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t23 = d__1 * d__1;
		    t24 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t34 = t5 * t8;
		    t36 = 2.611111111111111 - t5 * .09850555555555556 - t34 * 
			    .1357222222222222;
		    t41 = sigmaaa + sigmabb;
		    t45 = t14 + t34 * .349 - 11.;
		    t49 = rhoa * t10 * sigmaaa + t11 * sigmabb;
		    t52 = t25 * 36.46239897876478 * t23 + t30 * 
			    36.46239897876478 * t28 + t36 * sigma - (2.5 - t5 
			    * .01407222222222222 - t34 * .01938888888888889) *
			     1. * t41 - t45 * .1111111111111111 * t49;
		    t56 = t17 * .6666666666666667;
		    t57 = t23 * 1.;
		    t60 = t28 * 1.;
		    t63 = t22 * t52 - t17 * .6666666666666667 * sigma + (t56 
			    - t57) * sigmabb + (t56 - t60) * sigmaaa;
		    zk[i__] = t9 * -.19672 * t11 - t16 * .00649176 * t21 * 
			    t63;
		    t73 = t45 * t10;
/* Computing 2nd power */
		    d__1 = t7;
		    t84 = d__1 * d__1;
		    t85 = 1 / t84;
		    t91 = t85 * .02288509333333333 * rhoa * rhob / t4 / t17;
		    t92 = 1 / t17;
		    t93 = rhob * t92;
		    t95 = t9 * .19672 * t93;
/* Computing 2nd power */
		    d__1 = t17;
		    t96 = d__1 * d__1;
		    t98 = 1 / t96 / rho;
		    t102 = t98 * 5.48120936e-4 * t15 * t8 * t63;
		    t106 = t15 * 7.5520808e-4 * t85 * t98 * t63;
		    t111 = t16 * .02380312 / t19 / t96 * t63;
		    t113 = 1 / t4 / rho;
		    t115 = t113 * t8;
		    t119 = 1 / t19 / rho * t85;
		    t154 = t16 * .00649176 * t21 * (t22 * ((t113 * 
			    .03283518518518519 + t115 * .04524074074074074 - 
			    t119 * .01578901851851852) * sigma - (t113 * 
			    .004690740740740741 + t115 * .006462962962962963 
			    - t119 * .002255574074074074) * 1. * t41 - (t113 *
			     -.08443333333333333 - t115 * .1163333333333333 + 
			    t119 * .04060033333333333) * .1111111111111111 * 
			    t49 - t45 * .1111111111111111 * (rhoa * -1. * t92 
			    * sigmaaa - t93 * 1. * sigmabb)) - rho * 
			    1.333333333333333 * sigma + rho * 
			    1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa);
		    vrhoa[i__] = t8 * -.19672 * rhob * t10 - t16 * .00649176 *
			     t21 * (rhob * t52 + t22 * (t25 * 
			    97.23306394337274 * rhoa - t73 * 
			    .1111111111111111 * sigmaaa) - rhoa * 2. * 
			    sigmabb) - t91 + t95 - t102 - t106 + t111 - t154;
		    vrhob[i__] = t9 * -.19672 * t10 - t16 * .00649176 * t21 * 
			    (rhoa * t52 + t22 * (t30 * 97.23306394337274 * 
			    rhob - t73 * .1111111111111111 * sigmabb) - rhob *
			     2. * sigmaaa) - t91 + t95 - t102 - t106 + t111 - 
			    t154;
		    t170 = t5 * .01407222222222222;
		    t171 = t34 * .01938888888888889;
		    t185 = t16 * t21 * (t22 * t36 - t17 * .6666666666666667);
		    t186 = t185 * .00649176;
		    vsigmaaa[i__] = t16 * -.00649176 * t21 * (t22 * (t170 - 
			    2.5 + t171 - t45 * .1111111111111111 * rhoa * t10)
			     + t56 - t60) - t186;
		    vsigmaab[i__] = t185 * -.01298352;
		    vsigmabb[i__] = t16 * -.00649176 * t21 * (t22 * (t170 - 
			    2.5 + t171 - t45 * .1111111111111111 * rhob * t10)
			     + t56 - t57) - t186;
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
		    zk[i__] = 0.;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    v2rhob2[i__] = 0.;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    v2sigmabb2[i__] = 0.;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    zk[i__] = 0.;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    v2rhoa2[i__] = 0.;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    v2sigmaaa2[i__] = 0.;
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
		    t4 = pow_dd(&rho, &c_b2);
		    t5 = 1 / t4;
		    t7 = t5 * .349 + 1.;
		    t8 = 1 / t7;
		    t9 = t8 * rhoa;
		    t10 = 1 / rho;
		    t11 = rhob * t10;
		    t14 = t5 * .2533;
		    t15 = exp(-t14);
		    t16 = t15 * t8;
/* Computing 2nd power */
		    d__1 = rho;
		    t17 = d__1 * d__1;
		    t18 = t17 * rho;
/* Computing 2nd power */
		    d__1 = t4;
		    t19 = d__1 * d__1;
		    t21 = 1 / t19 / t18;
		    t22 = rhoa * rhob;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t23 = d__1 * d__1;
		    t24 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = d__1 * d__1;
		    t27 = t25 * 36.46239897876478 * t23;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t32 = t30 * 36.46239897876478 * t28;
		    t34 = t5 * t8;
		    t36 = 2.611111111111111 - t5 * .09850555555555556 - t34 * 
			    .1357222222222222;
		    t37 = t36 * sigma;
		    t41 = sigmaaa + sigmabb;
		    t43 = (2.5 - t5 * .01407222222222222 - t34 * 
			    .01938888888888889) * 1. * t41;
		    t45 = t14 + t34 * .349 - 11.;
		    t49 = rhoa * t10 * sigmaaa + t11 * sigmabb;
		    t51 = t45 * .1111111111111111 * t49;
		    t52 = t27 + t32 + t37 - t43 - t51;
		    t56 = t17 * .6666666666666667;
		    t57 = t23 * 1.;
		    t60 = t28 * 1.;
		    t63 = t22 * t52 - t17 * .6666666666666667 * sigma + (t56 
			    - t57) * sigmabb + (t56 - t60) * sigmaaa;
		    zk[i__] = t9 * -.19672 * t11 - t16 * .00649176 * t21 * 
			    t63;
		    t67 = t8 * rhob;
		    t71 = t25 * rhoa;
		    t73 = t45 * t10;
		    t76 = t71 * 97.23306394337274 - t73 * .1111111111111111 * 
			    sigmaaa;
		    t80 = rhob * t52 + t22 * t76 - rhoa * 2. * sigmabb;
/* Computing 2nd power */
		    d__1 = t7;
		    t84 = d__1 * d__1;
		    t85 = 1 / t84;
		    t86 = t85 * rhoa;
		    t88 = 1 / t4 / t17;
		    t91 = t86 * .02288509333333333 * rhob * t88;
		    t92 = 1 / t17;
		    t93 = rhob * t92;
		    t95 = t9 * .19672 * t93;
/* Computing 2nd power */
		    d__1 = t17;
		    t96 = d__1 * d__1;
		    t97 = t96 * rho;
		    t98 = 1 / t97;
		    t99 = t98 * t15;
		    t100 = t8 * t63;
		    t102 = t99 * 5.48120936e-4 * t100;
		    t103 = t15 * t85;
		    t106 = t103 * 7.5520808e-4 * t98 * t63;
		    t108 = 1 / t19 / t96;
		    t111 = t16 * .02380312 * t108 * t63;
		    t113 = 1 / t4 / rho;
		    t115 = t113 * t8;
		    t119 = 1 / t19 / rho * t85;
		    t121 = t113 * .03283518518518519 + t115 * 
			    .04524074074074074 - t119 * .01578901851851852;
		    t132 = t113 * -.08443333333333333 - t115 * 
			    .1163333333333333 + t119 * .04060033333333333;
		    t140 = rhoa * -1. * t92 * sigmaaa - t93 * 1. * sigmabb;
		    t143 = t121 * sigma - (t113 * .004690740740740741 + t115 *
			     .006462962962962963 - t119 * .002255574074074074)
			     * 1. * t41 - t132 * .1111111111111111 * t49 - 
			    t45 * .1111111111111111 * t140;
		    t151 = t22 * t143 - rho * 1.333333333333333 * sigma + rho 
			    * 1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa;
		    t154 = t16 * .00649176 * t21 * t151;
		    vrhoa[i__] = t67 * -.19672 * t10 - t16 * .00649176 * t21 *
			     t80 - t91 + t95 - t102 - t106 + t111 - t154;
		    t158 = t30 * rhob;
		    t162 = t158 * 97.23306394337274 - t73 * .1111111111111111 
			    * sigmabb;
		    t166 = rhoa * t52 + t22 * t162 - rhob * 2. * sigmaaa;
		    vrhob[i__] = t9 * -.19672 * t10 - t16 * .00649176 * t21 * 
			    t166 - t91 + t95 - t102 - t106 + t111 - t154;
		    t170 = t5 * .01407222222222222;
		    t171 = t34 * .01938888888888889;
		    t172 = t45 * rhoa;
		    t175 = t170 - 2.5 + t171 - t172 * .1111111111111111 * t10;
		    t177 = t22 * t175 + t56 - t60;
		    t183 = t22 * t36 - t17 * .6666666666666667;
		    t185 = t16 * t21 * t183;
		    t186 = t185 * .00649176;
		    vsigmaaa[i__] = t16 * -.00649176 * t21 * t177 - t186;
		    vsigmaab[i__] = t185 * -.01298352;
		    t188 = t45 * rhob;
		    t191 = t170 - 2.5 + t171 - t188 * .1111111111111111 * t10;
		    t193 = t22 * t191 + t56 - t57;
		    vsigmabb[i__] = t16 * -.00649176 * t21 * t193 - t186;
		    t198 = t85 * rhob * t88;
		    t200 = t67 * t92;
		    t203 = t99 * t8 * t80;
		    t206 = t103 * t98 * t80;
		    t209 = t16 * t108 * t80;
		    t212 = t132 * t10;
		    t215 = t45 * t92;
		    t222 = t16 * t21 * (rhob * t143 + t22 * (t212 * 
			    -.1111111111111111 * sigmaaa + t215 * 
			    .1111111111111111 * sigmaaa));
		    t237 = t16 * .1110812266666667 / t19 / t97 * t63;
		    t239 = t88 * t8;
		    t243 = 1 / t19 / t17 * t85;
		    t245 = 1 / t18;
		    t247 = 1 / t84 / t7;
		    t248 = t245 * t247;
		    t271 = rhob * t245;
		    t285 = t16 * .00649176 * t21 * (t22 * ((t88 * 
			    -.04378024691358025 - t239 * .06032098765432099 + 
			    t243 * .03157803703703704 - t248 * 
			    .003673578308641975) * sigma - (t88 * 
			    -.006254320987654321 - t239 * .008617283950617284 
			    + t243 * .004511148148148148 - t248 * 
			    5.247969012345679e-4) * 1. * t41 - (t88 * 
			    .1125777777777778 + t239 * .1551111111111111 - 
			    t243 * .08120066666666667 + t248 * 
			    .009446344222222222) * .1111111111111111 * t49 - 
			    t132 * .2222222222222222 * t140 - t45 * 
			    .1111111111111111 * (rhoa * 2. * t245 * sigmaaa + 
			    t271 * 2. * sigmabb)) - sigma * 1.333333333333333 
			    + sigmabb * 1.333333333333333 + sigmaaa * 
			    1.333333333333333);
		    t286 = t96 * t17;
		    t287 = 1 / t286;
		    t290 = t103 * .006545136693333333 * t287 * t63;
		    t293 = t287 * .004750381445333333 * t15 * t100;
		    t298 = t86 * .07628364444444444 * rhob / t4 / t18;
		    t300 = rhob * t21;
		    t302 = t247 * .005324598382222222 * rhoa * t300;
		    t305 = t99 * .001096241872 * t8 * t151;
		    t307 = 1 / t4 / t286;
		    t308 = t307 * t15;
		    t310 = t308 * 4.627967769626667e-5 * t100;
		    t313 = t308 * 1.275294711093333e-4 * t85 * t63;
		    t316 = t103 * .00151041616 * t98 * t151;
		    t320 = t15 * 1.757117466133333e-4 * t247 * t307 * t63;
		    t323 = t16 * .04760624 * t108 * t151;
		    t325 = t9 * .39344 * t271;
		    v2rhoa2[i__] = t198 * -.04577018666666667 + t200 * .39344 
			    - t203 * .001096241872 - t206 * .00151041616 + 
			    t209 * .04760624 - t222 * .01298352 - t16 * 
			    .00649176 * t21 * (rhob * 2. * t76 + t71 * 
			    162.0551065722879 * rhob - sigmabb * 2.) - t237 - 
			    t285 + t290 + t293 + t298 - t302 - t305 - t310 - 
			    t313 - t316 - t320 + t323 - t325;
		    t335 = t86 * t88;
		    t337 = t9 * t92;
		    t340 = t99 * t8 * t166;
		    t343 = t103 * t98 * t166;
		    t346 = t16 * t108 * t166;
		    t357 = t16 * t21 * (rhoa * t143 + t22 * (t212 * 
			    -.1111111111111111 * sigmabb + t215 * 
			    .1111111111111111 * sigmabb));
		    v2rhob2[i__] = -t237 - t285 + t290 + t293 + t298 - t302 - 
			    t305 - t310 - t313 - t316 - t320 + t323 - t16 * 
			    .00649176 * t21 * (rhoa * 2. * t162 + rhoa * 
			    162.0551065722879 * t158 - sigmaaa * 2.) - t335 * 
			    .04577018666666667 + t337 * .39344 - t340 * 
			    .001096241872 - t343 * .00151041616 + t346 * 
			    .04760624 - t357 * .01298352 - t325;
		    t377 = t337 * .19672 - t335 * .02288509333333333 - t340 * 
			    5.48120936e-4 - t343 * 7.5520808e-4 + t346 * 
			    .02380312 - t8 * .19672 * t10 - t16 * .00649176 * 
			    t21 * (t27 + t32 + t37 - t43 - t51 + rhob * t162 
			    + rhoa * t76) - t357 * .00649176 - t203 * 
			    5.48120936e-4 - t206 * 7.5520808e-4 + t209 * 
			    .02380312 - t222 * .00649176 - t237;
		    t380 = -t285 + t290 + t293 + t298 - t302 - t305 - t310 - 
			    t313 - t316 - t320 + t323 - t325 - t198 * 
			    .02288509333333333 + t200 * .19672;
		    v2rhoab[i__] = t377 + t380;
		    t383 = t22 * .1111111111111111 * t73;
		    t390 = t99 * 5.48120936e-4 * t8 * t177;
		    t393 = t103 * 7.5520808e-4 * t98 * t177;
		    t396 = t16 * .02380312 * t108 * t177;
		    t397 = t113 * .004690740740740741;
		    t398 = t115 * .006462962962962963;
		    t399 = t119 * .002255574074074074;
		    t407 = rho * 1.333333333333333;
		    t411 = t16 * .00649176 * t21 * (t22 * (-t397 - t398 + 
			    t399 - t132 * .1111111111111111 * rhoa * t10 + 
			    t172 * .1111111111111111 * t92) + t407);
		    t413 = t16 * t300 * t36;
		    t414 = t413 * .00649176;
		    t416 = t99 * t8 * t183;
		    t417 = t416 * 5.48120936e-4;
		    t419 = t103 * t98 * t183;
		    t420 = t419 * 7.5520808e-4;
		    t422 = t16 * t108 * t183;
		    t423 = t422 * .02380312;
		    t428 = t16 * t21 * (t22 * t121 - rho * 1.333333333333333);
		    t429 = t428 * .00649176;
		    v2rhoasigmaaa[i__] = t16 * -.00649176 * t21 * (rhob * 
			    t175 - t383) - t390 - t393 + t396 - t411 - t414 - 
			    t417 - t420 + t423 - t429;
		    t431 = t416 * .001096241872;
		    t432 = t419 * .00151041616;
		    t433 = t422 * .04760624;
		    t434 = t428 * .01298352;
		    v2rhoasigmaab[i__] = t413 * -.01298352 - t431 - t432 + 
			    t433 - t434;
		    t443 = t99 * 5.48120936e-4 * t8 * t193;
		    t446 = t103 * 7.5520808e-4 * t98 * t193;
		    t449 = t16 * .02380312 * t108 * t193;
		    t460 = t16 * .00649176 * t21 * (t22 * (-t397 - t398 + 
			    t399 - t132 * .1111111111111111 * rhob * t10 + 
			    t188 * .1111111111111111 * t92) + t407);
		    v2rhoasigmabb[i__] = t16 * -.00649176 * t21 * (rhob * 
			    t191 - rhoa * 2.) - t443 - t446 + t449 - t460 - 
			    t414 - t417 - t420 + t423 - t429;
		    t469 = t16 * t21 * rhoa * t36;
		    t470 = t469 * .00649176;
		    v2rhobsigmaaa[i__] = t16 * -.00649176 * t21 * (rhoa * 
			    t175 - rhob * 2.) - t390 - t393 + t396 - t411 - 
			    t470 - t417 - t420 + t423 - t429;
		    v2rhobsigmaab[i__] = t469 * -.01298352 - t431 - t432 + 
			    t433 - t434;
		    v2rhobsigmabb[i__] = t16 * -.00649176 * t21 * (rhoa * 
			    t191 - t383) - t443 - t446 + t449 - t460 - t470 - 
			    t417 - t420 + t423 - t429;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    v2sigmabb2[i__] = 0.;
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
} /* uks_c_lyp__ */

/* Subroutine */ int rks_c_lyp__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t5, t6, t9, t10, t20, t12, t11, t14, 
	    t22, t16, t30, t33, t45, t38, t54, t49, t60, t61, t64, t66, t81, 
	    t83, t85, t86, t13, t17, t18, t23, t28, t32, t47, t48, t51, t56, 
	    t65, t67, t68, t71, t76, t88, t99, t103, t121, t141, t150, t107, 
	    t116, t126, t127, t137, t165, t166, t190, t192, t188, t193, t195, 
	    t196, rho, sigma;


/*     C. Lee, W. Yang, and R.G. Parr */
/*     Development of the Colle-Salvetti correlation-energy formula into */
/*     a functional of the electron density */
/*     Phys. Rev. B37 (1988) 785-789 */


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
		t3 = 1 / t2;
		t6 = 1 / (t3 * .349 + 1.);
		t9 = t3 * .2533;
		t10 = exp(-t9);
/* Computing 2nd power */
		d__1 = rho;
		t12 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t14 = d__1 * d__1;
		t20 = t3 * t6;
		zk[i__] = t6 * -.04918 * rho - t10 * .00649176 * t6 / t14 / 
			t12 / rho * (t12 * .25 * (t14 * 11.48493600075277 * 
			t12 + (2.611111111111111 - t3 * .09850555555555556 - 
			t20 * .1357222222222222) * sigma - (2.5 - t3 * 
			.01407222222222222 - t20 * .01938888888888889) * .5 * 
			sigma - (t9 + t20 * .349 - 11.) * .02777777777777778 *
			 sigma) - t12 * .4583333333333333 * sigma);
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
		t3 = 1 / t2;
		t5 = t3 * .349 + 1.;
		t6 = 1 / t5;
		t9 = t3 * .2533;
		t10 = exp(-t9);
		t11 = t10 * t6;
/* Computing 2nd power */
		d__1 = rho;
		t12 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t14 = d__1 * d__1;
		t16 = 1 / t14 / t12 / rho;
		t20 = t3 * t6;
		t22 = 2.611111111111111 - t3 * .09850555555555556 - t20 * 
			.1357222222222222;
		t30 = t9 + t20 * .349 - 11.;
		t33 = t14 * 11.48493600075277 * t12 + t22 * sigma - (2.5 - t3 
			* .01407222222222222 - t20 * .01938888888888889) * .5 
			* sigma - t30 * .02777777777777778 * sigma;
		t38 = t12 * .25 * t33 - t12 * .4583333333333333 * sigma;
		zk[i__] = t6 * -.04918 * rho - t11 * .00649176 * t16 * t38;
		t45 = t14 * rho;
		t49 = t30 / rho * sigma;
		t54 = rho * sigma;
/* Computing 2nd power */
		d__1 = t5;
		t60 = d__1 * d__1;
		t61 = 1 / t60;
/* Computing 2nd power */
		d__1 = t12;
		t64 = d__1 * d__1;
		t66 = 1 / t64 / rho;
		t81 = 1 / t2 / rho;
		t83 = t81 * t6;
		t85 = 1 / t45;
		t86 = t85 * t61;
		vrhoa[i__] = t6 * -.04918 - t11 * .00649176 * t16 * (rho * .5 
			* t33 + t12 * .25 * (t45 * 30.62649600200738 - t49 * 
			.02777777777777778) - t54 * .25) - t61 * 
			.005721273333333333 * t3 - t66 * 5.48120936e-4 * t10 *
			 t6 * t38 - t10 * 7.5520808e-4 * t61 * t66 * t38 + 
			t11 * .02380312 / t14 / t64 * t38 - t11 * .00649176 * 
			t16 * (t12 * .25 * ((t81 * .03283518518518519 + t83 * 
			.04524074074074074 - t86 * .01578901851851852) * 
			sigma - (t81 * .004690740740740741 + t83 * 
			.006462962962962963 - t86 * .002255574074074074) * .5 
			* sigma - (t81 * -.08443333333333333 - t83 * 
			.1163333333333333 + t86 * .04060033333333333) * 
			.02777777777777778 * sigma + t49 * .02777777777777778)
			 - t54 * .6666666666666667);
		vsigmaaa[i__] = t11 * 7.213066666666667e-4 * t85 - t11 * 
			.02596704 * t16 * (t12 * .25 * t22 - t12 * 
			.6666666666666667);
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
		t3 = 1 / t2;
		t5 = t3 * .349 + 1.;
		t6 = 1 / t5;
		t9 = t3 * .2533;
		t10 = exp(-t9);
		t11 = t10 * t6;
/* Computing 2nd power */
		d__1 = rho;
		t12 = d__1 * d__1;
		t13 = t12 * rho;
/* Computing 2nd power */
		d__1 = t2;
		t14 = d__1 * d__1;
		t16 = 1 / t14 / t13;
		t17 = t14 * t12;
		t18 = t17 * 11.48493600075277;
		t20 = t3 * t6;
		t22 = 2.611111111111111 - t3 * .09850555555555556 - t20 * 
			.1357222222222222;
		t23 = t22 * sigma;
		t28 = (2.5 - t3 * .01407222222222222 - t20 * 
			.01938888888888889) * .5 * sigma;
		t30 = t9 + t20 * .349 - 11.;
		t32 = t30 * .02777777777777778 * sigma;
		t33 = t18 + t23 - t28 - t32;
		t38 = t12 * .25 * t33 - t12 * .4583333333333333 * sigma;
		zk[i__] = t6 * -.04918 * rho - t11 * .00649176 * t16 * t38;
		t45 = t14 * rho;
		t47 = 1 / rho;
		t48 = t30 * t47;
		t49 = t48 * sigma;
		t51 = t45 * 30.62649600200738 - t49 * .02777777777777778;
		t54 = rho * sigma;
		t56 = rho * .5 * t33 + t12 * .25 * t51 - t54 * .25;
/* Computing 2nd power */
		d__1 = t5;
		t60 = d__1 * d__1;
		t61 = 1 / t60;
/* Computing 2nd power */
		d__1 = t12;
		t64 = d__1 * d__1;
		t65 = t64 * rho;
		t66 = 1 / t65;
		t67 = t66 * t10;
		t68 = t6 * t38;
		t71 = t10 * t61;
		t76 = 1 / t14 / t64;
		t81 = 1 / t2 / rho;
		t83 = t81 * t6;
		t85 = 1 / t45;
		t86 = t85 * t61;
		t88 = t81 * .03283518518518519 + t83 * .04524074074074074 - 
			t86 * .01578901851851852;
		t99 = t81 * -.08443333333333333 - t83 * .1163333333333333 + 
			t86 * .04060033333333333;
		t103 = t88 * sigma - (t81 * .004690740740740741 + t83 * 
			.006462962962962963 - t86 * .002255574074074074) * .5 
			* sigma - t99 * .02777777777777778 * sigma + t49 * 
			.02777777777777778;
		t107 = t12 * .25 * t103 - t54 * .6666666666666667;
		vrhoa[i__] = t6 * -.04918 - t11 * .00649176 * t16 * t56 - t61 
			* .005721273333333333 * t3 - t67 * 5.48120936e-4 * 
			t68 - t71 * 7.5520808e-4 * t66 * t38 + t11 * 
			.02380312 * t76 * t38 - t11 * .00649176 * t16 * t107;
		t116 = t12 * .25 * t22 - t12 * .6666666666666667;
		vsigmaaa[i__] = t11 * 7.213066666666667e-4 * t85 - t11 * 
			.02596704 * t16 * t116;
		t121 = 1 / t60 / t5;
		t126 = t64 * t12;
		t127 = 1 / t126;
		t137 = t99 * t47 * sigma;
		t141 = t30 / t12 * sigma;
		t150 = rho * t51;
		t165 = 1 / t2 / t126;
		t166 = t165 * t10;
		t188 = 1 / t2 / t12;
		t190 = t188 * t6;
		t192 = 1 / t17;
		t193 = t192 * t61;
		t195 = 1 / t13;
		t196 = t195 * t121;
		s1 = t121 * -.002662299191111111 * t85 - t61 * 
			.007628364444444444 * t81 + t127 * 
			.009500762890666667 * t10 * t68 + t11 * .09521248 * 
			t76 * t56 - t11 * .02596704 * t16 * (rho * .5 * t103 
			+ t12 * .25 * (t137 * -.02777777777777778 + t141 * 
			.02777777777777778)) - t11 * .00649176 * t16 * (t18 + 
			t23 - t28 - t32 + t150) + t11 * .09521248 * t76 * 
			t107 - t67 * .002192483744 * t6 * t56 - t71 * 
			.00302083232 * t66 * t56;
		s2 = s1 - t166 * 2.550589422186667e-4 * t61 * t38 - t71 * 
			.00302083232 * t66 * t107 - t10 * 
			3.514234932266667e-4 * t121 * t165 * t38 - t67 * 
			.002192483744 * t6 * t107;
		v2rhoa2[i__] = s2 - t166 * 9.255935539253333e-5 * t68 - t11 * 
			.2221624533333333 / t14 / t65 * t38 - t11 * .01298352 
			* t16 * (t12 * .25 * ((t188 * -.04378024691358025 - 
			t190 * .06032098765432099 + t193 * .03157803703703704 
			- t196 * .003673578308641975) * sigma - (t188 * 
			-.006254320987654321 - t190 * .008617283950617284 + 
			t193 * .004511148148148148 - t196 * 
			5.247969012345679e-4) * .5 * sigma - (t188 * 
			.1125777777777778 + t190 * .1551111111111111 - t193 * 
			.08120066666666667 + t196 * .009446344222222222) * 
			.02777777777777778 * sigma + t137 * 
			.05555555555555556 - t141 * .05555555555555556) - 
			sigma * .6666666666666667) - t11 * .00649176 * t16 * (
			t150 * 1. + t17 * 25.52208000167282 - sigma * .5) + 
			t71 * .01309027338666667 * t127 * t38;
		v2rhoasigmaaa[i__] = t11 * -.00649176 * t16 * (rho * 
			-.9444444444444444 - rho * .02777777777777778 * t30) 
			+ t195 * 6.090232622222222e-5 * t10 * t6 + t71 * 
			8.391200888888889e-5 * t195 + t11 * 
			.009978075555555556 * t192 - t11 * .01298352 * t16 * (
			t12 * .25 * (t83 * -3e-23 + t86 * 1e-23 + t48 * 
			.05555555555555556) + rho * 1.333333333333333) - t11 *
			 .01298352 * t192 * t22 - t67 * .002192483744 * t6 * 
			t116 - t71 * .00302083232 * t66 * t116 + t11 * 
			.09521248 * t76 * t116 - t11 * .02596704 * t16 * (t12 
			* .25 * t88 - rho * 1.333333333333333);
		v2sigmaaa2[i__] = 0.;
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
} /* rks_c_lyp__ */


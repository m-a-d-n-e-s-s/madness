/* xc_xc_edf1.f -- translated by f2c (version 20050501).
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

/* :XC_EDF1subrstart */
/*    Generated: Wed Jul 23 17:14:37 GMT 2003 */
/* Subroutine */ int uks_xc_edf1__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14,
	     t15, t30, t31, t33, t28, t29, t34, t35, t37, t38, t39, t52, t53, 
	    t56, t58, t59, t62, t63, t65, t67, t74, t96, t16, t20, t21, t24, 
	    t32, t36, t40, t41, t49, t50, t17, t18, t22, t23, t26, t42, t46, 
	    t47, t55, t57, t64, t69, t70, t76, t81, t85, t89, t92, t97, t100, 
	    t103, t110, t113, t114, t118, t121, t122, t123, t131, t132, t149, 
	    t160, t161, t167, t168, t169, t171, t172, t174, t178, t182, t187, 
	    t189, t191, t195, t230, t234, t237, t238, t242, t245, t246, t247, 
	    t255, t256, t286, t288, t304, t305, t319, t320, t326, t328, t25, 
	    t43, t44, t54, t77, t82, t91, t99, t101, t109, t143, t150, t154, 
	    t159, t66, rho, t71, t72, t83, t107, t116, t120, t125, t126, t135,
	     t136, t139, t147, t152, t156, t162, t164, t173, t175, t176, t179,
	     t184, t197, t208, t216, t219, t227, t231, t233, t240, t244, t249,
	     t250, t259, t260, t263, t270, t274, t278, t284, t290, t291, t298,
	     t299, t302, t306, t309, t311, t317, t324, t330, t331, t338, t339,
	     t342, t344, t347, t349, t362, t365, t371, t372, t379, t381, t383,
	     t389, t390, t401, t402, t416, t418, t420, t423, t424, t426, t427,
	     t429, t432, t435, t439, t442, t444, t448, t450, t451, t474, t488,
	     t490, t493, t500, t502, t504, t507, t510, rhoa, rhob, t512, t515,
	     t522, t525, t527, t530, t535, t540, t541, t542, t557, t559, t562,
	     t565, t567, t570, t576, t577, t584, t586, t588, t594, t595, t606,
	     t607, t621, t632, t634, t640, t656, t661, t669, t671, t677, t687,
	     t707, t714, t717, t720, t721, t722, t723, t731, t735, t737, t738,
	     t740, t741, t743, t744, t746, t747, t752, t753, t755, t756, t757,
	     t758, t767, t770, t773, t784, t793, t794, t800, t808, t810, t816,
	     t826, t851, t858, t862, t867, t875, t889, t896, t900, t905, t913,
	     sigma, sigmaaa, sigmaab, sigmabb;


/*     R.D. Adamson, P.M.W. Gill, and J.A. Pople */
/*     Empirical density functionals */
/*     Chem. Phys. Lett. 284 (1998) 6-11 */


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
/* Computing 2nd power */
		    d__1 = rhob;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t7 = sigmabb / t4 / t2;
		    t8 = sqrt(sigmabb);
		    t9 = t3 * rhob;
		    t11 = t8 / t9;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    zk[i__] = (-.9593273689405774 - t7 * .03640595 / (t13 * 
			    .021 + 1.) + t7 * .035481306 / (t13 * .0252 + 1.))
			     * t9;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t7 = sigmaaa / t4 / t2;
		    t8 = sqrt(sigmaaa);
		    t9 = t3 * rhoa;
		    t11 = t8 / t9;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    zk[i__] = (-.9593273689405774 - t7 * .03640595 / (t13 * 
			    .021 + 1.) + t7 * .035481306 / (t13 * .0252 + 1.))
			     * t9;
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
		    t5 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t5;
		    t6 = d__1 * d__1;
		    t7 = t6 * t4;
		    t9 = sigmaaa / t7;
		    t10 = sqrt(sigmaaa);
		    t11 = t5 * rhoa;
		    t13 = t10 / t11;
/* Computing 2nd power */
		    d__1 = t13;
		    t14 = log(t13 + sqrt(d__1 * d__1 + 1));
		    t15 = t13 * t14;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t31 = t30 * t28;
		    t33 = sigmabb / t31;
		    t34 = sqrt(sigmabb);
		    t35 = t29 * rhob;
		    t37 = t34 / t35;
/* Computing 2nd power */
		    d__1 = t37;
		    t38 = log(t37 + sqrt(d__1 * d__1 + 1));
		    t39 = t37 * t38;
		    t52 = pow_dd(&rho, &c_b2);
		    t53 = 1 / t52;
		    t56 = 1 / (t53 * .3505 + 1.);
		    t58 = 1 / rho;
		    t59 = rhob * t58;
		    t62 = t53 * .25;
		    t63 = exp(-t62);
/* Computing 2nd power */
		    d__1 = rho;
		    t65 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t52;
		    t67 = d__1 * d__1;
		    t74 = t53 * t56;
		    t96 = t65 * .6666666666666667;
		    zk[i__] = (-.9593273689405774 - t9 * .03640595 / (t15 * 
			    .021 + 1.) + t9 * .035481306 / (t15 * .0252 + 1.))
			     * t11 + (-.9593273689405774 - t33 * .03640595 / (
			    t39 * .021 + 1.) + t33 * .035481306 / (t39 * 
			    .0252 + 1.)) * t35 - t56 * .22 * rhoa * t59 - t63 
			    * .00869 * t56 / t67 / t65 / rho * (rhoa * rhob * 
			    (t7 * 36.46239897876478 + t31 * 36.46239897876478 
			    + (2.611111111111111 - t53 * .09722222222222222 - 
			    t74 * .1363055555555556) * sigma - (2.5 - t53 * 
			    .01388888888888889 - t74 * .01947222222222222) * 
			    1. * (sigmaaa + sigmabb) - (t62 + t74 * .3505 - 
			    11.) * .1111111111111111 * (rhoa * t58 * sigmaaa 
			    + t59 * sigmabb)) - t65 * .6666666666666667 * 
			    sigma + (t96 - t4 * 1.) * sigmabb + (t96 - t28 * 
			    1.) * sigmaaa);
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
/* Computing 2nd power */
		    d__1 = rhob;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t6 = 1 / t4 / t2;
		    t7 = sigmabb * t6;
		    t8 = sqrt(sigmabb);
		    t9 = t3 * rhob;
		    t10 = 1 / t9;
		    t11 = t8 * t10;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    t15 = t13 * .021 + 1.;
		    t16 = 1 / t15;
		    t20 = t13 * .0252 + 1.;
		    t21 = 1 / t20;
		    t24 = -.9593273689405774 - t7 * .03640595 * t16 + t7 * 
			    .035481306 * t21;
		    zk[i__] = t24 * t9;
		    vrhoa[i__] = 0.;
		    t28 = sigmabb / t4 / t2 / rhob;
/* Computing 2nd power */
		    d__1 = t15;
		    t31 = d__1 * d__1;
		    t32 = 1 / t31;
		    t36 = t8 / t3 / t2 * t12;
		    t39 = sqrt(t7 + 1.);
		    t40 = 1 / t39;
		    t41 = t28 * t40;
/* Computing 2nd power */
		    d__1 = t20;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    vrhob[i__] = (t28 * .09708253333333333 * t16 + t7 * 
			    .03640595 * t32 * (t36 * -.028 - t41 * .028) - 
			    t28 * .094616816 * t21 - t7 * .035481306 * t50 * (
			    t36 * -.0336 - t41 * .0336)) * t9 + t24 * 
			    1.333333333333333 * t3;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t65 = 1 / t8 * t10 * t12;
		    t67 = t6 * t40;
		    vsigmabb[i__] = (t6 * -.03640595 * t16 + t7 * .03640595 * 
			    t32 * (t65 * .0105 + t67 * .0105) + t6 * 
			    .035481306 * t21 - t7 * .035481306 * t50 * (t65 * 
			    .0126 + t67 * .0126)) * t9;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t6 = 1 / t4 / t2;
		    t7 = sigmaaa * t6;
		    t8 = sqrt(sigmaaa);
		    t9 = t3 * rhoa;
		    t10 = 1 / t9;
		    t11 = t8 * t10;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    t15 = t13 * .021 + 1.;
		    t16 = 1 / t15;
		    t20 = t13 * .0252 + 1.;
		    t21 = 1 / t20;
		    t24 = -.9593273689405774 - t7 * .03640595 * t16 + t7 * 
			    .035481306 * t21;
		    zk[i__] = t24 * t9;
		    t28 = sigmaaa / t4 / t2 / rhoa;
/* Computing 2nd power */
		    d__1 = t15;
		    t31 = d__1 * d__1;
		    t32 = 1 / t31;
		    t36 = t8 / t3 / t2 * t12;
		    t39 = sqrt(t7 + 1.);
		    t40 = 1 / t39;
		    t41 = t28 * t40;
/* Computing 2nd power */
		    d__1 = t20;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    vrhoa[i__] = (t28 * .09708253333333333 * t16 + t7 * 
			    .03640595 * t32 * (t36 * -.028 - t41 * .028) - 
			    t28 * .094616816 * t21 - t7 * .035481306 * t50 * (
			    t36 * -.0336 - t41 * .0336)) * t9 + t24 * 
			    1.333333333333333 * t3;
		    vrhob[i__] = 0.;
		    t65 = 1 / t8 * t10 * t12;
		    t67 = t6 * t40;
		    vsigmaaa[i__] = (t6 * -.03640595 * t16 + t7 * .03640595 * 
			    t32 * (t65 * .0105 + t67 * .0105) + t6 * 
			    .035481306 * t21 - t7 * .035481306 * t50 * (t65 * 
			    .0126 + t67 * .0126)) * t9;
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
		    t5 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t5;
		    t6 = d__1 * d__1;
		    t7 = t6 * t4;
		    t8 = 1 / t7;
		    t9 = sigmaaa * t8;
		    t10 = sqrt(sigmaaa);
		    t11 = t5 * rhoa;
		    t12 = 1 / t11;
		    t13 = t10 * t12;
/* Computing 2nd power */
		    d__1 = t13;
		    t14 = log(t13 + sqrt(d__1 * d__1 + 1));
		    t15 = t13 * t14;
		    t17 = t15 * .021 + 1.;
		    t18 = 1 / t17;
		    t22 = t15 * .0252 + 1.;
		    t23 = 1 / t22;
		    t26 = -.9593273689405774 - t9 * .03640595 * t18 + t9 * 
			    .035481306 * t23;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t31 = t30 * t28;
		    t32 = 1 / t31;
		    t33 = sigmabb * t32;
		    t34 = sqrt(sigmabb);
		    t35 = t29 * rhob;
		    t36 = 1 / t35;
		    t37 = t34 * t36;
/* Computing 2nd power */
		    d__1 = t37;
		    t38 = log(t37 + sqrt(d__1 * d__1 + 1));
		    t39 = t37 * t38;
		    t41 = t39 * .021 + 1.;
		    t42 = 1 / t41;
		    t46 = t39 * .0252 + 1.;
		    t47 = 1 / t46;
		    t50 = -.9593273689405774 - t33 * .03640595 * t42 + t33 * 
			    .035481306 * t47;
		    t52 = pow_dd(&rho, &c_b2);
		    t53 = 1 / t52;
		    t55 = t53 * .3505 + 1.;
		    t56 = 1 / t55;
		    t57 = t56 * rhoa;
		    t58 = 1 / rho;
		    t59 = rhob * t58;
		    t62 = t53 * .25;
		    t63 = exp(-t62);
		    t64 = t63 * t56;
/* Computing 2nd power */
		    d__1 = rho;
		    t65 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t52;
		    t67 = d__1 * d__1;
		    t69 = 1 / t67 / t65 / rho;
		    t70 = rhoa * rhob;
		    t74 = t53 * t56;
		    t76 = 2.611111111111111 - t53 * .09722222222222222 - t74 *
			     .1363055555555556;
		    t81 = sigmaaa + sigmabb;
		    t85 = t62 + t74 * .3505 - 11.;
		    t89 = rhoa * t58 * sigmaaa + t59 * sigmabb;
		    t92 = t7 * 36.46239897876478 + t31 * 36.46239897876478 + 
			    t76 * sigma - (2.5 - t53 * .01388888888888889 - 
			    t74 * .01947222222222222) * 1. * t81 - t85 * 
			    .1111111111111111 * t89;
		    t96 = t65 * .6666666666666667;
		    t97 = t4 * 1.;
		    t100 = t28 * 1.;
		    t103 = t70 * t92 - t65 * .6666666666666667 * sigma + (t96 
			    - t97) * sigmabb + (t96 - t100) * sigmaaa;
		    zk[i__] = t26 * t11 + t50 * t35 - t57 * .22 * t59 - t64 * 
			    .00869 * t69 * t103;
		    t110 = sigmaaa / t6 / t4 / rhoa;
/* Computing 2nd power */
		    d__1 = t17;
		    t113 = d__1 * d__1;
		    t114 = 1 / t113;
		    t118 = t10 / t5 / t4 * t14;
		    t121 = sqrt(t9 + 1.);
		    t122 = 1 / t121;
		    t123 = t110 * t122;
/* Computing 2nd power */
		    d__1 = t22;
		    t131 = d__1 * d__1;
		    t132 = 1 / t131;
		    t149 = t85 * t58;
/* Computing 2nd power */
		    d__1 = t55;
		    t160 = d__1 * d__1;
		    t161 = 1 / t160;
		    t167 = t161 * .02570333333333333 * rhoa * rhob / t52 / 
			    t65;
		    t168 = 1 / t65;
		    t169 = rhob * t168;
		    t171 = t57 * .22 * t169;
/* Computing 2nd power */
		    d__1 = t65;
		    t172 = d__1 * d__1;
		    t174 = 1 / t172 / rho;
		    t178 = t174 * 7.241666666666667e-4 * t63 * t56 * t103;
		    t182 = t63 * .001015281666666667 * t161 * t174 * t103;
		    t187 = t64 * .03186333333333333 / t67 / t172 * t103;
		    t189 = 1 / t52 / rho;
		    t191 = t189 * t56;
		    t195 = 1 / t67 / rho * t161;
		    t230 = t64 * .00869 * t69 * (t70 * ((t189 * 
			    .03240740740740741 + t191 * .04543518518518519 - 
			    t195 * .01592503240740741) * sigma - (t189 * 
			    .00462962962962963 + t191 * .006490740740740741 - 
			    t195 * .00227500462962963) * 1. * t81 - (t189 * 
			    -.08333333333333333 - t191 * .1168333333333333 + 
			    t195 * .04095008333333333) * .1111111111111111 * 
			    t89 - t85 * .1111111111111111 * (rhoa * -1. * 
			    t168 * sigmaaa - t169 * 1. * sigmabb)) - rho * 
			    1.333333333333333 * sigma + rho * 
			    1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa);
		    vrhoa[i__] = (t110 * .09708253333333333 * t18 + t9 * 
			    .03640595 * t114 * (t118 * -.028 - t123 * .028) - 
			    t110 * .094616816 * t23 - t9 * .035481306 * t132 *
			     (t118 * -.0336 - t123 * .0336)) * t11 + t26 * 
			    1.333333333333333 * t5 - t56 * .22 * rhob * t58 - 
			    t64 * .00869 * t69 * (rhob * t92 + t70 * (t6 * 
			    97.23306394337274 * rhoa - t149 * 
			    .1111111111111111 * sigmaaa) - rhoa * 2. * 
			    sigmabb) - t167 + t171 - t178 - t182 + t187 - 
			    t230;
		    t234 = sigmabb / t30 / t28 / rhob;
/* Computing 2nd power */
		    d__1 = t41;
		    t237 = d__1 * d__1;
		    t238 = 1 / t237;
		    t242 = t34 / t29 / t28 * t38;
		    t245 = sqrt(t33 + 1.);
		    t246 = 1 / t245;
		    t247 = t234 * t246;
/* Computing 2nd power */
		    d__1 = t46;
		    t255 = d__1 * d__1;
		    t256 = 1 / t255;
		    vrhob[i__] = (t234 * .09708253333333333 * t42 + t33 * 
			    .03640595 * t238 * (t242 * -.028 - t247 * .028) - 
			    t234 * .094616816 * t47 - t33 * .035481306 * t256 
			    * (t242 * -.0336 - t247 * .0336)) * t35 + t50 * 
			    1.333333333333333 * t29 - t57 * .22 * t58 - t64 * 
			    .00869 * t69 * (rhoa * t92 + t70 * (t30 * 
			    97.23306394337274 * rhob - t149 * 
			    .1111111111111111 * sigmabb) - rhob * 2. * 
			    sigmaaa) - t167 + t171 - t178 - t182 + t187 - 
			    t230;
		    t286 = 1 / t10 * t12 * t14;
		    t288 = t8 * t122;
		    t304 = t53 * .01388888888888889;
		    t305 = t74 * .01947222222222222;
		    t319 = t64 * t69 * (t70 * t76 - t65 * .6666666666666667);
		    t320 = t319 * .00869;
		    vsigmaaa[i__] = (t8 * -.03640595 * t18 + t9 * .03640595 * 
			    t114 * (t286 * .0105 + t288 * .0105) + t8 * 
			    .035481306 * t23 - t9 * .035481306 * t132 * (t286 
			    * .0126 + t288 * .0126)) * t11 - t64 * .00869 * 
			    t69 * (t70 * (t304 - 2.5 + t305 - t85 * 
			    .1111111111111111 * rhoa * t58) + t96 - t100) - 
			    t320;
		    vsigmaab[i__] = t319 * -.01738;
		    t326 = 1 / t34 * t36 * t38;
		    t328 = t32 * t246;
		    vsigmabb[i__] = (t32 * -.03640595 * t42 + t33 * .03640595 
			    * t238 * (t326 * .0105 + t328 * .0105) + t32 * 
			    .035481306 * t47 - t33 * .035481306 * t256 * (
			    t326 * .0126 + t328 * .0126)) * t35 - t64 * 
			    .00869 * t69 * (t70 * (t304 - 2.5 + t305 - t85 * 
			    .1111111111111111 * rhob * t58) + t96 - t97) - 
			    t320;
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
/* Computing 2nd power */
		    d__1 = rhob;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t6 = 1 / t4 / t2;
		    t7 = sigmabb * t6;
		    t8 = sqrt(sigmabb);
		    t9 = t3 * rhob;
		    t10 = 1 / t9;
		    t11 = t8 * t10;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    t15 = t13 * .021 + 1.;
		    t16 = 1 / t15;
		    t20 = t13 * .0252 + 1.;
		    t21 = 1 / t20;
		    t24 = -.9593273689405774 - t7 * .03640595 * t16 + t7 * 
			    .035481306 * t21;
		    zk[i__] = t24 * t9;
		    vrhoa[i__] = 0.;
		    t25 = t2 * rhob;
		    t28 = sigmabb / t4 / t25;
/* Computing 2nd power */
		    d__1 = t15;
		    t31 = d__1 * d__1;
		    t32 = 1 / t31;
		    t36 = t8 / t3 / t2 * t12;
		    t38 = t7 + 1.;
		    t39 = sqrt(t38);
		    t40 = 1 / t39;
		    t41 = t28 * t40;
		    t43 = t36 * -.028 - t41 * .028;
		    t44 = t32 * t43;
/* Computing 2nd power */
		    d__1 = t20;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t36 * -.0336 - t41 * .0336;
		    t54 = t50 * t53;
		    t57 = t28 * .09708253333333333 * t16 + t7 * .03640595 * 
			    t44 - t28 * .094616816 * t21 - t7 * .035481306 * 
			    t54;
		    vrhob[i__] = t57 * t9 + t24 * 1.333333333333333 * t3;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t65 = 1 / t8 * t10 * t12;
		    t67 = t6 * t40;
		    t69 = t65 * .0105 + t67 * .0105;
		    t77 = t65 * .0126 + t67 * .0126;
		    vsigmabb[i__] = (t6 * -.03640595 * t16 + t7 * .03640595 * 
			    t32 * t69 + t6 * .035481306 * t21 - t7 * 
			    .035481306 * t50 * t77) * t9;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t2;
		    t82 = d__1 * d__1;
		    t85 = sigmabb / t4 / t82;
		    t91 = 1 / t31 / t15;
/* Computing 2nd power */
		    d__1 = t43;
		    t92 = d__1 * d__1;
		    t99 = t8 / t3 / t25 * t12;
		    t101 = t85 * t40;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t103 = d__1 * d__1;
		    t109 = 1 / t39 / t38;
		    t110 = t103 / t3 / t82 / t25 * t109;
		    t121 = 1 / t49 / t20;
/* Computing 2nd power */
		    d__1 = t53;
		    t122 = d__1 * d__1;
		    v2rhob2[i__] = (t85 * -.3559692888888889 * t16 - t28 * 
			    .1941650666666667 * t44 - t7 * .0728119 * t91 * 
			    t92 + t7 * .03640595 * t32 * (t99 * 
			    .06533333333333333 + t101 * .14 - t110 * 
			    .03733333333333333) + t85 * .3469283253333333 * 
			    t21 + t28 * .189233632 * t54 + t7 * .070962612 * 
			    t121 * t122 - t7 * .035481306 * t50 * (t99 * 
			    .0784 + t101 * .168 - t110 * .0448)) * t9 + t57 * 
			    2.666666666666667 * t3 + t24 * .4444444444444444 /
			     t4;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t69;
		    t143 = d__1 * d__1;
		    t150 = 1 / t8 / sigmabb * t10 * t12;
		    t154 = 1 / sigmabb * t6 * t40;
		    t159 = 1 / t3 / t82 / rhob * t109;
/* Computing 2nd power */
		    d__1 = t77;
		    t168 = d__1 * d__1;
		    v2sigmabb2[i__] = (t6 * .0728119 * t32 * t69 - t7 * 
			    .0728119 * t91 * t143 + t7 * .03640595 * t32 * (
			    t150 * -.00525 + t154 * .00525 - t159 * .00525) - 
			    t6 * .070962612 * t50 * t77 + t7 * .070962612 * 
			    t121 * t168 - t7 * .035481306 * t50 * (t150 * 
			    -.0063 + t154 * .0063 - t159 * .0063)) * t9;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t2 = d__1 * d__1;
		    t3 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t3;
		    t4 = d__1 * d__1;
		    t6 = 1 / t4 / t2;
		    t7 = sigmaaa * t6;
		    t8 = sqrt(sigmaaa);
		    t9 = t3 * rhoa;
		    t10 = 1 / t9;
		    t11 = t8 * t10;
/* Computing 2nd power */
		    d__1 = t11;
		    t12 = log(t11 + sqrt(d__1 * d__1 + 1));
		    t13 = t11 * t12;
		    t15 = t13 * .021 + 1.;
		    t16 = 1 / t15;
		    t20 = t13 * .0252 + 1.;
		    t21 = 1 / t20;
		    t24 = -.9593273689405774 - t7 * .03640595 * t16 + t7 * 
			    .035481306 * t21;
		    zk[i__] = t24 * t9;
		    t25 = t2 * rhoa;
		    t28 = sigmaaa / t4 / t25;
/* Computing 2nd power */
		    d__1 = t15;
		    t31 = d__1 * d__1;
		    t32 = 1 / t31;
		    t36 = t8 / t3 / t2 * t12;
		    t38 = t7 + 1.;
		    t39 = sqrt(t38);
		    t40 = 1 / t39;
		    t41 = t28 * t40;
		    t43 = t36 * -.028 - t41 * .028;
		    t44 = t32 * t43;
/* Computing 2nd power */
		    d__1 = t20;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t36 * -.0336 - t41 * .0336;
		    t54 = t50 * t53;
		    t57 = t28 * .09708253333333333 * t16 + t7 * .03640595 * 
			    t44 - t28 * .094616816 * t21 - t7 * .035481306 * 
			    t54;
		    vrhoa[i__] = t57 * t9 + t24 * 1.333333333333333 * t3;
		    vrhob[i__] = 0.;
		    t65 = 1 / t8 * t10 * t12;
		    t67 = t6 * t40;
		    t69 = t65 * .0105 + t67 * .0105;
		    t77 = t65 * .0126 + t67 * .0126;
		    vsigmaaa[i__] = (t6 * -.03640595 * t16 + t7 * .03640595 * 
			    t32 * t69 + t6 * .035481306 * t21 - t7 * 
			    .035481306 * t50 * t77) * t9;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t2;
		    t82 = d__1 * d__1;
		    t85 = sigmaaa / t4 / t82;
		    t91 = 1 / t31 / t15;
/* Computing 2nd power */
		    d__1 = t43;
		    t92 = d__1 * d__1;
		    t99 = t8 / t3 / t25 * t12;
		    t101 = t85 * t40;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t103 = d__1 * d__1;
		    t109 = 1 / t39 / t38;
		    t110 = t103 / t3 / t82 / t25 * t109;
		    t121 = 1 / t49 / t20;
/* Computing 2nd power */
		    d__1 = t53;
		    t122 = d__1 * d__1;
		    v2rhoa2[i__] = (t85 * -.3559692888888889 * t16 - t28 * 
			    .1941650666666667 * t44 - t7 * .0728119 * t91 * 
			    t92 + t7 * .03640595 * t32 * (t99 * 
			    .06533333333333333 + t101 * .14 - t110 * 
			    .03733333333333333) + t85 * .3469283253333333 * 
			    t21 + t28 * .189233632 * t54 + t7 * .070962612 * 
			    t121 * t122 - t7 * .035481306 * t50 * (t99 * 
			    .0784 + t101 * .168 - t110 * .0448)) * t9 + t57 * 
			    2.666666666666667 * t3 + t24 * .4444444444444444 /
			     t4;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t69;
		    t143 = d__1 * d__1;
		    t150 = 1 / t8 / sigmaaa * t10 * t12;
		    t154 = 1 / sigmaaa * t6 * t40;
		    t159 = 1 / t3 / t82 / rhoa * t109;
/* Computing 2nd power */
		    d__1 = t77;
		    t168 = d__1 * d__1;
		    v2sigmaaa2[i__] = (t6 * .0728119 * t32 * t69 - t7 * 
			    .0728119 * t91 * t143 + t7 * .03640595 * t32 * (
			    t150 * -.00525 + t154 * .00525 - t159 * .00525) - 
			    t6 * .070962612 * t50 * t77 + t7 * .070962612 * 
			    t121 * t168 - t7 * .035481306 * t50 * (t150 * 
			    -.0063 + t154 * .0063 - t159 * .0063)) * t9;
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
		    t5 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = t5;
		    t6 = d__1 * d__1;
		    t7 = t6 * t4;
		    t8 = 1 / t7;
		    t9 = sigmaaa * t8;
		    t10 = sqrt(sigmaaa);
		    t11 = t5 * rhoa;
		    t12 = 1 / t11;
		    t13 = t10 * t12;
/* Computing 2nd power */
		    d__1 = t13;
		    t14 = log(t13 + sqrt(d__1 * d__1 + 1));
		    t15 = t13 * t14;
		    t17 = t15 * .021 + 1.;
		    t18 = 1 / t17;
		    t22 = t15 * .0252 + 1.;
		    t23 = 1 / t22;
		    t26 = -.9593273689405774 - t9 * .03640595 * t18 + t9 * 
			    .035481306 * t23;
/* Computing 2nd power */
		    d__1 = rhob;
		    t28 = d__1 * d__1;
		    t29 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = t29;
		    t30 = d__1 * d__1;
		    t31 = t30 * t28;
		    t32 = 1 / t31;
		    t33 = sigmabb * t32;
		    t34 = sqrt(sigmabb);
		    t35 = t29 * rhob;
		    t36 = 1 / t35;
		    t37 = t34 * t36;
/* Computing 2nd power */
		    d__1 = t37;
		    t38 = log(t37 + sqrt(d__1 * d__1 + 1));
		    t39 = t37 * t38;
		    t41 = t39 * .021 + 1.;
		    t42 = 1 / t41;
		    t46 = t39 * .0252 + 1.;
		    t47 = 1 / t46;
		    t50 = -.9593273689405774 - t33 * .03640595 * t42 + t33 * 
			    .035481306 * t47;
		    t52 = pow_dd(&rho, &c_b2);
		    t53 = 1 / t52;
		    t55 = t53 * .3505 + 1.;
		    t56 = 1 / t55;
		    t57 = t56 * rhoa;
		    t58 = 1 / rho;
		    t59 = rhob * t58;
		    t62 = t53 * .25;
		    t63 = exp(-t62);
		    t64 = t63 * t56;
/* Computing 2nd power */
		    d__1 = rho;
		    t65 = d__1 * d__1;
		    t66 = t65 * rho;
/* Computing 2nd power */
		    d__1 = t52;
		    t67 = d__1 * d__1;
		    t69 = 1 / t67 / t66;
		    t70 = rhoa * rhob;
		    t71 = t7 * 36.46239897876478;
		    t72 = t31 * 36.46239897876478;
		    t74 = t53 * t56;
		    t76 = 2.611111111111111 - t53 * .09722222222222222 - t74 *
			     .1363055555555556;
		    t77 = t76 * sigma;
		    t81 = sigmaaa + sigmabb;
		    t83 = (2.5 - t53 * .01388888888888889 - t74 * 
			    .01947222222222222) * 1. * t81;
		    t85 = t62 + t74 * .3505 - 11.;
		    t89 = rhoa * t58 * sigmaaa + t59 * sigmabb;
		    t91 = t85 * .1111111111111111 * t89;
		    t92 = t71 + t72 + t77 - t83 - t91;
		    t96 = t65 * .6666666666666667;
		    t97 = t4 * 1.;
		    t100 = t28 * 1.;
		    t103 = t70 * t92 - t65 * .6666666666666667 * sigma + (t96 
			    - t97) * sigmabb + (t96 - t100) * sigmaaa;
		    zk[i__] = t26 * t11 + t50 * t35 - t57 * .22 * t59 - t64 * 
			    .00869 * t69 * t103;
		    t107 = t4 * rhoa;
		    t109 = 1 / t6 / t107;
		    t110 = sigmaaa * t109;
/* Computing 2nd power */
		    d__1 = t17;
		    t113 = d__1 * d__1;
		    t114 = 1 / t113;
		    t116 = 1 / t5 / t4;
		    t118 = t10 * t116 * t14;
		    t120 = t9 + 1.;
		    t121 = sqrt(t120);
		    t122 = 1 / t121;
		    t123 = t110 * t122;
		    t125 = t118 * -.028 - t123 * .028;
		    t126 = t114 * t125;
/* Computing 2nd power */
		    d__1 = t22;
		    t131 = d__1 * d__1;
		    t132 = 1 / t131;
		    t135 = t118 * -.0336 - t123 * .0336;
		    t136 = t132 * t135;
		    t139 = t110 * .09708253333333333 * t18 + t9 * .03640595 * 
			    t126 - t110 * .094616816 * t23 - t9 * .035481306 *
			     t136;
		    t143 = t56 * rhob;
		    t147 = t6 * rhoa;
		    t149 = t85 * t58;
		    t152 = t147 * 97.23306394337274 - t149 * 
			    .1111111111111111 * sigmaaa;
		    t156 = rhob * t92 + t70 * t152 - rhoa * 2. * sigmabb;
/* Computing 2nd power */
		    d__1 = t55;
		    t160 = d__1 * d__1;
		    t161 = 1 / t160;
		    t162 = t161 * rhoa;
		    t164 = 1 / t52 / t65;
		    t167 = t162 * .02570333333333333 * rhob * t164;
		    t168 = 1 / t65;
		    t169 = rhob * t168;
		    t171 = t57 * .22 * t169;
/* Computing 2nd power */
		    d__1 = t65;
		    t172 = d__1 * d__1;
		    t173 = t172 * rho;
		    t174 = 1 / t173;
		    t175 = t174 * t63;
		    t176 = t56 * t103;
		    t178 = t175 * 7.241666666666667e-4 * t176;
		    t179 = t63 * t161;
		    t182 = t179 * .001015281666666667 * t174 * t103;
		    t184 = 1 / t67 / t172;
		    t187 = t64 * .03186333333333333 * t184 * t103;
		    t189 = 1 / t52 / rho;
		    t191 = t189 * t56;
		    t195 = 1 / t67 / rho * t161;
		    t197 = t189 * .03240740740740741 + t191 * 
			    .04543518518518519 - t195 * .01592503240740741;
		    t208 = t189 * -.08333333333333333 - t191 * 
			    .1168333333333333 + t195 * .04095008333333333;
		    t216 = rhoa * -1. * t168 * sigmaaa - t169 * 1. * sigmabb;
		    t219 = t197 * sigma - (t189 * .00462962962962963 + t191 * 
			    .006490740740740741 - t195 * .00227500462962963) *
			     1. * t81 - t208 * .1111111111111111 * t89 - t85 *
			     .1111111111111111 * t216;
		    t227 = t70 * t219 - rho * 1.333333333333333 * sigma + rho 
			    * 1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa;
		    t230 = t64 * .00869 * t69 * t227;
		    vrhoa[i__] = t139 * t11 + t26 * 1.333333333333333 * t5 - 
			    t143 * .22 * t58 - t64 * .00869 * t69 * t156 - 
			    t167 + t171 - t178 - t182 + t187 - t230;
		    t231 = t28 * rhob;
		    t233 = 1 / t30 / t231;
		    t234 = sigmabb * t233;
/* Computing 2nd power */
		    d__1 = t41;
		    t237 = d__1 * d__1;
		    t238 = 1 / t237;
		    t240 = 1 / t29 / t28;
		    t242 = t34 * t240 * t38;
		    t244 = t33 + 1.;
		    t245 = sqrt(t244);
		    t246 = 1 / t245;
		    t247 = t234 * t246;
		    t249 = t242 * -.028 - t247 * .028;
		    t250 = t238 * t249;
/* Computing 2nd power */
		    d__1 = t46;
		    t255 = d__1 * d__1;
		    t256 = 1 / t255;
		    t259 = t242 * -.0336 - t247 * .0336;
		    t260 = t256 * t259;
		    t263 = t234 * .09708253333333333 * t42 + t33 * .03640595 *
			     t250 - t234 * .094616816 * t47 - t33 * 
			    .035481306 * t260;
		    t270 = t30 * rhob;
		    t274 = t270 * 97.23306394337274 - t149 * 
			    .1111111111111111 * sigmabb;
		    t278 = rhoa * t92 + t70 * t274 - rhob * 2. * sigmaaa;
		    vrhob[i__] = t263 * t35 + t50 * 1.333333333333333 * t29 - 
			    t57 * .22 * t58 - t64 * .00869 * t69 * t278 - 
			    t167 + t171 - t178 - t182 + t187 - t230;
		    t284 = 1 / t10;
		    t286 = t284 * t12 * t14;
		    t288 = t8 * t122;
		    t290 = t286 * .0105 + t288 * .0105;
		    t291 = t114 * t290;
		    t298 = t286 * .0126 + t288 * .0126;
		    t299 = t132 * t298;
		    t302 = t8 * -.03640595 * t18 + t9 * .03640595 * t291 + t8 
			    * .035481306 * t23 - t9 * .035481306 * t299;
		    t304 = t53 * .01388888888888889;
		    t305 = t74 * .01947222222222222;
		    t306 = t85 * rhoa;
		    t309 = t304 - 2.5 + t305 - t306 * .1111111111111111 * t58;
		    t311 = t70 * t309 + t96 - t100;
		    t317 = t70 * t76 - t65 * .6666666666666667;
		    t319 = t64 * t69 * t317;
		    t320 = t319 * .00869;
		    vsigmaaa[i__] = t302 * t11 - t64 * .00869 * t69 * t311 - 
			    t320;
		    vsigmaab[i__] = t319 * -.01738;
		    t324 = 1 / t34;
		    t326 = t324 * t36 * t38;
		    t328 = t32 * t246;
		    t330 = t326 * .0105 + t328 * .0105;
		    t331 = t238 * t330;
		    t338 = t326 * .0126 + t328 * .0126;
		    t339 = t256 * t338;
		    t342 = t32 * -.03640595 * t42 + t33 * .03640595 * t331 + 
			    t32 * .035481306 * t47 - t33 * .035481306 * t339;
		    t344 = t85 * rhob;
		    t347 = t304 - 2.5 + t305 - t344 * .1111111111111111 * t58;
		    t349 = t70 * t347 + t96 - t97;
		    vsigmabb[i__] = t342 * t35 - t64 * .00869 * t69 * t349 - 
			    t320;
/* Computing 2nd power */
		    d__1 = t4;
		    t362 = d__1 * d__1;
		    t365 = sigmaaa / t6 / t362;
		    t371 = 1 / t113 / t17;
/* Computing 2nd power */
		    d__1 = t125;
		    t372 = d__1 * d__1;
		    t379 = t10 / t5 / t107 * t14;
		    t381 = t365 * t122;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t383 = d__1 * d__1;
		    t389 = 1 / t121 / t120;
		    t390 = t383 / t5 / t362 / t107 * t389;
		    t401 = 1 / t131 / t22;
/* Computing 2nd power */
		    d__1 = t135;
		    t402 = d__1 * d__1;
		    t416 = 1 / t160 / t55;
		    t418 = rhob * t69;
		    t420 = t416 * .006006012222222222 * rhoa * t418;
		    t423 = t175 * .001448333333333333 * t56 * t227;
		    t424 = t172 * t65;
		    t426 = 1 / t52 / t424;
		    t427 = t426 * t63;
		    t429 = t427 * 6.034722222222222e-5 * t176;
		    t432 = t427 * 1.692136111111111e-4 * t161 * t103;
		    t435 = t179 * .002030563333333333 * t174 * t227;
		    t439 = t63 * 2.372374827777778e-4 * t416 * t426 * t103;
		    t442 = t64 * .06372666666666667 * t184 * t227;
		    t444 = t164 * t56;
		    t448 = 1 / t67 / t65 * t161;
		    t450 = 1 / t66;
		    t451 = t450 * t416;
		    t474 = rhob * t450;
		    t488 = t64 * .00869 * t69 * (t70 * ((t164 * 
			    -.04320987654320988 - t444 * .06058024691358025 + 
			    t448 * .03185006481481481 - t451 * 
			    .003721149239197531) * sigma - (t164 * 
			    -.00617283950617284 - t444 * .008654320987654321 
			    + t448 * .004550009259259259 - t451 * 
			    5.315927484567901e-4) * 1. * t81 - (t164 * 
			    .1111111111111111 + t444 * .1557777777777778 - 
			    t448 * .08190016666666667 + t451 * 
			    .009568669472222222) * .1111111111111111 * t89 - 
			    t208 * .2222222222222222 * t216 - t85 * 
			    .1111111111111111 * (rhoa * 2. * t450 * sigmaaa + 
			    t474 * 2. * sigmabb)) - sigma * 1.333333333333333 
			    + sigmabb * 1.333333333333333 + sigmaaa * 
			    1.333333333333333);
		    t490 = t208 * t58;
		    t493 = t85 * t168;
		    t500 = t64 * t69 * (rhob * t219 + t70 * (t490 * 
			    -.1111111111111111 * sigmaaa + t493 * 
			    .1111111111111111 * sigmaaa));
		    t502 = t64 * -.00869 * t69 * (rhob * 2. * t152 + t147 * 
			    162.0551065722879 * rhob - sigmabb * 2.) + (t365 *
			     -.3559692888888889 * t18 - t110 * 
			    .1941650666666667 * t126 - t9 * .0728119 * t371 * 
			    t372 + t9 * .03640595 * t114 * (t379 * 
			    .06533333333333333 + t381 * .14 - t390 * 
			    .03733333333333333) + t365 * .3469283253333333 * 
			    t23 + t110 * .189233632 * t136 + t9 * .070962612 *
			     t401 * t402 - t9 * .035481306 * t132 * (t379 * 
			    .0784 + t381 * .168 - t390 * .0448)) * t11 - t420 
			    - t423 - t429 - t432 - t435 - t439 + t442 - t488 
			    - t500 * .01738;
		    t504 = t179 * t174 * t156;
		    t507 = t64 * t184 * t156;
		    t510 = t161 * rhob * t164;
		    t512 = t143 * t168;
		    t515 = t175 * t56 * t156;
		    t522 = 1 / t424;
		    t525 = t179 * .008799107777777778 * t522 * t103;
		    t527 = t57 * .44 * t474;
		    t530 = t522 * .006276111111111111 * t63 * t176;
		    t535 = t162 * .08567777777777778 * rhob / t52 / t66;
		    t540 = t64 * .1486955555555556 / t67 / t173 * t103;
		    t541 = t504 * -.002030563333333333 + t507 * 
			    .06372666666666667 - t510 * .05140666666666667 + 
			    t512 * .44 - t515 * .001448333333333333 + t139 * 
			    2.666666666666667 * t5 + t26 * .4444444444444444 /
			     t6 + t525 - t527 + t530 + t535 - t540;
		    v2rhoa2[i__] = t502 + t541;
		    t542 = -t420 - t423 - t429 - t432 - t435 - t439 + t442 - 
			    t488 + t525 - t527 + t530;
		    t557 = t162 * t164;
		    t559 = t57 * t168;
		    t562 = t175 * t56 * t278;
		    t565 = t179 * t174 * t278;
/* Computing 2nd power */
		    d__1 = t28;
		    t567 = d__1 * d__1;
		    t570 = sigmabb / t30 / t567;
		    t576 = 1 / t237 / t41;
/* Computing 2nd power */
		    d__1 = t249;
		    t577 = d__1 * d__1;
		    t584 = t34 / t29 / t231 * t38;
		    t586 = t570 * t246;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t588 = d__1 * d__1;
		    t594 = 1 / t245 / t244;
		    t595 = t588 / t29 / t567 / t231 * t594;
		    t606 = 1 / t255 / t46;
/* Computing 2nd power */
		    d__1 = t259;
		    t607 = d__1 * d__1;
		    t621 = t64 * t184 * t278;
		    t632 = t64 * t69 * (rhoa * t219 + t70 * (t490 * 
			    -.1111111111111111 * sigmabb + t493 * 
			    .1111111111111111 * sigmabb));
		    t634 = t535 + t263 * 2.666666666666667 * t29 + t50 * 
			    .4444444444444444 / t30 - t64 * .00869 * t69 * (
			    rhoa * 2. * t274 + rhoa * 162.0551065722879 * 
			    t270 - sigmaaa * 2.) - t557 * .05140666666666667 
			    + t559 * .44 - t562 * .001448333333333333 - t565 *
			     .002030563333333333 + (t570 * -.3559692888888889 
			    * t42 - t234 * .1941650666666667 * t250 - t33 * 
			    .0728119 * t576 * t577 + t33 * .03640595 * t238 * 
			    (t584 * .06533333333333333 + t586 * .14 - t595 * 
			    .03733333333333333) + t570 * .3469283253333333 * 
			    t47 + t234 * .189233632 * t260 + t33 * .070962612 
			    * t606 * t607 - t33 * .035481306 * t256 * (t584 * 
			    .0784 + t586 * .168 - t595 * .0448)) * t35 - t540 
			    + t621 * .06372666666666667 - t632 * .01738;
		    v2rhob2[i__] = t542 + t634;
		    t640 = -t420 - t423 - t435 - t439 + t442 - t429 - t432 - 
			    t488 - t500 * .00869 - t504 * .001015281666666667 
			    + t507 * .03186333333333333 - t515 * 
			    7.241666666666667e-4 - t510 * .02570333333333333;
		    t656 = t512 * .22 + t525 - t527 + t530 + t535 - t540 - 
			    t557 * .02570333333333333 + t559 * .22 - t562 * 
			    7.241666666666667e-4 - t565 * .001015281666666667 
			    + t621 * .03186333333333333 - t632 * .00869 - t64 
			    * .00869 * t69 * (t71 + t72 + t77 - t83 - t91 + 
			    rhob * t274 + rhoa * t152) - t56 * .22 * t58;
		    v2rhoab[i__] = t640 + t656;
		    t661 = t8 * t114;
		    t669 = t284 * t116 * t14;
		    t671 = t109 * t122;
		    t677 = sigmaaa / t5 / t362 / t4 * t389;
		    t687 = t8 * t132;
		    t707 = t70 * .1111111111111111 * t149;
		    t714 = t175 * 7.241666666666667e-4 * t56 * t311;
		    t717 = t179 * .001015281666666667 * t174 * t311;
		    t720 = t64 * .03186333333333333 * t184 * t311;
		    t721 = t189 * .00462962962962963;
		    t722 = t191 * .006490740740740741;
		    t723 = t195 * .00227500462962963;
		    t731 = rho * 1.333333333333333;
		    t735 = t64 * .00869 * t69 * (t70 * (-t721 - t722 + t723 - 
			    t208 * .1111111111111111 * rhoa * t58 + t306 * 
			    .1111111111111111 * t168) + t731);
		    t737 = t64 * t418 * t76;
		    t738 = t737 * .00869;
		    t740 = t175 * t56 * t317;
		    t741 = t740 * 7.241666666666667e-4;
		    t743 = t179 * t174 * t317;
		    t744 = t743 * .001015281666666667;
		    t746 = t64 * t184 * t317;
		    t747 = t746 * .03186333333333333;
		    t752 = t64 * t69 * (t70 * t197 - rho * 1.333333333333333);
		    t753 = t752 * .00869;
		    v2rhoasigmaaa[i__] = (t109 * .09708253333333333 * t18 - 
			    t110 * .09708253333333333 * t291 + t661 * 
			    .03640595 * t125 - t9 * .0728119 * t371 * t125 * 
			    t290 + t9 * .03640595 * t114 * (t669 * -.014 - 
			    t671 * .042 + t677 * .014) - t109 * .094616816 * 
			    t23 + t110 * .094616816 * t299 - t687 * 
			    .035481306 * t135 + t9 * .070962612 * t401 * t135 
			    * t298 - t9 * .035481306 * t132 * (t669 * -.0168 
			    - t671 * .0504 + t677 * .0168)) * t11 + t302 * 
			    1.333333333333333 * t5 - t64 * .00869 * t69 * (
			    rhob * t309 - t707) - t714 - t717 + t720 - t735 - 
			    t738 - t741 - t744 + t747 - t753;
		    t755 = t740 * .001448333333333333;
		    t756 = t743 * .002030563333333333;
		    t757 = t746 * .06372666666666667;
		    t758 = t752 * .01738;
		    v2rhoasigmaab[i__] = t737 * -.01738 - t755 - t756 + t757 
			    - t758;
		    t767 = t175 * 7.241666666666667e-4 * t56 * t349;
		    t770 = t179 * .001015281666666667 * t174 * t349;
		    t773 = t64 * .03186333333333333 * t184 * t349;
		    t784 = t64 * .00869 * t69 * (t70 * (-t721 - t722 + t723 - 
			    t208 * .1111111111111111 * rhob * t58 + t344 * 
			    .1111111111111111 * t168) + t731);
		    v2rhoasigmabb[i__] = t64 * -.00869 * t69 * (rhob * t347 - 
			    rhoa * 2.) - t767 - t770 + t773 - t784 - t738 - 
			    t741 - t744 + t747 - t753;
		    t793 = t64 * t69 * rhoa * t76;
		    t794 = t793 * .00869;
		    v2rhobsigmaaa[i__] = t64 * -.00869 * t69 * (rhoa * t309 - 
			    rhob * 2.) - t714 - t717 + t720 - t735 - t794 - 
			    t741 - t744 + t747 - t753;
		    v2rhobsigmaab[i__] = t793 * -.01738 - t755 - t756 + t757 
			    - t758;
		    t800 = t32 * t238;
		    t808 = t324 * t240 * t38;
		    t810 = t233 * t246;
		    t816 = sigmabb / t29 / t567 / t28 * t594;
		    t826 = t32 * t256;
		    v2rhobsigmabb[i__] = (t233 * .09708253333333333 * t42 - 
			    t234 * .09708253333333333 * t331 + t800 * 
			    .03640595 * t249 - t33 * .0728119 * t576 * t249 * 
			    t330 + t33 * .03640595 * t238 * (t808 * -.014 - 
			    t810 * .042 + t816 * .014) - t233 * .094616816 * 
			    t47 + t234 * .094616816 * t339 - t826 * 
			    .035481306 * t259 + t33 * .070962612 * t606 * 
			    t259 * t338 - t33 * .035481306 * t256 * (t808 * 
			    -.0168 - t810 * .0504 + t816 * .0168)) * t35 + 
			    t342 * 1.333333333333333 * t29 - t64 * .00869 * 
			    t69 * (rhoa * t347 - t707) - t767 - t770 + t773 - 
			    t784 - t794 - t741 - t744 + t747 - t753;
/* Computing 2nd power */
		    d__1 = t290;
		    t851 = d__1 * d__1;
		    t858 = 1 / t10 / sigmaaa * t12 * t14;
		    t862 = 1 / sigmaaa * t8 * t122;
		    t867 = 1 / t5 / t362 / rhoa * t389;
/* Computing 2nd power */
		    d__1 = t298;
		    t875 = d__1 * d__1;
		    v2sigmaaa2[i__] = (t661 * .0728119 * t290 - t9 * .0728119 
			    * t371 * t851 + t9 * .03640595 * t114 * (t858 * 
			    -.00525 + t862 * .00525 - t867 * .00525) - t687 * 
			    .070962612 * t298 + t9 * .070962612 * t401 * t875 
			    - t9 * .035481306 * t132 * (t858 * -.0063 + t862 *
			     .0063 - t867 * .0063)) * t11;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t330;
		    t889 = d__1 * d__1;
		    t896 = 1 / t34 / sigmabb * t36 * t38;
		    t900 = 1 / sigmabb * t32 * t246;
		    t905 = 1 / t29 / t567 / rhob * t594;
/* Computing 2nd power */
		    d__1 = t338;
		    t913 = d__1 * d__1;
		    v2sigmabb2[i__] = (t800 * .0728119 * t330 - t33 * 
			    .0728119 * t576 * t889 + t33 * .03640595 * t238 * 
			    (t896 * -.00525 + t900 * .00525 - t905 * .00525) 
			    - t826 * .070962612 * t338 + t33 * .070962612 * 
			    t606 * t913 - t33 * .035481306 * t256 * (t896 * 
			    -.0063 + t900 * .0063 - t905 * .0063)) * t35;
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
} /* uks_xc_edf1__ */

/* Subroutine */ int rks_xc_edf1__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal s1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t21, t13, 
	    t14, t31, t34, t35, t42, t28, t16, t17, t22, t25, t30, t36, t39, 
	    t44, t52, t55, t60, t64, t67, t68, t72, t76, t77, t78, t86, t87, 
	    t37, t40, t45, t50, t54, t70, t75, t80, t81, t90, t91, t94, t102, 
	    t111, t121, t123, t106, t140, t117, t118, t141, t138, t170, t172, 
	    t104, t105, t108, t113, t122, t124, t125, t128, t133, t143, t154, 
	    t158, t162, t168, t174, t175, t182, t183, t186, t194, t205, t207, 
	    t209, t211, t212, t231, t235, t261, t267, t268, t275, t277, t279, 
	    t285, t286, t297, t298, t321, t323, t333, t345, t359, t365, t370, 
	    t378, t380, t383, t393, t458, t465, t469, t473, t481, rho, sigma;


/*     R.D. Adamson, P.M.W. Gill, and J.A. Pople */
/*     Empirical density functionals */
/*     Chem. Phys. Lett. 284 (1998) 6-11 */


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
/* Computing 2nd power */
		d__1 = rho;
		t2 = d__1 * d__1;
		t3 = pow_dd(&rho, &c_b2);
/* Computing 2nd power */
		d__1 = t3;
		t4 = d__1 * d__1;
		t5 = t4 * t2;
		t7 = sigma / t5;
		t8 = sqrt(sigma);
		t9 = t3 * rho;
		t11 = t8 / t9;
/* Computing 2nd power */
		d__1 = t11;
		t13 = log(t11 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		t14 = t11 * t13;
		t28 = 1 / t3;
		t31 = 1 / (t28 * .3505 + 1.);
		t34 = t28 * .25;
		t35 = exp(-t34);
		t42 = t28 * t31;
		zk[i__] = (-.9593273689405774 - t7 * .05779084332790167 / (
			t14 * .02645834204779234 + 1.) + t7 * 
			.05632306246960559 / (t14 * .0317500104573508 + 1.)) *
			 .7937005259840997 * t9 - t31 * .055 * rho - t35 * 
			.00869 * t31 / t4 / t2 / rho * (t2 * .25 * (t5 * 
			11.48493600075277 + (2.611111111111111 - t28 * 
			.09722222222222222 - t42 * .1363055555555556) * sigma 
			- (2.5 - t28 * .01388888888888889 - t42 * 
			.01947222222222222) * .5 * sigma - (t34 + t42 * .3505 
			- 11.) * .02777777777777778 * sigma) - t2 * 
			.4583333333333333 * sigma);
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
/* Computing 2nd power */
		d__1 = rho;
		t2 = d__1 * d__1;
		t3 = pow_dd(&rho, &c_b2);
/* Computing 2nd power */
		d__1 = t3;
		t4 = d__1 * d__1;
		t5 = t4 * t2;
		t6 = 1 / t5;
		t7 = sigma * t6;
		t8 = sqrt(sigma);
		t9 = t3 * rho;
		t10 = 1 / t9;
		t11 = t8 * t10;
/* Computing 2nd power */
		d__1 = t11;
		t13 = log(t11 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		t14 = t11 * t13;
		t16 = t14 * .02645834204779234 + 1.;
		t17 = 1 / t16;
		t21 = t14 * .0317500104573508 + 1.;
		t22 = 1 / t21;
		t25 = -.9593273689405774 - t7 * .05779084332790167 * t17 + t7 
			* .05632306246960559 * t22;
		t28 = 1 / t3;
		t30 = t28 * .3505 + 1.;
		t31 = 1 / t30;
		t34 = t28 * .25;
		t35 = exp(-t34);
		t36 = t35 * t31;
		t39 = 1 / t4 / t2 / rho;
		t42 = t28 * t31;
		t44 = 2.611111111111111 - t28 * .09722222222222222 - t42 * 
			.1363055555555556;
		t52 = t34 + t42 * .3505 - 11.;
		t55 = t5 * 11.48493600075277 + t44 * sigma - (2.5 - t28 * 
			.01388888888888889 - t42 * .01947222222222222) * .5 * 
			sigma - t52 * .02777777777777778 * sigma;
		t60 = t2 * .25 * t55 - t2 * .4583333333333333 * sigma;
		zk[i__] = t25 * .7937005259840997 * t9 - t31 * .055 * rho - 
			t36 * .00869 * t39 * t60;
		t64 = sigma * t39;
/* Computing 2nd power */
		d__1 = t16;
		t67 = d__1 * d__1;
		t68 = 1 / t67;
		t72 = t8 / t3 / t2 * t13;
		t76 = sqrt(t7 * 1.587401051968199 + 1.);
		t77 = 1 / t76;
		t78 = t64 * t77;
/* Computing 2nd power */
		d__1 = t21;
		t86 = d__1 * d__1;
		t87 = 1 / t86;
		t102 = t4 * rho;
		t106 = t52 / rho * sigma;
		t111 = rho * sigma;
/* Computing 2nd power */
		d__1 = t30;
		t117 = d__1 * d__1;
		t118 = 1 / t117;
/* Computing 2nd power */
		d__1 = t2;
		t121 = d__1 * d__1;
		t123 = 1 / t121 / rho;
		t138 = t10 * t31;
		t140 = 1 / t102;
		t141 = t140 * t118;
		s1 = (t64 * .3082178310821422 * t17 + t7 * .05779084332790167 
			* t68 * (t72 * -.0705555787941129 - t78 * 
			.08889445891021917) - t64 * .3003896665045631 * t22 - 
			t7 * .05632306246960559 * t87 * (t72 * 
			-.08466669455293548 - t78 * .106673350692263)) * 
			.3968502629920499 * t9 + t25 * 1.0582673679788 * t3 - 
			t31 * .055 - t36 * .00869 * t39 * (rho * .5 * t55 + 
			t2 * .25 * (t102 * 30.62649600200738 - t106 * 
			.02777777777777778) - t111 * .25);
		vrhoa[i__] = s1 - t118 * .006425833333333333 * t28 - t123 * 
			7.241666666666667e-4 * t35 * t31 * t60 - t35 * 
			.001015281666666667 * t118 * t123 * t60 + t36 * 
			.03186333333333333 / t4 / t121 * t60 - t36 * .00869 * 
			t39 * (t2 * .25 * ((t10 * .03240740740740741 + t138 * 
			.04543518518518519 - t141 * .01592503240740741) * 
			sigma - (t10 * .00462962962962963 + t138 * 
			.006490740740740741 - t141 * .00227500462962963) * .5 
			* sigma - (t10 * -.08333333333333333 - t138 * 
			.1168333333333333 + t141 * .04095008333333333) * 
			.02777777777777778 * sigma + t106 * 
			.02777777777777778) - t111 * .6666666666666667);
		t170 = 1 / t8 * t10 * t13;
		t172 = t6 * t77;
		vsigmaaa[i__] = (t6 * -.2311633733116067 * t17 + t7 * 
			.05779084332790167 * t68 * (t170 * .05291668409558467 
			+ t172 * .06667084418266438) + t6 * .2252922498784224 
			* t22 - t7 * .05632306246960559 * t87 * (t170 * 
			.06350002091470161 + t172 * .08000501301919725)) * 
			.7937005259840997 * t9 + t36 * 9.655555555555556e-4 * 
			t140 - t36 * .03476 * t39 * (t2 * .25 * t44 - t2 * 
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
/* Computing 2nd power */
		d__1 = rho;
		t2 = d__1 * d__1;
		t3 = pow_dd(&rho, &c_b2);
/* Computing 2nd power */
		d__1 = t3;
		t4 = d__1 * d__1;
		t5 = t4 * t2;
		t6 = 1 / t5;
		t7 = sigma * t6;
		t8 = sqrt(sigma);
		t9 = t3 * rho;
		t10 = 1 / t9;
		t11 = t8 * t10;
/* Computing 2nd power */
		d__1 = t11;
		t13 = log(t11 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		t14 = t11 * t13;
		t16 = t14 * .02645834204779234 + 1.;
		t17 = 1 / t16;
		t21 = t14 * .0317500104573508 + 1.;
		t22 = 1 / t21;
		t25 = -.9593273689405774 - t7 * .05779084332790167 * t17 + t7 
			* .05632306246960559 * t22;
		t28 = 1 / t3;
		t30 = t28 * .3505 + 1.;
		t31 = 1 / t30;
		t34 = t28 * .25;
		t35 = exp(-t34);
		t36 = t35 * t31;
		t37 = t2 * rho;
		t39 = 1 / t4 / t37;
		t40 = t5 * 11.48493600075277;
		t42 = t28 * t31;
		t44 = 2.611111111111111 - t28 * .09722222222222222 - t42 * 
			.1363055555555556;
		t45 = t44 * sigma;
		t50 = (2.5 - t28 * .01388888888888889 - t42 * 
			.01947222222222222) * .5 * sigma;
		t52 = t34 + t42 * .3505 - 11.;
		t54 = t52 * .02777777777777778 * sigma;
		t55 = t40 + t45 - t50 - t54;
		t60 = t2 * .25 * t55 - t2 * .4583333333333333 * sigma;
		zk[i__] = t25 * .7937005259840997 * t9 - t31 * .055 * rho - 
			t36 * .00869 * t39 * t60;
		t64 = sigma * t39;
/* Computing 2nd power */
		d__1 = t16;
		t67 = d__1 * d__1;
		t68 = 1 / t67;
		t70 = 1 / t3 / t2;
		t72 = t8 * t70 * t13;
		t75 = t7 * 1.587401051968199 + 1.;
		t76 = sqrt(t75);
		t77 = 1 / t76;
		t78 = t64 * t77;
		t80 = t72 * -.0705555787941129 - t78 * .08889445891021917;
		t81 = t68 * t80;
/* Computing 2nd power */
		d__1 = t21;
		t86 = d__1 * d__1;
		t87 = 1 / t86;
		t90 = t72 * -.08466669455293548 - t78 * .106673350692263;
		t91 = t87 * t90;
		t94 = t64 * .3082178310821422 * t17 + t7 * .05779084332790167 
			* t81 - t64 * .3003896665045631 * t22 - t7 * 
			.05632306246960559 * t91;
		t102 = t4 * rho;
		t104 = 1 / rho;
		t105 = t52 * t104;
		t106 = t105 * sigma;
		t108 = t102 * 30.62649600200738 - t106 * .02777777777777778;
		t111 = rho * sigma;
		t113 = rho * .5 * t55 + t2 * .25 * t108 - t111 * .25;
/* Computing 2nd power */
		d__1 = t30;
		t117 = d__1 * d__1;
		t118 = 1 / t117;
/* Computing 2nd power */
		d__1 = t2;
		t121 = d__1 * d__1;
		t122 = t121 * rho;
		t123 = 1 / t122;
		t124 = t123 * t35;
		t125 = t31 * t60;
		t128 = t35 * t118;
		t133 = 1 / t4 / t121;
		t138 = t10 * t31;
		t140 = 1 / t102;
		t141 = t140 * t118;
		t143 = t10 * .03240740740740741 + t138 * .04543518518518519 - 
			t141 * .01592503240740741;
		t154 = t10 * -.08333333333333333 - t138 * .1168333333333333 + 
			t141 * .04095008333333333;
		t158 = t143 * sigma - (t10 * .00462962962962963 + t138 * 
			.006490740740740741 - t141 * .00227500462962963) * .5 
			* sigma - t154 * .02777777777777778 * sigma + t106 * 
			.02777777777777778;
		t162 = t2 * .25 * t158 - t111 * .6666666666666667;
		vrhoa[i__] = t94 * .3968502629920499 * t9 + t25 * 
			1.0582673679788 * t3 - t31 * .055 - t36 * .00869 * 
			t39 * t113 - t118 * .006425833333333333 * t28 - t124 *
			 7.241666666666667e-4 * t125 - t128 * 
			.001015281666666667 * t123 * t60 + t36 * 
			.03186333333333333 * t133 * t60 - t36 * .00869 * t39 *
			 t162;
		t168 = 1 / t8;
		t170 = t168 * t10 * t13;
		t172 = t6 * t77;
		t174 = t170 * .05291668409558467 + t172 * .06667084418266438;
		t175 = t68 * t174;
		t182 = t170 * .06350002091470161 + t172 * .08000501301919725;
		t183 = t87 * t182;
		t186 = t6 * -.2311633733116067 * t17 + t7 * 
			.05779084332790167 * t175 + t6 * .2252922498784224 * 
			t22 - t7 * .05632306246960559 * t183;
		t194 = t2 * .25 * t44 - t2 * .6666666666666667;
		vsigmaaa[i__] = t186 * .7937005259840997 * t9 + t36 * 
			9.655555555555556e-4 * t140 - t36 * .03476 * t39 * 
			t194;
		t205 = t70 * t31;
		t207 = t6 * t118;
		t209 = 1 / t37;
		t211 = 1 / t117 / t30;
		t212 = t209 * t211;
		t231 = t154 * t104 * sigma;
		t235 = t52 / t2 * sigma;
		t261 = sigma * t133;
		t267 = 1 / t67 / t16;
/* Computing 2nd power */
		d__1 = t80;
		t268 = d__1 * d__1;
		t275 = t8 / t3 / t37 * t13;
		t277 = t261 * t77;
/* Computing 2nd power */
		d__1 = sigma;
		t279 = d__1 * d__1;
		t285 = 1 / t76 / t75;
		t286 = t279 / t3 / t121 / t37 * t285;
		t297 = 1 / t86 / t21;
/* Computing 2nd power */
		d__1 = t90;
		t298 = d__1 * d__1;
		t321 = t121 * t2;
		t323 = 1 / t3 / t321;
		t333 = 1 / t321;
		t345 = rho * t108;
		t359 = t323 * t35;
		t365 = t35 * -4.744749655555556e-4 * t211 * t323 * t60 + t36 *
			 .1274533333333333 * t133 * t113 - t124 * 
			.002896666666666667 * t31 * t113 + t128 * 
			.01759821555555556 * t333 * t60 + t333 * 
			.01255222222222222 * t35 * t125 - t36 * 
			.2973911111111111 / t4 / t122 * t60 - t36 * .00869 * 
			t39 * (t40 + t45 - t50 - t54 + t345) - t36 * .00869 * 
			t39 * (t345 + t5 * 25.52208000167282 - sigma * .5) - 
			t124 * .002896666666666667 * t31 * t162 - t359 * 
			1.206944444444444e-4 * t125 - t359 * 
			3.384272222222222e-4 * t118 * t60;
		s1 = t128 * -.004061126666666667 * t123 * t162 + t36 * 
			.1274533333333333 * t133 * t162 - t36 * .01738 * t39 *
			 (t2 * .25 * ((t70 * -.04320987654320988 - t205 * 
			.06058024691358025 + t207 * .03185006481481481 - t212 
			* .003721149239197531) * sigma - (t70 * 
			-.00617283950617284 - t205 * .008654320987654321 + 
			t207 * .004550009259259259 - t212 * 
			5.315927484567901e-4) * .5 * sigma - (t70 * 
			.1111111111111111 + t205 * .1557777777777778 - t207 * 
			.08190016666666667 + t212 * .009568669472222222) * 
			.02777777777777778 * sigma + t231 * 
			.05555555555555556 - t235 * .05555555555555556) - 
			sigma * .6666666666666667) - t36 * .03476 * t39 * (
			rho * .5 * t158 + t2 * .25 * (t231 * 
			-.02777777777777778 + t235 * .02777777777777778)) - 
			t128 * .004061126666666667 * t123 * t113;
		v2rhoa2[i__] = s1 - t211 * .003003006111111111 * t140 + (t261 
			* -2.260264094602376 * t17 - t64 * .6164356621642845 *
			 t81 - t7 * .1155816866558033 * t267 * t268 + t7 * 
			.05779084332790167 * t68 * (t275 * .3292593677058602 
			+ t277 * .8889445891021917 - t286 * .3762964202352688)
			 + t261 * 2.202857554366796 * t22 + t64 * 
			.6007793330091263 * t91 + t7 * .1126461249392112 * 
			t297 * t298 - t7 * .05632306246960559 * t87 * (t275 * 
			.3951112412470322 + t277 * 1.06673350692263 - t286 * 
			.4515557042823225)) * .3968502629920499 * t9 + t94 * 
			2.116534735957599 * t3 + t25 * .7055115786525331 / t4 
			- t118 * .008567777777777778 * t10 + t365;
		t370 = t6 * t68;
		t378 = t168 * t70 * t13;
		t380 = t39 * t77;
		t383 = sigma * t323 * t285;
		t393 = t6 * t87;
		s1 = (t39 * 1.232871324328569 * t17 - t64 * .3082178310821422 
			* t175 + t370 * .2311633733116067 * t80 - t7 * 
			.1155816866558033 * t267 * t80 * t174 + t7 * 
			.05779084332790167 * t68 * (t378 * -.1411111575882258 
			- t380 * .533366753461315 + t383 * .2822223151764516) 
			- t39 * 1.201558666018253 * t22 + t64 * 
			.3003896665045631 * t183 - t393 * .2252922498784224 * 
			t90 + t7 * .1126461249392112 * t297 * t90 * t182 - t7 
			* .05632306246960559 * t87 * (t378 * 
			-.169333389105871 - t380 * .640040104153578 + t383 * 
			.3386667782117419)) * .3968502629920499 * t9 + t186 * 
			1.0582673679788 * t3 - t36 * .00869 * t39 * (rho * 
			-.9444444444444444 - rho * .02777777777777778 * t52) 
			+ t209 * 8.046296296296296e-5 * t35 * t31 + t128 * 
			1.128090740740741e-4 * t209 + t36 * 
			.01335685185185185 * t6;
		v2rhoasigmaaa[i__] = s1 - t36 * .01738 * t39 * (t2 * .25 * (
			t10 * -1e-23 - t138 * 1e-23 + t105 * 
			.05555555555555556) + rho * 1.333333333333333) - t36 *
			 .01738 * t6 * t44 - t124 * .002896666666666667 * t31 
			* t194 - t128 * .004061126666666667 * t123 * t194 + 
			t36 * .1274533333333333 * t133 * t194 - t36 * .03476 *
			 t39 * (t2 * .25 * t143 - rho * 1.333333333333333);
/* Computing 2nd power */
		d__1 = t174;
		t458 = d__1 * d__1;
		t465 = 1 / t8 / sigma * t10 * t13;
		t469 = 1 / sigma * t6 * t77;
		t473 = 1 / t3 / t122 * t285;
/* Computing 2nd power */
		d__1 = t182;
		t481 = d__1 * d__1;
		v2sigmaaa2[i__] = (t370 * .4623267466232134 * t174 - t7 * 
			.1155816866558033 * t267 * t458 + t7 * 
			.05779084332790167 * t68 * (t465 * -.1058333681911693 
			+ t469 * .1333416883653288 - t473 * .2116667363823387)
			 - t393 * .4505844997568447 * t182 + t7 * 
			.1126461249392112 * t297 * t481 - t7 * 
			.05632306246960559 * t87 * (t465 * -.1270000418294032 
			+ t469 * .1600100260383945 - t473 * .2540000836588064)
			) * .7937005259840997 * t9;
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
} /* rks_xc_edf1__ */


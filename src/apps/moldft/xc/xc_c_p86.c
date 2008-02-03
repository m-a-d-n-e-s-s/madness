/* xc_c_p86.f -- translated by f2c (version 20050501).
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

static doublereal c_b3 = .33333333333333331;
static doublereal c_b4 = .16666666666666666;
static doublereal c_b60 = 0.;


/* ----------------------------------------------------------------------- */

doublereal piecewise_(logical *o, doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val;


/*     Fix for a Maple problem. The Maple Fortran generator is not */
/*     able to transform its piecewise command to Fortran. Therefore */
/*     this function was written to take care of this. */

    if (*o) {
	ret_val = *a;
    } else {
	ret_val = *b;
    }
    return ret_val;
} /* piecewise_ */

/* :C_P86subrstart */
/*    Generated: Tue Mar  9 13:39:23 GMT 2004 */
/* Subroutine */ int uks_c_p86__(integer *ideriv, integer *npt, doublereal *
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
    logical L__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), sqrt(
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t6;
    static logical t7;
    static doublereal t8, t9, t20, t12, t22, t31, t41, t34, t13, t18, t23, 
	    t24, t39, t27, t28, t30, t36, t46, t48, t50, t59, t62, t67, t69, 
	    t74, t75, t76, t79, t80, t81, t84, t29, t33, t40, t43, t44, t47, 
	    t49, t53, t54, t57, t58, t71, t73, t78, t82, t91, t11, t16, t19, 
	    t21, t42, t113, t52, t56, t61, t64, t68, t83, t85, t86, t93, t98, 
	    t100, t103, t108, t112, t114, t115, t118, t119, t120, t122, t123, 
	    t127, t128, t144, t150, t152, t165, t166, t167, t172, t176, t185, 
	    t197, t200, t206, t216, t217, t222, t228, t233, t236, t32, t51, 
	    t60, t77, t88, t96, t99, t101, t102, t105, t109, t110, t124, t129,
	     t130, t133, t134, t138, t139, t141, t158, t163, rho, t173, t175, 
	    t196, t208, t72, t116, t117, t125, t132, t135, t149, t159, t168, 
	    t169, t171, t174, t177, t178, t182, t186, t190, t193, t194, t195, 
	    t198, t202, t204, t207, t210, t213, t214, t224, t225, t230, t231, 
	    t232, t235, t240, t244, t252, t253, t254, t255, t257, t258, t259, 
	    t262, t264, t274, t285, t297, t298, t302, t306, t308, t311, t314, 
	    t315, t317, t319, t320, t321, t324, t325, t328, t333, t334, t335, 
	    t338, t339, t343, t346, t349, t352, t353, t355, t358, t364, t375, 
	    t377, t379, t383, t397, t406, t407, t411, t414, t416, t419, t420, 
	    t434, t437, t439, rhoa, rhob, t442, t444, t450, t453, t478, t484, 
	    t487, t489, t506, t507, t508, t510, t516, t525, t526, t530, t532, 
	    t534, t536, t540, t544, t546, t560, t565, t567, t573, t574, t576, 
	    t578, t580, t590, t591, t595, t596, t599, t600, t603, t604, t605, 
	    t606, t610, t611, t613, t614, t616, t617, t618, t619, t622, t623, 
	    t624, t625, t626, t627, t628, t629, t630, t632, t634, t639, t641, 
	    t645, sigma, sigmaaa, sigmaab, sigmabb;
    extern doublereal piecewise_(logical *, doublereal *, doublereal *);


/*     J.P. Perdew */
/*     Density-functional approximation for the correlation energy of */
/*     the inhomogeneous electron gas */
/*     Phys. Rev. B33 (1986) 8822-8824 */


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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t6 = pow_dd(&t2, &c_b4);
		    t12 = log(t4);
		    L__1 = 1. <= t4;
		    d__1 = -.0843 / (t6 * 1.101176160755631 + 1. + t3 * 
			    .1619735131738333);
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&L__1, &d__1, &d__2);
		    t20 = sqrt(sigmabb);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t31 = (t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6) / (t3 * 5.411317332115466 + 
			    1. + t22 * .1816419932959077 + t2 * 
			    .01763993811759022) + .001667;
		    t34 = pow_dd(&rhob, &c_b4);
		    t39 = exp(t20 * -8.1290825e-4 / t31 / t34 / rhob);
		    t41 = pow_dd(&rhob, &c_b3);
		    zk[i__] = rhob * t18 + t39 * .7937005259840997 * t31 * 
			    sigmabb / t41 / rhob;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t6 = pow_dd(&t2, &c_b4);
		    t12 = log(t4);
		    L__1 = 1. <= t4;
		    d__1 = -.0843 / (t6 * 1.101176160755631 + 1. + t3 * 
			    .1619735131738333);
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&L__1, &d__1, &d__2);
		    t20 = sqrt(sigmaaa);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t31 = (t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6) / (t3 * 5.411317332115466 + 
			    1. + t22 * .1816419932959077 + t2 * 
			    .01763993811759022) + .001667;
		    t34 = pow_dd(&rhoa, &c_b4);
		    t39 = exp(t20 * -8.1290825e-4 / t31 / t34 / rhoa);
		    t41 = pow_dd(&rhoa, &c_b3);
		    zk[i__] = rhoa * t18 + t39 * .7937005259840997 * t31 * 
			    sigmaaa / t41 / rhoa;
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b3);
		    t6 = t5 * .6203504908994;
		    t8 = pow_dd(&t4, &c_b4);
		    t13 = .1423 / (t8 * .8292885914166397 + 1. + t5 * 
			    .20682485366586);
		    t22 = (rhoa - rhob * 1.) * t4;
		    t23 = t22 + 1.;
		    t24 = pow_dd(&t23, &c_b3);
		    t27 = 1. - t22 * 1.;
		    t28 = pow_dd(&t27, &c_b3);
		    t30 = t24 * t23 + t28 * t27 - 2.;
		    t34 = log(t6);
		    t36 = t5 * t34;
		    L__1 = 1. <= t6;
		    d__1 = -t13 + (-.0843 / (t8 * 1.101176160755631 + 1. + t5 
			    * .1619735131738333) + t13) * 1.923661050931536 * 
			    t30;
		    d__2 = t34 * .0311 - .048 + t36 * .0012407009817988 - t5 *
			     .00719606569443304 + (t34 * -.01555 + .0211 - 
			    t36 * 8.0645563816922e-4 + t5 * 
			    .00421838333811592) * 1.923661050931536 * t30;
		    t46 = piecewise_(&L__1, &d__1, &d__2);
		    t48 = sqrt(sigma);
/* Computing 2nd power */
		    d__1 = t5;
		    t50 = d__1 * d__1;
		    t59 = (t5 * .01443307452126544 + .002568 + t50 * 
			    2.843543831490386e-6) / (t5 * 5.411317332115466 + 
			    1. + t50 * .1816419932959077 + t4 * 
			    .01763993811759022) + .001667;
		    t62 = pow_dd(&rho, &c_b4);
		    t67 = exp(t48 * -8.1290825e-4 / t59 / t62 / rho);
		    t69 = pow_dd(&rho, &c_b3);
		    t74 = t22 * .5 + .5;
		    t75 = pow_dd(&t74, &c_b3);
/* Computing 2nd power */
		    d__1 = t75;
		    t76 = d__1 * d__1;
		    t79 = .5 - t22 * .5;
		    t80 = pow_dd(&t79, &c_b3);
/* Computing 2nd power */
		    d__1 = t80;
		    t81 = d__1 * d__1;
		    t84 = sqrt(t76 * t74 + t81 * t79);
		    zk[i__] = rho * t46 + t67 * .7937005259840997 * t59 * 
			    sigma / t69 / rho / t84;
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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t7 = 1. <= t4;
		    t6 = pow_dd(&t2, &c_b4);
		    t9 = t6 * 1.101176160755631 + 1. + t3 * .1619735131738333;
		    t12 = log(t4);
		    d__1 = -.0843 / t9;
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&t7, &d__1, &d__2);
		    t20 = sqrt(sigmabb);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t24 = t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6;
		    t28 = t3 * 5.411317332115466 + 1. + t22 * 
			    .1816419932959077 + t2 * .01763993811759022;
		    t29 = 1 / t28;
		    t31 = t24 * t29 + .001667;
		    t33 = t20 / t31;
		    t34 = pow_dd(&rhob, &c_b4);
		    t36 = 1 / t34 / rhob;
		    t39 = exp(t33 * -8.1290825e-4 * t36);
		    t40 = t39 * t31;
		    t41 = pow_dd(&rhob, &c_b3);
		    t43 = 1 / t41 / rhob;
		    t44 = sigmabb * t43;
		    zk[i__] = rhob * t18 + t40 * .7937005259840997 * t44;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t9;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t49 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t50 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t57 = 1 / t22;
		    t58 = t57 * t54;
		    d__1 = .0843 / t47 * (-.1835293601259385 / t50 / t6 * t54 
			    - t58 * .05399117105794445);
		    d__2 = t2 * -.005183333333333333 - t57 * 
			    1.447484478765267e-4 * t12 * t54 - t3 * 
			    1.447484478765267e-4 * t2 + t58 * 
			    9.9256078543904e-4;
		    t71 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t31;
		    t73 = d__1 * d__1;
		    t78 = 1 / t3 * t54;
/* Computing 2nd power */
		    d__1 = t28;
		    t82 = d__1 * d__1;
		    t91 = (t58 * -.004811024840421814 - t78 * 
			    1.895695887660258e-6) * t29 - t24 * 1. / t82 * (
			    t58 * -1.803772444038489 - t78 * 
			    .1210946621972718 - t54 * .01763993811759022);
		    vrhob[i__] = t18 + rhob * t71 + (t20 * 8.1290825e-4 / t73 
			    * t36 * t91 + t33 * 9.483929583333333e-4 / t34 / 
			    t53) * .7937005259840997 * t39 * t31 * sigmabb * 
			    t43 + t39 * .7937005259840997 * t91 * t44 - t40 * 
			    1.0582673679788 * sigmabb / t41 / t53;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t113 = sqrt(rhob);
		    vsigmabb[i__] = t20 * -3.22602852800907e-4 / t113 / t53 * 
			    t39 + t40 * .7937005259840997 * t43;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t7 = 1. <= t4;
		    t6 = pow_dd(&t2, &c_b4);
		    t9 = t6 * 1.101176160755631 + 1. + t3 * .1619735131738333;
		    t12 = log(t4);
		    d__1 = -.0843 / t9;
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&t7, &d__1, &d__2);
		    t20 = sqrt(sigmaaa);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t24 = t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6;
		    t28 = t3 * 5.411317332115466 + 1. + t22 * 
			    .1816419932959077 + t2 * .01763993811759022;
		    t29 = 1 / t28;
		    t31 = t24 * t29 + .001667;
		    t33 = t20 / t31;
		    t34 = pow_dd(&rhoa, &c_b4);
		    t36 = 1 / t34 / rhoa;
		    t39 = exp(t33 * -8.1290825e-4 * t36);
		    t40 = t39 * t31;
		    t41 = pow_dd(&rhoa, &c_b3);
		    t43 = 1 / t41 / rhoa;
		    t44 = sigmaaa * t43;
		    zk[i__] = rhoa * t18 + t40 * .7937005259840997 * t44;
/* Computing 2nd power */
		    d__1 = t9;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t49 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t50 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t57 = 1 / t22;
		    t58 = t57 * t54;
		    d__1 = .0843 / t47 * (-.1835293601259385 / t50 / t6 * t54 
			    - t58 * .05399117105794445);
		    d__2 = t2 * -.005183333333333333 - t57 * 
			    1.447484478765267e-4 * t12 * t54 - t3 * 
			    1.447484478765267e-4 * t2 + t58 * 
			    9.9256078543904e-4;
		    t71 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t31;
		    t73 = d__1 * d__1;
		    t78 = 1 / t3 * t54;
/* Computing 2nd power */
		    d__1 = t28;
		    t82 = d__1 * d__1;
		    t91 = (t58 * -.004811024840421814 - t78 * 
			    1.895695887660258e-6) * t29 - t24 * 1. / t82 * (
			    t58 * -1.803772444038489 - t78 * 
			    .1210946621972718 - t54 * .01763993811759022);
		    vrhoa[i__] = t18 + rhoa * t71 + (t20 * 8.1290825e-4 / t73 
			    * t36 * t91 + t33 * 9.483929583333333e-4 / t34 / 
			    t53) * .7937005259840997 * t39 * t31 * sigmaaa * 
			    t43 + t39 * .7937005259840997 * t91 * t44 - t40 * 
			    1.0582673679788 * sigmaaa / t41 / t53;
		    vrhob[i__] = 0.;
		    t113 = sqrt(rhoa);
		    vsigmaaa[i__] = t20 * -3.22602852800907e-4 / t113 / t53 * 
			    t39 + t40 * .7937005259840997 * t43;
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b3);
		    t6 = t5 * .6203504908994;
		    t7 = 1. <= t6;
		    t8 = pow_dd(&t4, &c_b4);
		    t11 = t8 * .8292885914166397 + 1. + t5 * .20682485366586;
		    t13 = .1423 / t11;
		    t16 = t8 * 1.101176160755631 + 1. + t5 * 
			    .1619735131738333;
		    t19 = -.0843 / t16 + t13;
		    t21 = rhoa - rhob * 1.;
		    t22 = t21 * t4;
		    t23 = t22 + 1.;
		    t24 = pow_dd(&t23, &c_b3);
		    t27 = 1. - t22 * 1.;
		    t28 = pow_dd(&t27, &c_b3);
		    t30 = t24 * t23 + t28 * t27 - 2.;
		    t34 = log(t6);
		    t36 = t5 * t34;
		    t42 = t34 * -.01555 + .0211 - t36 * 8.0645563816922e-4 + 
			    t5 * .00421838333811592;
		    d__1 = -t13 + t19 * 1.923661050931536 * t30;
		    d__2 = t34 * .0311 - .048 + t36 * .0012407009817988 - t5 *
			     .00719606569443304 + t42 * 1.923661050931536 * 
			    t30;
		    t46 = piecewise_(&t7, &d__1, &d__2);
		    t48 = sqrt(sigma);
/* Computing 2nd power */
		    d__1 = t5;
		    t50 = d__1 * d__1;
		    t52 = t5 * .01443307452126544 + .002568 + t50 * 
			    2.843543831490386e-6;
		    t56 = t5 * 5.411317332115466 + 1. + t50 * 
			    .1816419932959077 + t4 * .01763993811759022;
		    t57 = 1 / t56;
		    t59 = t52 * t57 + .001667;
		    t61 = t48 / t59;
		    t62 = pow_dd(&rho, &c_b4);
		    t64 = 1 / t62 / rho;
		    t67 = exp(t61 * -8.1290825e-4 * t64);
		    t68 = t67 * t59;
		    t69 = pow_dd(&rho, &c_b3);
		    t71 = 1 / t69 / rho;
		    t74 = t22 * .5 + .5;
		    t75 = pow_dd(&t74, &c_b3);
/* Computing 2nd power */
		    d__1 = t75;
		    t76 = d__1 * d__1;
		    t79 = .5 - t22 * .5;
		    t80 = pow_dd(&t79, &c_b3);
/* Computing 2nd power */
		    d__1 = t80;
		    t81 = d__1 * d__1;
		    t83 = t76 * t74 + t81 * t79;
		    t84 = sqrt(t83);
		    t85 = 1 / t84;
		    t86 = sigma * t71 * t85;
		    zk[i__] = rho * t46 + t68 * .7937005259840997 * t86;
		    t93 = t24 * 1.333333333333333 * t4 - t28 * 
			    1.333333333333333 * t4;
		    d__1 = t19 * 1.923661050931536 * t93;
		    d__2 = t42 * 1.923661050931536 * t93;
		    t98 = piecewise_(&t7, &d__1, &d__2);
		    t100 = t68 * sigma;
		    t103 = t71 / t84 / t83;
		    t108 = t76 * .8333333333333333 * t4 - t81 * 
			    .8333333333333333 * t4;
/* Computing 2nd power */
		    d__1 = t11;
		    t112 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t8;
		    t114 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t114;
		    t115 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rho;
		    t118 = d__1 * d__1;
		    t119 = 1 / t118;
		    t120 = 1 / t115 / t8 * t119;
		    t122 = 1 / t50;
		    t123 = t122 * t119;
		    t127 = .1423 / t112 * (t120 * -.1382147652361066 - t123 * 
			    .06894161788861999);
/* Computing 2nd power */
		    d__1 = t16;
		    t128 = d__1 * d__1;
		    t144 = t24 * -1.333333333333333 * t21 * t119 + t28 * 
			    1.333333333333333 * t21 * t119;
		    t150 = t122 * t34 * t119;
		    t152 = t5 * t4;
		    d__1 = t127 + (.0843 / t128 * (t120 * -.1835293601259385 
			    - t123 * .05399117105794445) - t127) * 
			    1.923661050931536 * t30 + t19 * 1.923661050931536 
			    * t144;
		    d__2 = t4 * -.01036666666666667 - t150 * 
			    4.135669939329333e-4 - t152 * 
			    4.135669939329333e-4 + t123 * .002398688564811013 
			    + (t4 * .005183333333333333 + t150 * 
			    2.688185460564067e-4 + t152 * 
			    2.688185460564067e-4 - t123 * .001406127779371973)
			     * 1.923661050931536 * t30 + t42 * 
			    1.923661050931536 * t144;
		    t165 = piecewise_(&t7, &d__1, &d__2);
		    t166 = rho * t165;
/* Computing 2nd power */
		    d__1 = t59;
		    t167 = d__1 * d__1;
		    t172 = 1 / t5 * t119;
/* Computing 2nd power */
		    d__1 = t56;
		    t176 = d__1 * d__1;
		    t185 = (t123 * -.004811024840421814 - t172 * 
			    1.895695887660258e-6) * t57 - t52 * 1. / t176 * (
			    t123 * -1.803772444038489 - t172 * 
			    .1210946621972718 - t119 * .01763993811759022);
		    t197 = (t48 * 8.1290825e-4 / t167 * t64 * t185 + t61 * 
			    9.483929583333333e-4 / t62 / t118) * 
			    .7937005259840997 * t67 * t59 * t86;
		    t200 = t67 * .7937005259840997 * t185 * t86;
		    t206 = t68 * 1.0582673679788 * sigma / t69 / t118 * t85;
		    t216 = t100 * .3968502629920499 * t103 * (t76 * 
			    -.8333333333333333 * t21 * t119 + t81 * 
			    .8333333333333333 * t21 * t119);
		    vrhoa[i__] = rho * t98 - t100 * .3968502629920499 * t103 *
			     t108 + t46 + t166 + t197 + t200 - t206 - t216;
		    t217 = -t93;
		    d__1 = t19 * 1.923661050931536 * t217;
		    d__2 = t42 * 1.923661050931536 * t217;
		    t222 = piecewise_(&t7, &d__1, &d__2);
		    vrhob[i__] = rho * t222 + t100 * .3968502629920499 * t103 
			    * t108 + t46 + t166 + t197 + t200 - t206 - t216;
		    t228 = sqrt(rho);
		    t233 = t48 / t228 / t118 * t67 * t85;
		    t236 = t68 * t71 * t85;
		    vsigmaaa[i__] = t233 * -3.22602852800907e-4 + t236 * 
			    .7937005259840997;
		    vsigmaab[i__] = t233 * -6.45205705601814e-4 + t236 * 
			    1.587401051968199;
		    vsigmabb[i__] = vsigmaaa[i__];
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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t7 = 1. <= t4;
		    t6 = pow_dd(&t2, &c_b4);
		    t9 = t6 * 1.101176160755631 + 1. + t3 * .1619735131738333;
		    t12 = log(t4);
		    d__1 = -.0843 / t9;
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&t7, &d__1, &d__2);
		    t20 = sqrt(sigmabb);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t24 = t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6;
		    t28 = t3 * 5.411317332115466 + 1. + t22 * 
			    .1816419932959077 + t2 * .01763993811759022;
		    t29 = 1 / t28;
		    t31 = t24 * t29 + .001667;
		    t32 = 1 / t31;
		    t33 = t20 * t32;
		    t34 = pow_dd(&rhob, &c_b4);
		    t36 = 1 / t34 / rhob;
		    t39 = exp(t33 * -8.1290825e-4 * t36);
		    t40 = t39 * t31;
		    t41 = pow_dd(&rhob, &c_b3);
		    t43 = 1 / t41 / rhob;
		    t44 = sigmabb * t43;
		    zk[i__] = rhob * t18 + t40 * .7937005259840997 * t44;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t9;
		    t47 = d__1 * d__1;
		    t48 = 1 / t47;
/* Computing 2nd power */
		    d__1 = t6;
		    t49 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t50 = d__1 * d__1;
		    t51 = t50 * t6;
		    t52 = 1 / t51;
/* Computing 2nd power */
		    d__1 = rhob;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t57 = 1 / t22;
		    t58 = t57 * t54;
		    t60 = t52 * -.1835293601259385 * t54 - t58 * 
			    .05399117105794445;
		    t64 = t57 * t12;
		    t67 = t3 * t2;
		    d__1 = t48 * .0843 * t60;
		    d__2 = t2 * -.005183333333333333 - t64 * 
			    1.447484478765267e-4 * t54 - t67 * 
			    1.447484478765267e-4 + t58 * 9.9256078543904e-4;
		    t71 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t31;
		    t73 = d__1 * d__1;
		    t75 = t20 / t73;
		    t77 = 1 / t3;
		    t78 = t77 * t54;
		    t80 = t58 * -.004811024840421814 - t78 * 
			    1.895695887660258e-6;
/* Computing 2nd power */
		    d__1 = t28;
		    t82 = d__1 * d__1;
		    t83 = 1 / t82;
		    t84 = t24 * t83;
		    t88 = t58 * -1.803772444038489 - t78 * .1210946621972718 
			    - t54 * .01763993811759022;
		    t91 = t80 * t29 - t84 * 1. * t88;
		    t96 = 1 / t34 / t53;
		    t99 = t75 * 8.1290825e-4 * t36 * t91 + t33 * 
			    9.483929583333333e-4 * t96;
		    t100 = t99 * t39;
		    t101 = t31 * sigmabb;
		    t102 = t101 * t43;
		    t105 = t39 * t91;
		    t109 = 1 / t41 / t53;
		    t110 = sigmabb * t109;
		    vrhob[i__] = t18 + rhob * t71 + t100 * .7937005259840997 *
			     t102 + t105 * .7937005259840997 * t44 - t40 * 
			    1.0582673679788 * t110;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t113 = sqrt(rhob);
		    t115 = 1 / t113 / t53;
		    vsigmabb[i__] = t20 * -3.22602852800907e-4 * t115 * t39 + 
			    t40 * .7937005259840997 * t43;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t60;
		    t124 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t53;
		    t129 = d__1 * d__1;
		    t130 = 1 / t129;
		    t133 = t53 * rhob;
		    t134 = 1 / t133;
		    t138 = 1 / t22 / t2;
		    t139 = t138 * t130;
		    t141 = t57 * t134;
		    d__1 = -.1686 / t47 / t9 * t124 + t48 * .0843 * (
			    -.1529411334382821 / t51 / t2 * t130 + t52 * 
			    .367058720251877 * t134 - t139 * 
			    .03599411403862963 + t141 * .1079823421158889);
		    d__2 = t54 * .005183333333333333 - t138 * 
			    9.649896525101778e-5 * t12 * t130 - t141 * 
			    .001888622605627062 + t64 * 2.894968957530533e-4 *
			     t134 + t3 * 1.447484478765267e-4 * t54 + t139 * 
			    6.617071902926934e-4;
		    t158 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t91;
		    t163 = d__1 * d__1;
		    t173 = 1 / t67 * t130;
		    t175 = t77 * t134;
/* Computing 2nd power */
		    d__1 = t88;
		    t185 = d__1 * d__1;
		    t196 = (t139 * -.003207349893614542 + t141 * 
			    .009622049680843627 - t173 * 6.318986292200858e-7 
			    + t175 * 3.791391775320515e-6) * t29 - t80 * 2. * 
			    t83 * t88 + t24 * 2. / t82 / t28 * t185 - t84 * 
			    1. * (t139 * -1.202514962692326 + t141 * 
			    3.607544888076978 - t173 * .04036488739909061 + 
			    t175 * .2421893243945437 + t134 * 
			    .03527987623518044);
/* Computing 2nd power */
		    d__1 = t99;
		    t208 = d__1 * d__1;
		    v2rhob2[i__] = t71 * 2. + rhob * t158 + (t20 * 
			    -.0016258165 / t73 / t31 * t36 * t163 - t75 * 
			    .001896785916666667 * t96 * t91 + t75 * 
			    8.1290825e-4 * t36 * t196 - t33 * 
			    .002054851409722222 / t34 / t133) * 
			    .7937005259840997 * t39 * t102 + t208 * 
			    .7937005259840997 * t39 * t102 + t100 * 
			    1.587401051968199 * t91 * sigmabb * t43 - t100 * 
			    2.116534735957599 * t101 * t109 + t39 * 
			    .7937005259840997 * t196 * t44 - t105 * 
			    2.116534735957599 * t110 + t40 * 
			    2.469290525283866 * sigmabb / t41 / t133;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t41;
		    t233 = d__1 * d__1;
		    v2sigmabb2[i__] = -4.839042792013605e-4 / t20 * t115 * 
			    t39 + 1.311232602576965e-7 / t233 / t133 * t32 * 
			    t39;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b3);
		    t4 = t3 * .6203504908994;
		    t7 = 1. <= t4;
		    t6 = pow_dd(&t2, &c_b4);
		    t9 = t6 * 1.101176160755631 + 1. + t3 * .1619735131738333;
		    t12 = log(t4);
		    d__1 = -.0843 / t9;
		    d__2 = t12 * .01555 - .0269 + t3 * 4.3424534362958e-4 * 
			    t12 - t3 * .00297768235631712;
		    t18 = piecewise_(&t7, &d__1, &d__2);
		    t20 = sqrt(sigmaaa);
/* Computing 2nd power */
		    d__1 = t3;
		    t22 = d__1 * d__1;
		    t24 = t3 * .01443307452126544 + .002568 + t22 * 
			    2.843543831490386e-6;
		    t28 = t3 * 5.411317332115466 + 1. + t22 * 
			    .1816419932959077 + t2 * .01763993811759022;
		    t29 = 1 / t28;
		    t31 = t24 * t29 + .001667;
		    t32 = 1 / t31;
		    t33 = t20 * t32;
		    t34 = pow_dd(&rhoa, &c_b4);
		    t36 = 1 / t34 / rhoa;
		    t39 = exp(t33 * -8.1290825e-4 * t36);
		    t40 = t39 * t31;
		    t41 = pow_dd(&rhoa, &c_b3);
		    t43 = 1 / t41 / rhoa;
		    t44 = sigmaaa * t43;
		    zk[i__] = rhoa * t18 + t40 * .7937005259840997 * t44;
/* Computing 2nd power */
		    d__1 = t9;
		    t47 = d__1 * d__1;
		    t48 = 1 / t47;
/* Computing 2nd power */
		    d__1 = t6;
		    t49 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t50 = d__1 * d__1;
		    t51 = t50 * t6;
		    t52 = 1 / t51;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t57 = 1 / t22;
		    t58 = t57 * t54;
		    t60 = t52 * -.1835293601259385 * t54 - t58 * 
			    .05399117105794445;
		    t64 = t57 * t12;
		    t67 = t3 * t2;
		    d__1 = t48 * .0843 * t60;
		    d__2 = t2 * -.005183333333333333 - t64 * 
			    1.447484478765267e-4 * t54 - t67 * 
			    1.447484478765267e-4 + t58 * 9.9256078543904e-4;
		    t71 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t31;
		    t73 = d__1 * d__1;
		    t75 = t20 / t73;
		    t77 = 1 / t3;
		    t78 = t77 * t54;
		    t80 = t58 * -.004811024840421814 - t78 * 
			    1.895695887660258e-6;
/* Computing 2nd power */
		    d__1 = t28;
		    t82 = d__1 * d__1;
		    t83 = 1 / t82;
		    t84 = t24 * t83;
		    t88 = t58 * -1.803772444038489 - t78 * .1210946621972718 
			    - t54 * .01763993811759022;
		    t91 = t80 * t29 - t84 * 1. * t88;
		    t96 = 1 / t34 / t53;
		    t99 = t75 * 8.1290825e-4 * t36 * t91 + t33 * 
			    9.483929583333333e-4 * t96;
		    t100 = t99 * t39;
		    t101 = t31 * sigmaaa;
		    t102 = t101 * t43;
		    t105 = t39 * t91;
		    t109 = 1 / t41 / t53;
		    t110 = sigmaaa * t109;
		    vrhoa[i__] = t18 + rhoa * t71 + t100 * .7937005259840997 *
			     t102 + t105 * .7937005259840997 * t44 - t40 * 
			    1.0582673679788 * t110;
		    vrhob[i__] = 0.;
		    t113 = sqrt(rhoa);
		    t115 = 1 / t113 / t53;
		    vsigmaaa[i__] = t20 * -3.22602852800907e-4 * t115 * t39 + 
			    t40 * .7937005259840997 * t43;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t60;
		    t124 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t53;
		    t129 = d__1 * d__1;
		    t130 = 1 / t129;
		    t133 = t53 * rhoa;
		    t134 = 1 / t133;
		    t138 = 1 / t22 / t2;
		    t139 = t138 * t130;
		    t141 = t57 * t134;
		    d__1 = -.1686 / t47 / t9 * t124 + t48 * .0843 * (
			    -.1529411334382821 / t51 / t2 * t130 + t52 * 
			    .367058720251877 * t134 - t139 * 
			    .03599411403862963 + t141 * .1079823421158889);
		    d__2 = t54 * .005183333333333333 - t138 * 
			    9.649896525101778e-5 * t12 * t130 - t141 * 
			    .001888622605627062 + t64 * 2.894968957530533e-4 *
			     t134 + t3 * 1.447484478765267e-4 * t54 + t139 * 
			    6.617071902926934e-4;
		    t158 = piecewise_(&t7, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t91;
		    t163 = d__1 * d__1;
		    t173 = 1 / t67 * t130;
		    t175 = t77 * t134;
/* Computing 2nd power */
		    d__1 = t88;
		    t185 = d__1 * d__1;
		    t196 = (t139 * -.003207349893614542 + t141 * 
			    .009622049680843627 - t173 * 6.318986292200858e-7 
			    + t175 * 3.791391775320515e-6) * t29 - t80 * 2. * 
			    t83 * t88 + t24 * 2. / t82 / t28 * t185 - t84 * 
			    1. * (t139 * -1.202514962692326 + t141 * 
			    3.607544888076978 - t173 * .04036488739909061 + 
			    t175 * .2421893243945437 + t134 * 
			    .03527987623518044);
/* Computing 2nd power */
		    d__1 = t99;
		    t208 = d__1 * d__1;
		    v2rhoa2[i__] = t71 * 2. + rhoa * t158 + (t20 * 
			    -.0016258165 / t73 / t31 * t36 * t163 - t75 * 
			    .001896785916666667 * t96 * t91 + t75 * 
			    8.1290825e-4 * t36 * t196 - t33 * 
			    .002054851409722222 / t34 / t133) * 
			    .7937005259840997 * t39 * t102 + t208 * 
			    .7937005259840997 * t39 * t102 + t100 * 
			    1.587401051968199 * t91 * sigmaaa * t43 - t100 * 
			    2.116534735957599 * t101 * t109 + t39 * 
			    .7937005259840997 * t196 * t44 - t105 * 
			    2.116534735957599 * t110 + t40 * 
			    2.469290525283866 * sigmaaa / t41 / t133;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t41;
		    t233 = d__1 * d__1;
		    v2sigmaaa2[i__] = -4.839042792013605e-4 / t20 * t115 * 
			    t39 + 1.311232602576965e-7 / t233 / t133 * t32 * 
			    t39;
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b3);
		    t6 = t5 * .6203504908994;
		    t7 = 1. <= t6;
		    t8 = pow_dd(&t4, &c_b4);
		    t11 = t8 * .8292885914166397 + 1. + t5 * .20682485366586;
		    t13 = .1423 / t11;
		    t16 = t8 * 1.101176160755631 + 1. + t5 * 
			    .1619735131738333;
		    t19 = -.0843 / t16 + t13;
		    t21 = rhoa - rhob * 1.;
		    t22 = t21 * t4;
		    t23 = t22 + 1.;
		    t24 = pow_dd(&t23, &c_b3);
		    t27 = 1. - t22 * 1.;
		    t28 = pow_dd(&t27, &c_b3);
		    t30 = t24 * t23 + t28 * t27 - 2.;
		    t34 = log(t6);
		    t36 = t5 * t34;
		    t42 = t34 * -.01555 + .0211 - t36 * 8.0645563816922e-4 + 
			    t5 * .00421838333811592;
		    d__1 = -t13 + t19 * 1.923661050931536 * t30;
		    d__2 = t34 * .0311 - .048 + t36 * .0012407009817988 - t5 *
			     .00719606569443304 + t42 * 1.923661050931536 * 
			    t30;
		    t46 = piecewise_(&t7, &d__1, &d__2);
		    t48 = sqrt(sigma);
/* Computing 2nd power */
		    d__1 = t5;
		    t50 = d__1 * d__1;
		    t52 = t5 * .01443307452126544 + .002568 + t50 * 
			    2.843543831490386e-6;
		    t56 = t5 * 5.411317332115466 + 1. + t50 * 
			    .1816419932959077 + t4 * .01763993811759022;
		    t57 = 1 / t56;
		    t59 = t52 * t57 + .001667;
		    t60 = 1 / t59;
		    t61 = t48 * t60;
		    t62 = pow_dd(&rho, &c_b4);
		    t64 = 1 / t62 / rho;
		    t67 = exp(t61 * -8.1290825e-4 * t64);
		    t68 = t67 * t59;
		    t69 = pow_dd(&rho, &c_b3);
		    t71 = 1 / t69 / rho;
		    t72 = sigma * t71;
		    t74 = t22 * .5 + .5;
		    t75 = pow_dd(&t74, &c_b3);
/* Computing 2nd power */
		    d__1 = t75;
		    t76 = d__1 * d__1;
		    t79 = .5 - t22 * .5;
		    t80 = pow_dd(&t79, &c_b3);
/* Computing 2nd power */
		    d__1 = t80;
		    t81 = d__1 * d__1;
		    t83 = t76 * t74 + t81 * t79;
		    t84 = sqrt(t83);
		    t85 = 1 / t84;
		    t86 = t72 * t85;
		    zk[i__] = rho * t46 + t68 * .7937005259840997 * t86;
		    t93 = t24 * 1.333333333333333 * t4 - t28 * 
			    1.333333333333333 * t4;
		    d__1 = t19 * 1.923661050931536 * t93;
		    d__2 = t42 * 1.923661050931536 * t93;
		    t98 = piecewise_(&t7, &d__1, &d__2);
		    t100 = t68 * sigma;
		    t102 = 1 / t84 / t83;
		    t103 = t71 * t102;
		    t108 = t76 * .8333333333333333 * t4 - t81 * 
			    .8333333333333333 * t4;
		    t109 = t103 * t108;
/* Computing 2nd power */
		    d__1 = t11;
		    t112 = d__1 * d__1;
		    t113 = 1 / t112;
/* Computing 2nd power */
		    d__1 = t8;
		    t114 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t114;
		    t115 = d__1 * d__1;
		    t116 = t115 * t8;
		    t117 = 1 / t116;
/* Computing 2nd power */
		    d__1 = rho;
		    t118 = d__1 * d__1;
		    t119 = 1 / t118;
		    t120 = t117 * t119;
		    t122 = 1 / t50;
		    t123 = t122 * t119;
		    t125 = t120 * -.1382147652361066 - t123 * 
			    .06894161788861999;
		    t127 = t113 * .1423 * t125;
/* Computing 2nd power */
		    d__1 = t16;
		    t128 = d__1 * d__1;
		    t129 = 1 / t128;
		    t132 = t120 * -.1835293601259385 - t123 * 
			    .05399117105794445;
		    t135 = t129 * .0843 * t132 - t127;
		    t138 = t24 * t21;
		    t141 = t28 * t21;
		    t144 = t138 * -1.333333333333333 * t119 + t141 * 
			    1.333333333333333 * t119;
		    t149 = t122 * t34;
		    t150 = t149 * t119;
		    t152 = t5 * t4;
		    t159 = t4 * .005183333333333333 + t150 * 
			    2.688185460564067e-4 + t152 * 
			    2.688185460564067e-4 - t123 * .001406127779371973;
		    d__1 = t127 + t135 * 1.923661050931536 * t30 + t19 * 
			    1.923661050931536 * t144;
		    d__2 = t4 * -.01036666666666667 - t150 * 
			    4.135669939329333e-4 - t152 * 
			    4.135669939329333e-4 + t123 * .002398688564811013 
			    + t159 * 1.923661050931536 * t30 + t42 * 
			    1.923661050931536 * t144;
		    t165 = piecewise_(&t7, &d__1, &d__2);
		    t166 = rho * t165;
/* Computing 2nd power */
		    d__1 = t59;
		    t167 = d__1 * d__1;
		    t168 = 1 / t167;
		    t169 = t48 * t168;
		    t171 = 1 / t5;
		    t172 = t171 * t119;
		    t174 = t123 * -.004811024840421814 - t172 * 
			    1.895695887660258e-6;
/* Computing 2nd power */
		    d__1 = t56;
		    t176 = d__1 * d__1;
		    t177 = 1 / t176;
		    t178 = t52 * t177;
		    t182 = t123 * -1.803772444038489 - t172 * 
			    .1210946621972718 - t119 * .01763993811759022;
		    t185 = t174 * t57 - t178 * 1. * t182;
		    t186 = t64 * t185;
		    t190 = 1 / t62 / t118;
		    t193 = t169 * 8.1290825e-4 * t186 + t61 * 
			    9.483929583333333e-4 * t190;
		    t194 = t193 * t67;
		    t195 = t194 * t59;
		    t197 = t195 * .7937005259840997 * t86;
		    t198 = t67 * t185;
		    t200 = t198 * .7937005259840997 * t86;
		    t202 = 1 / t69 / t118;
		    t204 = sigma * t202 * t85;
		    t206 = t68 * 1.0582673679788 * t204;
		    t207 = t76 * t21;
		    t210 = t81 * t21;
		    t213 = t207 * -.8333333333333333 * t119 + t210 * 
			    .8333333333333333 * t119;
		    t214 = t103 * t213;
		    t216 = t100 * .3968502629920499 * t214;
		    vrhoa[i__] = rho * t98 - t100 * .3968502629920499 * t109 
			    + t46 + t166 + t197 + t200 - t206 - t216;
		    t217 = -t93;
		    d__1 = t19 * 1.923661050931536 * t217;
		    d__2 = t42 * 1.923661050931536 * t217;
		    t222 = piecewise_(&t7, &d__1, &d__2);
		    t224 = -t108;
		    t225 = t103 * t224;
		    vrhob[i__] = rho * t222 - t100 * .3968502629920499 * t225 
			    + t46 + t166 + t197 + t200 - t206 - t216;
		    t228 = sqrt(rho);
		    t230 = 1 / t228 / t118;
		    t231 = t48 * t230;
		    t232 = t67 * t85;
		    t233 = t231 * t232;
		    t235 = t71 * t85;
		    t236 = t68 * t235;
		    vsigmaaa[i__] = t233 * -3.22602852800907e-4 + t236 * 
			    .7937005259840997;
		    vsigmaab[i__] = t233 * -6.45205705601814e-4 + t236 * 
			    1.587401051968199;
		    vsigmabb[i__] = vsigmaaa[i__];
		    t240 = t165 * 2.;
/* Computing 2nd power */
		    d__1 = t185;
		    t244 = d__1 * d__1;
		    t252 = 1 / t50 / t4;
/* Computing 2nd power */
		    d__1 = t118;
		    t253 = d__1 * d__1;
		    t254 = 1 / t253;
		    t255 = t252 * t254;
		    t257 = t118 * rho;
		    t258 = 1 / t257;
		    t259 = t122 * t258;
		    t262 = 1 / t152 * t254;
		    t264 = t171 * t258;
/* Computing 2nd power */
		    d__1 = t182;
		    t274 = d__1 * d__1;
		    t285 = (t255 * -.003207349893614542 + t259 * 
			    .009622049680843627 - t262 * 6.318986292200858e-7 
			    + t264 * 3.791391775320515e-6) * t57 - t174 * 2. *
			     t177 * t182 + t52 * 2. / t176 / t56 * t274 - 
			    t178 * 1. * (t255 * -1.202514962692326 + t259 * 
			    3.607544888076978 - t262 * .04036488739909061 + 
			    t264 * .2421893243945437 + t258 * 
			    .03527987623518044);
		    t297 = (t48 * -.0016258165 / t167 / t59 * t64 * t244 - 
			    t169 * .001896785916666667 * t190 * t185 + t169 * 
			    8.1290825e-4 * t64 * t285 - t61 * 
			    .002054851409722222 / t62 / t257) * 
			    .7937005259840997 * t67 * t59 * t86;
/* Computing 2nd power */
		    d__1 = t193;
		    t298 = d__1 * d__1;
		    t302 = t298 * .7937005259840997 * t67 * t59 * t86;
		    t306 = t195 * .7937005259840997 * t72 * t102 * t213;
		    t308 = t195 * 2.116534735957599 * t204;
		    t311 = t194 * 1.587401051968199 * t185 * t86;
		    t314 = t67 * .7937005259840997 * t285 * t86;
		    t315 = t198 * sigma;
		    t317 = t315 * .7937005259840997 * t214;
		    t319 = t198 * 2.116534735957599 * t204;
/* Computing 2nd power */
		    d__1 = t24;
		    t320 = d__1 * d__1;
		    t321 = 1 / t320;
/* Computing 2nd power */
		    d__1 = t28;
		    t324 = d__1 * d__1;
		    t325 = 1 / t324;
		    t328 = t321 * .4444444444444444 * t119 + t325 * 
			    .4444444444444444 * t119;
		    d__1 = t19 * 1.923661050931536 * t328;
		    d__2 = t42 * 1.923661050931536 * t328;
		    t333 = piecewise_(&t7, &d__1, &d__2);
		    t334 = rho * t333;
/* Computing 2nd power */
		    d__1 = t83;
		    t335 = d__1 * d__1;
		    t338 = t71 / t84 / t335;
/* Computing 2nd power */
		    d__1 = t108;
		    t339 = d__1 * d__1;
		    t343 = 1 / t75;
		    t346 = 1 / t80;
		    t349 = t343 * .2777777777777778 * t119 + t346 * 
			    .2777777777777778 * t119;
		    t352 = t100 * .3968502629920499 * t103 * t349;
		    t353 = t240 + t297 + t302 - t306 - t308 + t311 + t314 - 
			    t317 - t319 + t334 + t100 * .5952753944880748 * 
			    t338 * t339 - t352;
		    t355 = t202 * t102;
		    t358 = t100 * 1.0582673679788 * t355 * t213;
		    t364 = t68 * 2.469290525283866 * sigma / t69 / t257 * t85;
		    t375 = t343 * -.2777777777777778 * t21 * t258 - t76 * 
			    .8333333333333333 * t119 - t346 * 
			    .2777777777777778 * t21 * t258 + t81 * 
			    .8333333333333333 * t119;
		    t377 = t100 * t103 * t375;
		    t379 = t315 * t109;
		    t383 = t195 * t72 * t102 * t108;
		    t397 = t321 * -.4444444444444444 * t21 * t258 - t24 * 
			    1.333333333333333 * t119 - t325 * 
			    .4444444444444444 * t21 * t258 + t28 * 
			    1.333333333333333 * t119;
		    d__1 = t135 * 1.923661050931536 * t93 + t19 * 
			    1.923661050931536 * t397;
		    d__2 = t159 * 1.923661050931536 * t93 + t42 * 
			    1.923661050931536 * t397;
		    t406 = piecewise_(&t7, &d__1, &d__2);
		    t407 = rho * t406;
		    t411 = t100 * t338 * t213 * t108;
		    t414 = t100 * t355 * t108;
/* Computing 2nd power */
		    d__1 = t213;
		    t416 = d__1 * d__1;
		    t419 = t100 * .5952753944880748 * t338 * t416;
/* Computing 2nd power */
		    d__1 = t21;
		    t420 = d__1 * d__1;
		    t434 = t100 * .3968502629920499 * t103 * (t343 * 
			    .2777777777777778 * t420 * t254 + t207 * 
			    1.666666666666667 * t258 + t346 * 
			    .2777777777777778 * t420 * t254 - t210 * 
			    1.666666666666667 * t258);
/* Computing 2nd power */
		    d__1 = t125;
		    t437 = d__1 * d__1;
		    t439 = .2846 / t112 / t11 * t437;
		    t442 = 1 / t116 / t4 * t254;
		    t444 = t117 * t258;
		    t450 = t113 * .1423 * (t442 * -.1151789710300888 + t444 * 
			    .2764295304722132 - t255 * .04596107859241333 + 
			    t259 * .13788323577724);
/* Computing 2nd power */
		    d__1 = t132;
		    t453 = d__1 * d__1;
		    t478 = t321 * .4444444444444444 * t420 * t254 + t138 * 
			    2.666666666666667 * t258 + t325 * 
			    .4444444444444444 * t420 * t254 - t141 * 
			    2.666666666666667 * t258;
		    t484 = t252 * t34 * t254;
		    t487 = t149 * t258;
		    t489 = t5 * t119;
		    d__1 = -t439 + t450 + (-.1686 / t128 / t16 * t453 + t129 *
			     .0843 * (t442 * -.1529411334382821 + t444 * 
			    .367058720251877 - t255 * .03599411403862963 + 
			    t259 * .1079823421158889) + t439 - t450) * 
			    1.923661050931536 * t30 + t135 * 
			    3.847322101863073 * t144 + t19 * 
			    1.923661050931536 * t478;
		    d__2 = t119 * .01036666666666667 - t484 * 
			    2.757113292886222e-4 - t259 * .004521665800333405 
			    + t487 * 8.271339878658667e-4 + t489 * 
			    4.135669939329333e-4 + t255 * .001599125709874009 
			    + (t119 * -.005183333333333333 + t484 * 
			    1.792123640376044e-4 + t259 * .002633043194706342 
			    - t487 * 5.376370921128133e-4 - t489 * 
			    2.688185460564067e-4 - t255 * 
			    9.374185195813156e-4) * 1.923661050931536 * t30 + 
			    t159 * 3.847322101863073 * t144 + t42 * 
			    1.923661050931536 * t478;
		    t506 = piecewise_(&t7, &d__1, &d__2);
		    t507 = rho * t506;
		    t508 = t98 * 2. + t358 + t364 - t377 * .7937005259840997 
			    - t379 * .7937005259840997 - t383 * 
			    .7937005259840997 + t407 * 2. + t411 * 
			    1.19055078897615 + t414 * 1.0582673679788 + t419 
			    - t434 + t507;
		    v2rhoa2[i__] = t353 + t508;
/* Computing 2nd power */
		    d__1 = t224;
		    t510 = d__1 * d__1;
		    t516 = -t397;
		    d__1 = t135 * 1.923661050931536 * t217 + t19 * 
			    1.923661050931536 * t516;
		    d__2 = t159 * 1.923661050931536 * t217 + t42 * 
			    1.923661050931536 * t516;
		    t525 = piecewise_(&t7, &d__1, &d__2);
		    t526 = rho * t525;
		    t530 = -t100 * t103 * t375;
		    t532 = t315 * t225;
		    t534 = t240 + t297 + t222 * 2. + t100 * .5952753944880748 
			    * t338 * t510 + t302 - t306 - t308 + t311 + t314 
			    + t526 * 2. - t530 * .7937005259840997 - t532 * 
			    .7937005259840997;
		    t536 = t100 * t355 * t224;
		    t540 = t100 * t338 * t213 * t224;
		    t544 = t195 * t72 * t102 * t224;
		    t546 = t536 * 1.0582673679788 + t540 * 1.19055078897615 - 
			    t544 * .7937005259840997 - t317 - t319 + t334 - 
			    t352 + t358 + t364 + t419 - t434 + t507;
		    v2rhob2[i__] = t534 + t546;
		    t560 = -t328;
		    d__1 = t19 * 1.923661050931536 * t560;
		    d__2 = t42 * 1.923661050931536 * t560;
		    t565 = piecewise_(&t7, &d__1, &d__2);
		    t567 = t526 - t530 * .3968502629920499 - t532 * 
			    .3968502629920499 + t536 * .5291336839893998 + 
			    t540 * .5952753944880748 - t544 * 
			    .3968502629920499 + t100 * .3968502629920499 * 
			    t103 * t349 + t100 * .5952753944880748 * t338 * 
			    t108 * t224 + t98 + rho * t565 + t240 + t222 + 
			    t297 + t302 - t306;
		    t573 = -t308 + t311 + t314 - t319 - t317 + t358 + t364 - 
			    t377 * .3968502629920499 - t379 * 
			    .3968502629920499 - t383 * .3968502629920499 + 
			    t407 + t411 * .5952753944880748 + t507 + t414 * 
			    .5291336839893998 - t434 + t419;
		    v2rhoab[i__] = t567 + t573;
		    t574 = t67 * t102;
		    t576 = t231 * t574 * t108;
		    t578 = t68 * t109;
		    t580 = 1 / t48;
		    t590 = (t580 * 4.06454125e-4 * t168 * t186 + t580 * 
			    4.741964791666667e-4 * t60 * t190) * t67 * t59 * 
			    t86;
		    t591 = t590 * .7937005259840997;
		    t595 = t193 * t48 * t230 * t67 * t85;
		    t596 = t595 * 3.22602852800907e-4;
		    t599 = t194 * t59 * t71 * t85;
		    t600 = t599 * .7937005259840997;
		    t603 = t61 * t230 * t198 * t85;
		    t604 = t603 * 3.22602852800907e-4;
		    t605 = t198 * t235;
		    t606 = t605 * .7937005259840997;
		    t610 = t48 / t228 / t257 * t232;
		    t611 = t610 * 4.30137137067876e-4;
		    t613 = t68 * t202 * t85;
		    t614 = t613 * 1.0582673679788;
		    t616 = t231 * t574 * t213;
		    t617 = t616 * 1.613014264004535e-4;
		    t618 = t68 * t214;
		    t619 = t618 * .3968502629920499;
		    v2rhoasigmaaa[i__] = t576 * 1.613014264004535e-4 - t578 * 
			    .3968502629920499 + t591 - t596 + t600 - t604 + 
			    t606 + t611 - t614 + t617 - t619;
		    t622 = t590 * 1.587401051968199;
		    t623 = t595 * 6.45205705601814e-4;
		    t624 = t599 * 1.587401051968199;
		    t625 = t603 * 6.45205705601814e-4;
		    t626 = t605 * 1.587401051968199;
		    t627 = t610 * 8.602742741357521e-4;
		    t628 = t613 * 2.116534735957599;
		    t629 = t616 * 3.22602852800907e-4;
		    t630 = t618 * .7937005259840997;
		    v2rhoasigmaab[i__] = t576 * 3.22602852800907e-4 - t578 * 
			    .7937005259840997 + t622 - t623 + t624 - t625 + 
			    t626 + t627 - t628 + t629 - t630;
		    v2rhoasigmabb[i__] = v2rhoasigmaaa[i__];
		    t632 = t231 * t574 * t224;
		    t634 = t68 * t225;
		    v2rhobsigmaaa[i__] = t632 * 1.613014264004535e-4 - t634 * 
			    .3968502629920499 + t591 - t596 + t600 - t604 + 
			    t606 + t611 - t614 + t617 - t619;
		    v2rhobsigmaab[i__] = t632 * 3.22602852800907e-4 - t634 * 
			    .7937005259840997 + t622 - t623 + t624 - t625 + 
			    t626 + t627 - t628 + t629 - t630;
		    v2rhobsigmabb[i__] = v2rhobsigmaaa[i__];
		    t639 = t580 * t230 * t232;
/* Computing 2nd power */
		    d__1 = t69;
		    t641 = d__1 * d__1;
		    t645 = 1 / t641 / t257 * t60 * t232;
		    v2sigmaaa2[i__] = t639 * -4.839042792013605e-4 + t645 * 
			    1.311232602576965e-7;
		    v2sigmaaaab[i__] = t639 * -9.678085584027211e-4 + t645 * 
			    2.622465205153929e-7;
		    v2sigmaaabb[i__] = v2sigmaaa2[i__];
		    v2sigmaab2[i__] = t639 * -.001935617116805442 + t645 * 
			    5.244930410307859e-7;
		    v2sigmaabbb[i__] = v2sigmaaaab[i__];
		    v2sigmabb2[i__] = v2sigmaaabb[i__];
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
} /* uks_c_p86__ */

/* Subroutine */ int rks_c_p86__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    logical L__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), sqrt(
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4;
    static logical t5;
    static doublereal t6, t9, t20, t12, t22, t31, t41, t34, t24, t18, t28, 
	    t29, t39, t33, t36, t40, t43, t44, t47, t49, t51, t52, t55, t56, 
	    t59, t60, t73, t75, t80, t84, t93, t11, t14, t32, t48, t50, t53, 
	    t54, t62, t66, t69, t76, t77, t79, t82, t85, t86, t90, t94, t98, 
	    t101, t102, t103, t104, t115, t107, t111, t112, t117, t129, t135, 
	    t137, t139, t147, t152, t153, t156, t160, t161, t163, t180, t185, 
	    t193, t203, t205, t215, t226, t247, t254, t266, t287, rho, sigma;
    extern doublereal piecewise_(logical *, doublereal *, doublereal *);


/*     J.P. Perdew */
/*     Density-functional approximation for the correlation energy of */
/*     the inhomogeneous electron gas */
/*     Phys. Rev. B33 (1986) 8822-8824 */


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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b3);
		t4 = t3 * .6203504908994;
		t6 = pow_dd(&t2, &c_b4);
		t12 = log(t4);
		L__1 = 1. <= t4;
		d__1 = -.1423 / (t6 * .8292885914166397 + 1. + t3 * 
			.20682485366586);
		d__2 = t12 * .0311 - .048 + t3 * .0012407009817988 * t12 - t3 
			* .00719606569443304;
		t18 = piecewise_(&L__1, &d__1, &d__2);
		t20 = sqrt(sigma);
/* Computing 2nd power */
		d__1 = t3;
		t22 = d__1 * d__1;
		t31 = (t3 * .01443307452126544 + .002568 + t22 * 
			2.843543831490386e-6) / (t3 * 5.411317332115466 + 1. 
			+ t22 * .1816419932959077 + t2 * .01763993811759022) 
			+ .001667;
		t34 = pow_dd(&rho, &c_b4);
		t39 = exp(t20 * -8.1290825e-4 / t31 / t34 / rho);
		t41 = pow_dd(&rho, &c_b3);
		zk[i__] = rho * t18 + t39 * 1. * t31 * sigma / t41 / rho;
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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b3);
		t4 = t3 * .6203504908994;
		t5 = 1. <= t4;
		t6 = pow_dd(&t2, &c_b4);
		t9 = t6 * .8292885914166397 + 1. + t3 * .20682485366586;
		t12 = log(t4);
		d__1 = -.1423 / t9;
		d__2 = t12 * .0311 - .048 + t3 * .0012407009817988 * t12 - t3 
			* .00719606569443304;
		t18 = piecewise_(&t5, &d__1, &d__2);
		t20 = sqrt(sigma);
/* Computing 2nd power */
		d__1 = t3;
		t22 = d__1 * d__1;
		t24 = t3 * .01443307452126544 + .002568 + t22 * 
			2.843543831490386e-6;
		t28 = t3 * 5.411317332115466 + 1. + t22 * .1816419932959077 + 
			t2 * .01763993811759022;
		t29 = 1 / t28;
		t31 = t24 * t29 + .001667;
		t33 = t20 / t31;
		t34 = pow_dd(&rho, &c_b4);
		t36 = 1 / t34 / rho;
		t39 = exp(t33 * -8.1290825e-4 * t36);
		t40 = t39 * t31;
		t41 = pow_dd(&rho, &c_b3);
		t43 = 1 / t41 / rho;
		t44 = sigma * t43;
		zk[i__] = rho * t18 + t40 * 1. * t44;
		t47 = piecewise_(&t5, &c_b60, &c_b60);
/* Computing 2nd power */
		d__1 = t9;
		t49 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t6;
		t51 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t51;
		t52 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = rho;
		t55 = d__1 * d__1;
		t56 = 1 / t55;
		t59 = 1 / t22;
		t60 = t59 * t56;
		d__1 = .1423 / t49 * (-.1382147652361066 / t52 / t6 * t56 - 
			t60 * .06894161788861999);
		d__2 = t2 * -.01036666666666667 - t59 * 4.135669939329333e-4 *
			 t12 * t56 - t3 * 4.135669939329333e-4 * t2 + t60 * 
			.002398688564811013;
		t73 = piecewise_(&t5, &d__1, &d__2);
/* Computing 2nd power */
		d__1 = t31;
		t75 = d__1 * d__1;
		t80 = 1 / t3 * t56;
/* Computing 2nd power */
		d__1 = t28;
		t84 = d__1 * d__1;
		t93 = (t60 * -.004811024840421814 - t80 * 
			1.895695887660258e-6) * t29 - t24 * 1. / t84 * (t60 * 
			-1.803772444038489 - t80 * .1210946621972718 - t56 * 
			.01763993811759022);
		vrhoa[i__] = rho * t47 + t18 + rho * t73 + (t20 * 
			8.1290825e-4 / t75 * t36 * t93 + t33 * 
			9.483929583333333e-4 / t34 / t55) * 1. * t39 * t31 * 
			sigma * t43 + t39 * 1. * t93 * t44 - t40 * 
			1.333333333333333 * sigma / t41 / t55;
		t115 = sqrt(rho);
		vsigmaaa[i__] = t20 * -.0016258165 / t115 / t55 * t39 + t40 * 
			4. * t43;
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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b3);
		t4 = t3 * .6203504908994;
		t5 = 1. <= t4;
		t6 = pow_dd(&t2, &c_b4);
		t9 = t6 * .8292885914166397 + 1. + t3 * .20682485366586;
		t11 = .1423 / t9;
		t12 = log(t4);
		t14 = t3 * t12;
		d__1 = -t11;
		d__2 = t12 * .0311 - .048 + t14 * .0012407009817988 - t3 * 
			.00719606569443304;
		t18 = piecewise_(&t5, &d__1, &d__2);
		t20 = sqrt(sigma);
/* Computing 2nd power */
		d__1 = t3;
		t22 = d__1 * d__1;
		t24 = t3 * .01443307452126544 + .002568 + t22 * 
			2.843543831490386e-6;
		t28 = t3 * 5.411317332115466 + 1. + t22 * .1816419932959077 + 
			t2 * .01763993811759022;
		t29 = 1 / t28;
		t31 = t24 * t29 + .001667;
		t32 = 1 / t31;
		t33 = t20 * t32;
		t34 = pow_dd(&rho, &c_b4);
		t36 = 1 / t34 / rho;
		t39 = exp(t33 * -8.1290825e-4 * t36);
		t40 = t39 * t31;
		t41 = pow_dd(&rho, &c_b3);
		t43 = 1 / t41 / rho;
		t44 = sigma * t43;
		zk[i__] = rho * t18 + t40 * 1. * t44;
		t47 = piecewise_(&t5, &c_b60, &c_b60);
		t48 = rho * t47;
/* Computing 2nd power */
		d__1 = t9;
		t49 = d__1 * d__1;
		t50 = 1 / t49;
/* Computing 2nd power */
		d__1 = t6;
		t51 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t51;
		t52 = d__1 * d__1;
		t53 = t52 * t6;
		t54 = 1 / t53;
/* Computing 2nd power */
		d__1 = rho;
		t55 = d__1 * d__1;
		t56 = 1 / t55;
		t59 = 1 / t22;
		t60 = t59 * t56;
		t62 = t54 * -.1382147652361066 * t56 - t60 * 
			.06894161788861999;
		t66 = t59 * t12;
		t69 = t3 * t2;
		d__1 = t50 * .1423 * t62;
		d__2 = t2 * -.01036666666666667 - t66 * 4.135669939329333e-4 *
			 t56 - t69 * 4.135669939329333e-4 + t60 * 
			.002398688564811013;
		t73 = piecewise_(&t5, &d__1, &d__2);
/* Computing 2nd power */
		d__1 = t31;
		t75 = d__1 * d__1;
		t76 = 1 / t75;
		t77 = t20 * t76;
		t79 = 1 / t3;
		t80 = t79 * t56;
		t82 = t60 * -.004811024840421814 - t80 * 1.895695887660258e-6;
/* Computing 2nd power */
		d__1 = t28;
		t84 = d__1 * d__1;
		t85 = 1 / t84;
		t86 = t24 * t85;
		t90 = t60 * -1.803772444038489 - t80 * .1210946621972718 - 
			t56 * .01763993811759022;
		t93 = t82 * t29 - t86 * 1. * t90;
		t94 = t36 * t93;
		t98 = 1 / t34 / t55;
		t101 = t77 * 8.1290825e-4 * t94 + t33 * 9.483929583333333e-4 *
			 t98;
		t102 = t101 * t39;
		t103 = t31 * sigma;
		t104 = t103 * t43;
		t107 = t39 * t93;
		t111 = 1 / t41 / t55;
		t112 = sigma * t111;
		vrhoa[i__] = t48 + t18 + rho * t73 + t102 * 1. * t104 + t107 *
			 1. * t44 - t40 * 1.333333333333333 * t112;
		t115 = sqrt(rho);
		t117 = 1 / t115 / t55;
		vsigmaaa[i__] = t20 * -.0016258165 * t117 * t39 + t40 * 4. * 
			t43;
		t129 = (-.0843 / (t6 * 1.101176160755631 + 1. + t3 * 
			.1619735131738333) + t11) * t56;
		t135 = (t12 * -.01555 + .0211 - t14 * 8.0645563816922e-4 + t3 
			* .00421838333811592) * t56;
		d__1 = t129 * 1.709920934161366;
		d__2 = t135 * 1.709920934161366;
		t137 = piecewise_(&t5, &d__1, &d__2);
		t139 = t55 * rho;
/* Computing 2nd power */
		d__1 = t62;
		t147 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t55;
		t152 = d__1 * d__1;
		t153 = 1 / t152;
		t156 = 1 / t139;
		t160 = 1 / t22 / t2;
		t161 = t160 * t153;
		t163 = t59 * t156;
		d__1 = -.2846 / t49 / t9 * t147 + t50 * .1423 * (
			-.1151789710300888 / t53 / t2 * t153 + t54 * 
			.2764295304722132 * t156 - t161 * .04596107859241333 
			+ t163 * .13788323577724);
		d__2 = t56 * .01036666666666667 - t160 * 2.757113292886222e-4 
			* t12 * t153 - t163 * .004521665800333405 + t66 * 
			8.271339878658667e-4 * t156 + t3 * 
			4.135669939329333e-4 * t56 + t161 * 
			.001599125709874009;
		t180 = piecewise_(&t5, &d__1, &d__2);
		d__1 = t129 * -1.709920934161366;
		d__2 = t135 * -1.709920934161366;
		t185 = piecewise_(&t5, &d__1, &d__2);
/* Computing 2nd power */
		d__1 = t93;
		t193 = d__1 * d__1;
		t203 = 1 / t69 * t153;
		t205 = t79 * t156;
/* Computing 2nd power */
		d__1 = t90;
		t215 = d__1 * d__1;
		t226 = (t161 * -.003207349893614542 + t163 * 
			.009622049680843627 - t203 * 6.318986292200858e-7 + 
			t205 * 3.791391775320515e-6) * t29 - t82 * 2. * t85 * 
			t90 + t24 * 2. / t84 / t28 * t215 - t86 * 1. * (t161 *
			 -1.202514962692326 + t163 * 3.607544888076978 - t203 
			* .04036488739909061 + t205 * .2421893243945437 + 
			t156 * .03527987623518044);
/* Computing 2nd power */
		d__1 = t101;
		t247 = d__1 * d__1;
		v2rhoa2[i__] = rho * t137 + t40 * 6.222222222222222 * sigma / 
			t41 / t139 + rho * 2. * t180 + rho * t185 - t102 * 
			5.333333333333333 * t103 * t111 + (t20 * -.0016258165 
			/ t75 / t31 * t36 * t193 - t77 * .001896785916666667 *
			 t98 * t93 + t77 * 8.1290825e-4 * t36 * t226 - t33 * 
			.002054851409722222 / t34 / t139) * 2. * t39 * t104 + 
			t102 * 4. * t93 * sigma * t43 + t39 * 2. * t226 * t44 
			- t107 * 5.333333333333333 * t112 + t247 * 2. * t39 * 
			t104 + t48 * 4. + t47 * 4. + t73 * 4.;
		t254 = 1 / t20;
		t266 = t117 * t39;
		v2rhoasigmaaa[i__] = (t254 * 4.06454125e-4 * t76 * t94 + t254 
			* 4.741964791666667e-4 * t32 * t98) * 4. * t39 * t104 
			- t101 * .0016258165 * t20 * t266 + t102 * 4. * t31 * 
			t43 - t33 * .0016258165 * t266 * t93 + t107 * 4. * 
			t43 + t20 * .002167755333333333 / t115 / t139 * t39 - 
			t40 * 5.333333333333333 * t111;
/* Computing 2nd power */
		d__1 = t41;
		t287 = d__1 * d__1;
		v2sigmaaa2[i__] = t254 * -.009754899 * t117 * t39 + 
			2.64327929167225e-6 / t287 / t139 * t32 * t39;
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
} /* rks_c_p86__ */


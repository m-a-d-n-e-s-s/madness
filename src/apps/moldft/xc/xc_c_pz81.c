/* xc_c_pz81.f -- translated by f2c (version 20050501).
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
static doublereal c_b3 = .16666666666666666;
static doublereal c_b31 = 0.;

/* :C_PZ81subrstart */
/*    Generated: Tue Mar  9 13:25:27 GMT 2004 */
/* Subroutine */ int uks_c_pz81__(integer *ideriv, integer *npt, doublereal *
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
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t1, t2, t3;
    static logical t4;
    static doublereal t5, t8, t10, t11, t20, t21, t31, t24, t25, t17, t27, 
	    t19, t33, t43, t18, t28, t29, t30, t13, t16, t39, t48, t53, t55, 
	    t57, t58, t61, t62, t63, t65, t66, t67, t71, t72, t88, t94, t96, 
	    t22, t23, t32, t36, t54, t82, t56, t59, t60, t69, t73, t76, t79, 
	    t110, t111, t85, t93, t103, t122, t116, t123, t109, t118, t119, 
	    t126, t131, t132, t138, t148, t157, t158, t160, t163, t165, t168, 
	    t169, t170, t172, t175, t176, t178, t182, t185, t200, t211, t217, 
	    t220, t222, t239, t240, t244, t253, t254, t256, t261, rho, rhoa, 
	    rhob;
    extern doublereal piecewise_(logical *, doublereal *, doublereal *);


/*     J.P. Perdew, and A. Zunger */
/*     Self-interaction correction to density-functional approximations */
/*     for many-electron systems */
/*     Phys. Rev. B23 (1981) 5048-5079 */


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
		    t1 = 1 / rhob;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t5 = pow_dd(&t1, &c_b3);
		    t11 = log(t3);
		    L__1 = 1. <= t3;
		    d__1 = -.0843 / (t5 * 1.101176160755631 + 1. + t2 * 
			    .1619735131738333);
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&L__1, &d__1, &d__2);
		    zk[i__] = rhob * t17;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t5 = pow_dd(&t1, &c_b3);
		    t11 = log(t3);
		    L__1 = 1. <= t3;
		    d__1 = -.0843 / (t5 * 1.101176160755631 + 1. + t2 * 
			    .1619735131738333);
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&L__1, &d__1, &d__2);
		    zk[i__] = rhoa * t17;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = 1 / rho;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t5 = pow_dd(&t1, &c_b3);
		    t10 = .1423 / (t5 * .8292885914166397 + 1. + t2 * 
			    .20682485366586);
		    t19 = (rhoa - rhob * 1.) * t1;
		    t20 = t19 + 1.;
		    t21 = pow_dd(&t20, &c_b2);
		    t24 = 1. - t19 * 1.;
		    t25 = pow_dd(&t24, &c_b2);
		    t27 = t21 * t20 + t25 * t24 - 2.;
		    t31 = log(t3);
		    t33 = t2 * t31;
		    L__1 = 1. <= t3;
		    d__1 = -t10 + (-.0843 / (t5 * 1.101176160755631 + 1. + t2 
			    * .1619735131738333) + t10) * 1.923661050931536 * 
			    t27;
		    d__2 = t31 * .0311 - .048 + t33 * .0012407009817988 - t2 *
			     .00719606569443304 + (t31 * -.01555 + .0211 - 
			    t33 * 8.0645563816922e-4 + t2 * 
			    .00421838333811592) * 1.923661050931536 * t27;
		    t43 = piecewise_(&L__1, &d__1, &d__2);
		    zk[i__] = rho * t43;
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
		    t1 = 1 / rhob;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * 1.101176160755631 + 1. + t2 * .1619735131738333;
		    t11 = log(t3);
		    d__1 = -.0843 / t8;
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rhob * t17;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t8;
		    t18 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t5;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t21 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t24 = d__1 * d__1;
		    t25 = 1 / t24;
/* Computing 2nd power */
		    d__1 = t2;
		    t28 = d__1 * d__1;
		    t29 = 1 / t28;
		    t30 = t29 * t25;
		    d__1 = .0843 / t18 * (-.1835293601259385 / t21 / t5 * t25 
			    - t30 * .05399117105794445);
		    d__2 = t1 * -.005183333333333333 - t29 * 
			    1.447484478765267e-4 * t11 * t25 - t2 * 
			    1.447484478765267e-4 * t1 + t30 * 
			    9.9256078543904e-4;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    vrhob[i__] = t17 + rhob * t43;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * 1.101176160755631 + 1. + t2 * .1619735131738333;
		    t11 = log(t3);
		    d__1 = -.0843 / t8;
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rhoa * t17;
/* Computing 2nd power */
		    d__1 = t8;
		    t18 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t5;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t21 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t24 = d__1 * d__1;
		    t25 = 1 / t24;
/* Computing 2nd power */
		    d__1 = t2;
		    t28 = d__1 * d__1;
		    t29 = 1 / t28;
		    t30 = t29 * t25;
		    d__1 = .0843 / t18 * (-.1835293601259385 / t21 / t5 * t25 
			    - t30 * .05399117105794445);
		    d__2 = t1 * -.005183333333333333 - t29 * 
			    1.447484478765267e-4 * t11 * t25 - t2 * 
			    1.447484478765267e-4 * t1 + t30 * 
			    9.9256078543904e-4;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    vrhoa[i__] = t17 + rhoa * t43;
		    vrhob[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = 1 / rho;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * .8292885914166397 + 1. + t2 * .20682485366586;
		    t10 = .1423 / t8;
		    t13 = t5 * 1.101176160755631 + 1. + t2 * 
			    .1619735131738333;
		    t16 = -.0843 / t13 + t10;
		    t18 = rhoa - rhob * 1.;
		    t19 = t18 * t1;
		    t20 = t19 + 1.;
		    t21 = pow_dd(&t20, &c_b2);
		    t24 = 1. - t19 * 1.;
		    t25 = pow_dd(&t24, &c_b2);
		    t27 = t21 * t20 + t25 * t24 - 2.;
		    t31 = log(t3);
		    t33 = t2 * t31;
		    t39 = t31 * -.01555 + .0211 - t33 * 8.0645563816922e-4 + 
			    t2 * .00421838333811592;
		    d__1 = -t10 + t16 * 1.923661050931536 * t27;
		    d__2 = t31 * .0311 - .048 + t33 * .0012407009817988 - t2 *
			     .00719606569443304 + t39 * 1.923661050931536 * 
			    t27;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rho * t43;
		    t48 = t21 * 1.333333333333333 * t1 - t25 * 
			    1.333333333333333 * t1;
		    d__1 = t16 * 1.923661050931536 * t48;
		    d__2 = t39 * 1.923661050931536 * t48;
		    t53 = piecewise_(&t4, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t8;
		    t55 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t5;
		    t57 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t57;
		    t58 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rho;
		    t61 = d__1 * d__1;
		    t62 = 1 / t61;
		    t63 = 1 / t58 / t5 * t62;
/* Computing 2nd power */
		    d__1 = t2;
		    t65 = d__1 * d__1;
		    t66 = 1 / t65;
		    t67 = t66 * t62;
		    t71 = .1423 / t55 * (t63 * -.1382147652361066 - t67 * 
			    .06894161788861999);
/* Computing 2nd power */
		    d__1 = t13;
		    t72 = d__1 * d__1;
		    t88 = t21 * -1.333333333333333 * t18 * t62 + t25 * 
			    1.333333333333333 * t18 * t62;
		    t94 = t66 * t31 * t62;
		    t96 = t2 * t1;
		    d__1 = t71 + (.0843 / t72 * (t63 * -.1835293601259385 - 
			    t67 * .05399117105794445) - t71) * 
			    1.923661050931536 * t27 + t16 * 1.923661050931536 
			    * t88;
		    d__2 = t1 * -.01036666666666667 - t94 * 
			    4.135669939329333e-4 - t96 * 4.135669939329333e-4 
			    + t67 * .002398688564811013 + (t1 * 
			    .005183333333333333 + t94 * 2.688185460564067e-4 
			    + t96 * 2.688185460564067e-4 - t67 * 
			    .001406127779371973) * 1.923661050931536 * t27 + 
			    t39 * 1.923661050931536 * t88;
		    t109 = piecewise_(&t4, &d__1, &d__2);
		    t110 = rho * t109;
		    vrhoa[i__] = rho * t53 + t43 + t110;
		    t111 = -t48;
		    d__1 = t16 * 1.923661050931536 * t111;
		    d__2 = t39 * 1.923661050931536 * t111;
		    t116 = piecewise_(&t4, &d__1, &d__2);
		    vrhob[i__] = rho * t116 + t43 + t110;
		}
/* rhoa,rhob */
	    } else {
/* rho */
		zk[i__] = 0.;
		vrhoa[i__] = 0.;
		vrhob[i__] = 0.;
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
		    t1 = 1 / rhob;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * 1.101176160755631 + 1. + t2 * .1619735131738333;
		    t11 = log(t3);
		    d__1 = -.0843 / t8;
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rhob * t17;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t8;
		    t18 = d__1 * d__1;
		    t19 = 1 / t18;
/* Computing 2nd power */
		    d__1 = t5;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t21 = d__1 * d__1;
		    t22 = t21 * t5;
		    t23 = 1 / t22;
/* Computing 2nd power */
		    d__1 = rhob;
		    t24 = d__1 * d__1;
		    t25 = 1 / t24;
/* Computing 2nd power */
		    d__1 = t2;
		    t28 = d__1 * d__1;
		    t29 = 1 / t28;
		    t30 = t29 * t25;
		    t32 = t23 * -.1835293601259385 * t25 - t30 * 
			    .05399117105794445;
		    t36 = t29 * t11;
		    d__1 = t19 * .0843 * t32;
		    d__2 = t1 * -.005183333333333333 - t36 * 
			    1.447484478765267e-4 * t25 - t2 * 
			    1.447484478765267e-4 * t1 + t30 * 
			    9.9256078543904e-4;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    vrhob[i__] = t17 + rhob * t43;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t32;
		    t48 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t24;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t58 = 1 / t24 / rhob;
		    t62 = 1 / t28 / t1;
		    t63 = t62 * t54;
		    t65 = t29 * t58;
		    d__1 = -.1686 / t18 / t8 * t48 + t19 * .0843 * (
			    -.1529411334382821 / t22 / t1 * t54 + t23 * 
			    .367058720251877 * t58 - t63 * .03599411403862963 
			    + t65 * .1079823421158889);
		    d__2 = t25 * .005183333333333333 - t62 * 
			    9.649896525101778e-5 * t11 * t54 - t65 * 
			    .001888622605627062 + t36 * 2.894968957530533e-4 *
			     t58 + t2 * 1.447484478765267e-4 * t25 + t63 * 
			    6.617071902926934e-4;
		    t82 = piecewise_(&t4, &d__1, &d__2);
		    v2rhob2[i__] = t43 * 2. + rhob * t82;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * 1.101176160755631 + 1. + t2 * .1619735131738333;
		    t11 = log(t3);
		    d__1 = -.0843 / t8;
		    d__2 = t11 * .01555 - .0269 + t2 * 4.3424534362958e-4 * 
			    t11 - t2 * .00297768235631712;
		    t17 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rhoa * t17;
/* Computing 2nd power */
		    d__1 = t8;
		    t18 = d__1 * d__1;
		    t19 = 1 / t18;
/* Computing 2nd power */
		    d__1 = t5;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t21 = d__1 * d__1;
		    t22 = t21 * t5;
		    t23 = 1 / t22;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t24 = d__1 * d__1;
		    t25 = 1 / t24;
/* Computing 2nd power */
		    d__1 = t2;
		    t28 = d__1 * d__1;
		    t29 = 1 / t28;
		    t30 = t29 * t25;
		    t32 = t23 * -.1835293601259385 * t25 - t30 * 
			    .05399117105794445;
		    t36 = t29 * t11;
		    d__1 = t19 * .0843 * t32;
		    d__2 = t1 * -.005183333333333333 - t36 * 
			    1.447484478765267e-4 * t25 - t2 * 
			    1.447484478765267e-4 * t1 + t30 * 
			    9.9256078543904e-4;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    vrhoa[i__] = t17 + rhoa * t43;
		    vrhob[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t32;
		    t48 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t24;
		    t53 = d__1 * d__1;
		    t54 = 1 / t53;
		    t58 = 1 / t24 / rhoa;
		    t62 = 1 / t28 / t1;
		    t63 = t62 * t54;
		    t65 = t29 * t58;
		    d__1 = -.1686 / t18 / t8 * t48 + t19 * .0843 * (
			    -.1529411334382821 / t22 / t1 * t54 + t23 * 
			    .367058720251877 * t58 - t63 * .03599411403862963 
			    + t65 * .1079823421158889);
		    d__2 = t25 * .005183333333333333 - t62 * 
			    9.649896525101778e-5 * t11 * t54 - t65 * 
			    .001888622605627062 + t36 * 2.894968957530533e-4 *
			     t58 + t2 * 1.447484478765267e-4 * t25 + t63 * 
			    6.617071902926934e-4;
		    t82 = piecewise_(&t4, &d__1, &d__2);
		    v2rhoa2[i__] = t43 * 2. + rhoa * t82;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = 1 / rho;
		    t2 = pow_dd(&t1, &c_b2);
		    t3 = t2 * .6203504908994;
		    t4 = 1. <= t3;
		    t5 = pow_dd(&t1, &c_b3);
		    t8 = t5 * .8292885914166397 + 1. + t2 * .20682485366586;
		    t10 = .1423 / t8;
		    t13 = t5 * 1.101176160755631 + 1. + t2 * 
			    .1619735131738333;
		    t16 = -.0843 / t13 + t10;
		    t18 = rhoa - rhob * 1.;
		    t19 = t18 * t1;
		    t20 = t19 + 1.;
		    t21 = pow_dd(&t20, &c_b2);
		    t24 = 1. - t19 * 1.;
		    t25 = pow_dd(&t24, &c_b2);
		    t27 = t21 * t20 + t25 * t24 - 2.;
		    t31 = log(t3);
		    t33 = t2 * t31;
		    t39 = t31 * -.01555 + .0211 - t33 * 8.0645563816922e-4 + 
			    t2 * .00421838333811592;
		    d__1 = -t10 + t16 * 1.923661050931536 * t27;
		    d__2 = t31 * .0311 - .048 + t33 * .0012407009817988 - t2 *
			     .00719606569443304 + t39 * 1.923661050931536 * 
			    t27;
		    t43 = piecewise_(&t4, &d__1, &d__2);
		    zk[i__] = rho * t43;
		    t48 = t21 * 1.333333333333333 * t1 - t25 * 
			    1.333333333333333 * t1;
		    d__1 = t16 * 1.923661050931536 * t48;
		    d__2 = t39 * 1.923661050931536 * t48;
		    t53 = piecewise_(&t4, &d__1, &d__2);
/* Computing 2nd power */
		    d__1 = t8;
		    t55 = d__1 * d__1;
		    t56 = 1 / t55;
/* Computing 2nd power */
		    d__1 = t5;
		    t57 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t57;
		    t58 = d__1 * d__1;
		    t59 = t58 * t5;
		    t60 = 1 / t59;
/* Computing 2nd power */
		    d__1 = rho;
		    t61 = d__1 * d__1;
		    t62 = 1 / t61;
		    t63 = t60 * t62;
/* Computing 2nd power */
		    d__1 = t2;
		    t65 = d__1 * d__1;
		    t66 = 1 / t65;
		    t67 = t66 * t62;
		    t69 = t63 * -.1382147652361066 - t67 * .06894161788861999;
		    t71 = t56 * .1423 * t69;
/* Computing 2nd power */
		    d__1 = t13;
		    t72 = d__1 * d__1;
		    t73 = 1 / t72;
		    t76 = t63 * -.1835293601259385 - t67 * .05399117105794445;
		    t79 = t73 * .0843 * t76 - t71;
		    t82 = t21 * t18;
		    t85 = t25 * t18;
		    t88 = t82 * -1.333333333333333 * t62 + t85 * 
			    1.333333333333333 * t62;
		    t93 = t66 * t31;
		    t94 = t93 * t62;
		    t96 = t2 * t1;
		    t103 = t1 * .005183333333333333 + t94 * 
			    2.688185460564067e-4 + t96 * 2.688185460564067e-4 
			    - t67 * .001406127779371973;
		    d__1 = t71 + t79 * 1.923661050931536 * t27 + t16 * 
			    1.923661050931536 * t88;
		    d__2 = t1 * -.01036666666666667 - t94 * 
			    4.135669939329333e-4 - t96 * 4.135669939329333e-4 
			    + t67 * .002398688564811013 + t103 * 
			    1.923661050931536 * t27 + t39 * 1.923661050931536 
			    * t88;
		    t109 = piecewise_(&t4, &d__1, &d__2);
		    t110 = rho * t109;
		    vrhoa[i__] = rho * t53 + t43 + t110;
		    t111 = -t48;
		    d__1 = t16 * 1.923661050931536 * t111;
		    d__2 = t39 * 1.923661050931536 * t111;
		    t116 = piecewise_(&t4, &d__1, &d__2);
		    vrhob[i__] = rho * t116 + t43 + t110;
/* Computing 2nd power */
		    d__1 = t21;
		    t118 = d__1 * d__1;
		    t119 = 1 / t118;
/* Computing 2nd power */
		    d__1 = t25;
		    t122 = d__1 * d__1;
		    t123 = 1 / t122;
		    t126 = t119 * .4444444444444444 * t62 + t123 * 
			    .4444444444444444 * t62;
		    d__1 = t16 * 1.923661050931536 * t126;
		    d__2 = t39 * 1.923661050931536 * t126;
		    t131 = piecewise_(&t4, &d__1, &d__2);
		    t132 = rho * t131;
		    t138 = 1 / t61 / rho;
		    t148 = t119 * -.4444444444444444 * t18 * t138 - t21 * 
			    1.333333333333333 * t62 - t123 * 
			    .4444444444444444 * t18 * t138 + t25 * 
			    1.333333333333333 * t62;
		    d__1 = t79 * 1.923661050931536 * t48 + t16 * 
			    1.923661050931536 * t148;
		    d__2 = t103 * 1.923661050931536 * t48 + t39 * 
			    1.923661050931536 * t148;
		    t157 = piecewise_(&t4, &d__1, &d__2);
		    t158 = rho * t157;
		    t160 = t109 * 2.;
/* Computing 2nd power */
		    d__1 = t69;
		    t163 = d__1 * d__1;
		    t165 = .2846 / t55 / t8 * t163;
/* Computing 2nd power */
		    d__1 = t61;
		    t168 = d__1 * d__1;
		    t169 = 1 / t168;
		    t170 = 1 / t59 / t1 * t169;
		    t172 = t60 * t138;
		    t175 = 1 / t65 / t1;
		    t176 = t175 * t169;
		    t178 = t66 * t138;
		    t182 = t56 * .1423 * (t170 * -.1151789710300888 + t172 * 
			    .2764295304722132 - t176 * .04596107859241333 + 
			    t178 * .13788323577724);
/* Computing 2nd power */
		    d__1 = t76;
		    t185 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t200 = d__1 * d__1;
		    t211 = t119 * .4444444444444444 * t200 * t169 + t82 * 
			    2.666666666666667 * t138 + t123 * 
			    .4444444444444444 * t200 * t169 - t85 * 
			    2.666666666666667 * t138;
		    t217 = t175 * t31 * t169;
		    t220 = t93 * t138;
		    t222 = t2 * t62;
		    d__1 = -t165 + t182 + (-.1686 / t72 / t13 * t185 + t73 * 
			    .0843 * (t170 * -.1529411334382821 + t172 * 
			    .367058720251877 - t176 * .03599411403862963 + 
			    t178 * .1079823421158889) + t165 - t182) * 
			    1.923661050931536 * t27 + t79 * 3.847322101863073 
			    * t88 + t16 * 1.923661050931536 * t211;
		    d__2 = t62 * .01036666666666667 - t217 * 
			    2.757113292886222e-4 - t178 * .004521665800333405 
			    + t220 * 8.271339878658667e-4 + t222 * 
			    4.135669939329333e-4 + t176 * .001599125709874009 
			    + (t62 * -.005183333333333333 + t217 * 
			    1.792123640376044e-4 + t178 * .002633043194706342 
			    - t220 * 5.376370921128133e-4 - t222 * 
			    2.688185460564067e-4 - t176 * 
			    9.374185195813156e-4) * 1.923661050931536 * t27 + 
			    t103 * 3.847322101863073 * t88 + t39 * 
			    1.923661050931536 * t211;
		    t239 = piecewise_(&t4, &d__1, &d__2);
		    t240 = rho * t239;
		    v2rhoa2[i__] = t132 + t53 * 2. + t158 * 2. + t160 + t240;
		    t244 = -t148;
		    d__1 = t79 * 1.923661050931536 * t111 + t16 * 
			    1.923661050931536 * t244;
		    d__2 = t103 * 1.923661050931536 * t111 + t39 * 
			    1.923661050931536 * t244;
		    t253 = piecewise_(&t4, &d__1, &d__2);
		    t254 = rho * t253;
		    v2rhob2[i__] = t132 + t116 * 2. + t254 * 2. + t160 + t240;
		    t256 = -t126;
		    d__1 = t16 * 1.923661050931536 * t256;
		    d__2 = t39 * 1.923661050931536 * t256;
		    t261 = piecewise_(&t4, &d__1, &d__2);
		    v2rhoab[i__] = rho * t261 + t116 + t254 + t53 + t158 + 
			    t160 + t240;
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
	    }
/* rho */
	}
    }
/* ideriv */
    return 0;
} /* uks_c_pz81__ */

/* Subroutine */ int rks_c_pz81__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    logical L__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t1, t2, t3;
    static logical t4;
    static doublereal t5, t8, t10, t11, t20, t22, t23, t30, t31, t17, t18, 
	    t26, t27, t32, t45, t13, t19, t21, t24, t25, t34, t38, t53, t59, 
	    t61, t68, t73, t74, t78, t82, t83, t85, t102, t107, rho;
    extern doublereal piecewise_(logical *, doublereal *, doublereal *);


/*     J.P. Perdew, and A. Zunger */
/*     Self-interaction correction to density-functional approximations */
/*     for many-electron systems */
/*     Phys. Rev. B23 (1981) 5048-5079 */


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
		t2 = pow_dd(&t1, &c_b2);
		t3 = t2 * .6203504908994;
		t5 = pow_dd(&t1, &c_b3);
		t11 = log(t3);
		L__1 = 1. <= t3;
		d__1 = -.1423 / (t5 * .8292885914166397 + 1. + t2 * 
			.20682485366586);
		d__2 = t11 * .0311 - .048 + t2 * .0012407009817988 * t11 - t2 
			* .00719606569443304;
		t17 = piecewise_(&L__1, &d__1, &d__2);
		zk[i__] = rho * t17;
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
		t2 = pow_dd(&t1, &c_b2);
		t3 = t2 * .6203504908994;
		t4 = 1. <= t3;
		t5 = pow_dd(&t1, &c_b3);
		t8 = t5 * .8292885914166397 + 1. + t2 * .20682485366586;
		t11 = log(t3);
		d__1 = -.1423 / t8;
		d__2 = t11 * .0311 - .048 + t2 * .0012407009817988 * t11 - t2 
			* .00719606569443304;
		t17 = piecewise_(&t4, &d__1, &d__2);
		zk[i__] = rho * t17;
		t18 = piecewise_(&t4, &c_b31, &c_b31);
/* Computing 2nd power */
		d__1 = t8;
		t20 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t5;
		t22 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t22;
		t23 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = rho;
		t26 = d__1 * d__1;
		t27 = 1 / t26;
/* Computing 2nd power */
		d__1 = t2;
		t30 = d__1 * d__1;
		t31 = 1 / t30;
		t32 = t31 * t27;
		d__1 = .1423 / t20 * (-.1382147652361066 / t23 / t5 * t27 - 
			t32 * .06894161788861999);
		d__2 = t1 * -.01036666666666667 - t31 * 4.135669939329333e-4 *
			 t11 * t27 - t2 * 4.135669939329333e-4 * t1 + t32 * 
			.002398688564811013;
		t45 = piecewise_(&t4, &d__1, &d__2);
		vrhoa[i__] = rho * t18 + t17 + rho * t45;
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
		t2 = pow_dd(&t1, &c_b2);
		t3 = t2 * .6203504908994;
		t4 = 1. <= t3;
		t5 = pow_dd(&t1, &c_b3);
		t8 = t5 * .8292885914166397 + 1. + t2 * .20682485366586;
		t10 = .1423 / t8;
		t11 = log(t3);
		t13 = t2 * t11;
		d__1 = -t10;
		d__2 = t11 * .0311 - .048 + t13 * .0012407009817988 - t2 * 
			.00719606569443304;
		t17 = piecewise_(&t4, &d__1, &d__2);
		zk[i__] = rho * t17;
		t18 = piecewise_(&t4, &c_b31, &c_b31);
		t19 = rho * t18;
/* Computing 2nd power */
		d__1 = t8;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
/* Computing 2nd power */
		d__1 = t5;
		t22 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t22;
		t23 = d__1 * d__1;
		t24 = t23 * t5;
		t25 = 1 / t24;
/* Computing 2nd power */
		d__1 = rho;
		t26 = d__1 * d__1;
		t27 = 1 / t26;
/* Computing 2nd power */
		d__1 = t2;
		t30 = d__1 * d__1;
		t31 = 1 / t30;
		t32 = t31 * t27;
		t34 = t25 * -.1382147652361066 * t27 - t32 * 
			.06894161788861999;
		t38 = t31 * t11;
		d__1 = t21 * .1423 * t34;
		d__2 = t1 * -.01036666666666667 - t38 * 4.135669939329333e-4 *
			 t27 - t2 * 4.135669939329333e-4 * t1 + t32 * 
			.002398688564811013;
		t45 = piecewise_(&t4, &d__1, &d__2);
		vrhoa[i__] = t19 + t17 + rho * t45;
		t53 = (-.0843 / (t5 * 1.101176160755631 + 1. + t2 * 
			.1619735131738333) + t10) * t27;
		t59 = (t11 * -.01555 + .0211 - t13 * 8.0645563816922e-4 + t2 *
			 .00421838333811592) * t27;
		d__1 = t53 * 1.709920934161366;
		d__2 = t59 * 1.709920934161366;
		t61 = piecewise_(&t4, &d__1, &d__2);
/* Computing 2nd power */
		d__1 = t34;
		t68 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t26;
		t73 = d__1 * d__1;
		t74 = 1 / t73;
		t78 = 1 / t26 / rho;
		t82 = 1 / t30 / t1;
		t83 = t82 * t74;
		t85 = t31 * t78;
		d__1 = -.2846 / t20 / t8 * t68 + t21 * .1423 * (
			-.1151789710300888 / t24 / t1 * t74 + t25 * 
			.2764295304722132 * t78 - t83 * .04596107859241333 + 
			t85 * .13788323577724);
		d__2 = t27 * .01036666666666667 - t82 * 2.757113292886222e-4 *
			 t11 * t74 - t85 * .004521665800333405 + t38 * 
			8.271339878658667e-4 * t78 + t2 * 
			4.135669939329333e-4 * t27 + t83 * 
			.001599125709874009;
		t102 = piecewise_(&t4, &d__1, &d__2);
		d__1 = t53 * -1.709920934161366;
		d__2 = t59 * -1.709920934161366;
		t107 = piecewise_(&t4, &d__1, &d__2);
		v2rhoa2[i__] = rho * t61 + t18 * 4. + t19 * 4. + t45 * 4. + 
			rho * 2. * t102 + rho * t107;
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
} /* rks_c_pz81__ */


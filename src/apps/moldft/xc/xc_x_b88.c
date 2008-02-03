/* xc_x_b88.f -- translated by f2c (version 20050501).
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

/* :X_B88subrstart */
/*    Generated: Wed Jan 29 09:08:35 GMT 2003 */
/* Subroutine */ int uks_x_b88__(integer *ideriv, integer *npt, doublereal *
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
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t21, t12, t23,
	     t24, t25, t13, t18, t19, t17, t29, t34, t37, t38, t14, t15, t22, 
	    t28, t33, t35, t39, t40, t45, t50, t53, t54, t62, t64, t68, t69, 
	    t74, t79, t82, t83, t20, t36, t41, t42, t60, t67, t75, t81, t87, 
	    t97, t44, t47, t52, t57, t58, t65, t73, t76, t86, t92, t98, t99, 
	    t104, t110, t111, t117, t124, t125, t132, t138, t144, t154, t161, 
	    t162, t169, t175, t181, t192, t218, t242, t266, rho, rhoa, rhob, 
	    sigma, sigmaaa, sigmaab, sigmabb;


/*     A.D. Becke */
/*     Density-functional exchange-energy approximation with correct */
/*     asymptotic behaviour */
/*     Phys. Rev. A38 (1988) 3098-3100 */


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
		    t5 = 1 / t3;
		    t7 = sqrt(sigmabb);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    zk[i__] = t3 * -.9305257363491 - t5 * .0042 * sigmabb / (
			    t8 * .0252 * t9 + 1.);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t5 = 1 / t3;
		    t7 = sqrt(sigmaaa);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    zk[i__] = t3 * -.9305257363491 - t5 * .0042 * sigmaaa / (
			    t8 * .0252 * t9 + 1.);
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
		    t7 = 1 / t5;
		    t9 = sqrt(sigmaaa);
		    t10 = t9 * t7;
/* Computing 2nd power */
		    d__1 = t10;
		    t11 = log(t10 + sqrt(d__1 * d__1 + 1));
		    t18 = pow_dd(&rhob, &c_b2);
		    t19 = t18 * rhob;
		    t21 = 1 / t19;
		    t23 = sqrt(sigmabb);
		    t24 = t23 * t21;
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = log(t24 + sqrt(d__1 * d__1 + 1));
		    zk[i__] = t5 * -.9305257363491 - t7 * .0042 * sigmaaa / (
			    t10 * .0252 * t11 + 1.) - t19 * .9305257363491 - 
			    t21 * .0042 * sigmabb / (t24 * .0252 * t25 + 1.);
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
		    t5 = 1 / t3;
		    t6 = t5 * sigmabb;
		    t7 = sqrt(sigmabb);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    t12 = t8 * .0252 * t9 + 1.;
		    t13 = 1 / t12;
		    zk[i__] = t3 * -.9305257363491 - t6 * .0042 * t13;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = rhob;
		    t17 = d__1 * d__1;
		    t19 = 1 / t2 / t17;
/* Computing 2nd power */
		    d__1 = t12;
		    t23 = d__1 * d__1;
		    t24 = 1 / t23;
/* Computing 2nd power */
		    d__1 = t2;
		    t29 = d__1 * d__1;
		    t34 = 1 / t29 / t17;
		    t37 = sqrt(sigmabb * t34 + 1.);
		    t38 = 1 / t37;
		    vrhob[i__] = t2 * -1.2407009817988 + t19 * .0056 * 
			    sigmabb * t13 + t6 * .0042 * t24 * (t7 * -.0336 * 
			    t19 * t9 - sigmabb * .0336 / t29 / t17 / rhob * 
			    t38);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = t5 * -.0042 * t13 + t6 * .0042 * t24 * (
			    .0126 / t7 * t5 * t9 + t34 * .0126 * t38);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t5 = 1 / t3;
		    t6 = t5 * sigmaaa;
		    t7 = sqrt(sigmaaa);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    t12 = t8 * .0252 * t9 + 1.;
		    t13 = 1 / t12;
		    zk[i__] = t3 * -.9305257363491 - t6 * .0042 * t13;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t17 = d__1 * d__1;
		    t19 = 1 / t2 / t17;
/* Computing 2nd power */
		    d__1 = t12;
		    t23 = d__1 * d__1;
		    t24 = 1 / t23;
/* Computing 2nd power */
		    d__1 = t2;
		    t29 = d__1 * d__1;
		    t34 = 1 / t29 / t17;
		    t37 = sqrt(sigmaaa * t34 + 1.);
		    t38 = 1 / t37;
		    vrhoa[i__] = t2 * -1.2407009817988 + t19 * .0056 * 
			    sigmaaa * t13 + t6 * .0042 * t24 * (t7 * -.0336 * 
			    t19 * t9 - sigmaaa * .0336 / t29 / t17 / rhoa * 
			    t38);
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = t5 * -.0042 * t13 + t6 * .0042 * t24 * (
			    .0126 / t7 * t5 * t9 + t34 * .0126 * t38);
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
		    t7 = 1 / t5;
		    t8 = t7 * sigmaaa;
		    t9 = sqrt(sigmaaa);
		    t10 = t9 * t7;
/* Computing 2nd power */
		    d__1 = t10;
		    t11 = log(t10 + sqrt(d__1 * d__1 + 1));
		    t14 = t10 * .0252 * t11 + 1.;
		    t15 = 1 / t14;
		    t18 = pow_dd(&rhob, &c_b2);
		    t19 = t18 * rhob;
		    t21 = 1 / t19;
		    t22 = t21 * sigmabb;
		    t23 = sqrt(sigmabb);
		    t24 = t23 * t21;
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = log(t24 + sqrt(d__1 * d__1 + 1));
		    t28 = t24 * .0252 * t25 + 1.;
		    t29 = 1 / t28;
		    zk[i__] = t5 * -.9305257363491 - t8 * .0042 * t15 - t19 * 
			    .9305257363491 - t22 * .0042 * t29;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t33 = d__1 * d__1;
		    t35 = 1 / t4 / t33;
/* Computing 2nd power */
		    d__1 = t14;
		    t39 = d__1 * d__1;
		    t40 = 1 / t39;
/* Computing 2nd power */
		    d__1 = t4;
		    t45 = d__1 * d__1;
		    t50 = 1 / t45 / t33;
		    t53 = sqrt(sigmaaa * t50 + 1.);
		    t54 = 1 / t53;
		    vrhoa[i__] = t4 * -1.2407009817988 + t35 * .0056 * 
			    sigmaaa * t15 + t8 * .0042 * t40 * (t9 * -.0336 * 
			    t35 * t11 - sigmaaa * .0336 / t45 / t33 / rhoa * 
			    t54);
/* Computing 2nd power */
		    d__1 = rhob;
		    t62 = d__1 * d__1;
		    t64 = 1 / t18 / t62;
/* Computing 2nd power */
		    d__1 = t28;
		    t68 = d__1 * d__1;
		    t69 = 1 / t68;
/* Computing 2nd power */
		    d__1 = t18;
		    t74 = d__1 * d__1;
		    t79 = 1 / t74 / t62;
		    t82 = sqrt(sigmabb * t79 + 1.);
		    t83 = 1 / t82;
		    vrhob[i__] = t18 * -1.2407009817988 + t64 * .0056 * 
			    sigmabb * t29 + t22 * .0042 * t69 * (t23 * -.0336 
			    * t64 * t25 - sigmabb * .0336 / t74 / t62 / rhob *
			     t83);
		    vsigmaaa[i__] = t7 * -.0042 * t15 + t8 * .0042 * t40 * (
			    .0126 / t9 * t7 * t11 + t50 * .0126 * t54);
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = t21 * -.0042 * t29 + t22 * .0042 * t69 * (
			    .0126 / t23 * t21 * t25 + t79 * .0126 * t83);
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
		    t5 = 1 / t3;
		    t6 = t5 * sigmabb;
		    t7 = sqrt(sigmabb);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    t12 = t8 * .0252 * t9 + 1.;
		    t13 = 1 / t12;
		    zk[i__] = t3 * -.9305257363491 - t6 * .0042 * t13;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = rhob;
		    t17 = d__1 * d__1;
		    t19 = 1 / t2 / t17;
		    t20 = t19 * sigmabb;
/* Computing 2nd power */
		    d__1 = t12;
		    t23 = d__1 * d__1;
		    t24 = 1 / t23;
		    t28 = t17 * rhob;
/* Computing 2nd power */
		    d__1 = t2;
		    t29 = d__1 * d__1;
		    t34 = 1 / t29 / t17;
		    t36 = sigmabb * t34 + 1.;
		    t37 = sqrt(t36);
		    t38 = 1 / t37;
		    t41 = t7 * -.0336 * t19 * t9 - sigmabb * .0336 / t29 / 
			    t28 * t38;
		    t42 = t24 * t41;
		    vrhob[i__] = t2 * -1.2407009817988 + t20 * .0056 * t13 + 
			    t6 * .0042 * t42;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t53 = .0126 / t7 * t5 * t9 + t34 * .0126 * t38;
		    vsigmabb[i__] = t5 * -.0042 * t13 + t6 * .0042 * t24 * 
			    t53;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t60 = 1 / t2 / t28;
		    t67 = 1 / t23 / t12;
/* Computing 2nd power */
		    d__1 = t41;
		    t68 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t17;
		    t75 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t81 = d__1 * d__1;
		    t87 = 1 / t37 / t36;
		    v2rhob2[i__] = -.4135669939329333 / t29 - t60 * 
			    .01306666666666667 * sigmabb * t13 - t20 * .0112 *
			     t42 - t6 * .0084 * t67 * t68 + t6 * .0042 * t24 *
			     (t7 * .0784 * t60 * t9 + sigmabb * .168 / t29 / 
			    t75 * t38 - t81 * .0448 / t2 / t75 / t28 * t87);
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t53;
		    t97 = d__1 * d__1;
		    v2sigmabb2[i__] = t5 * .0084 * t24 * t53 - t6 * .0084 * 
			    t67 * t97 + t6 * .0042 * t24 * (-.0063 / t7 / 
			    sigmabb * t5 * t9 + .0063 / sigmabb * t34 * t38 - 
			    .0063 / t2 / t75 / rhob * t87);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
		    t5 = 1 / t3;
		    t6 = t5 * sigmaaa;
		    t7 = sqrt(sigmaaa);
		    t8 = t7 * t5;
/* Computing 2nd power */
		    d__1 = t8;
		    t9 = log(t8 + sqrt(d__1 * d__1 + 1));
		    t12 = t8 * .0252 * t9 + 1.;
		    t13 = 1 / t12;
		    zk[i__] = t3 * -.9305257363491 - t6 * .0042 * t13;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t17 = d__1 * d__1;
		    t19 = 1 / t2 / t17;
		    t20 = t19 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t12;
		    t23 = d__1 * d__1;
		    t24 = 1 / t23;
		    t28 = t17 * rhoa;
/* Computing 2nd power */
		    d__1 = t2;
		    t29 = d__1 * d__1;
		    t34 = 1 / t29 / t17;
		    t36 = sigmaaa * t34 + 1.;
		    t37 = sqrt(t36);
		    t38 = 1 / t37;
		    t41 = t7 * -.0336 * t19 * t9 - sigmaaa * .0336 / t29 / 
			    t28 * t38;
		    t42 = t24 * t41;
		    vrhoa[i__] = t2 * -1.2407009817988 + t20 * .0056 * t13 + 
			    t6 * .0042 * t42;
		    vrhob[i__] = 0.;
		    t53 = .0126 / t7 * t5 * t9 + t34 * .0126 * t38;
		    vsigmaaa[i__] = t5 * -.0042 * t13 + t6 * .0042 * t24 * 
			    t53;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t60 = 1 / t2 / t28;
		    t67 = 1 / t23 / t12;
/* Computing 2nd power */
		    d__1 = t41;
		    t68 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t17;
		    t75 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t81 = d__1 * d__1;
		    t87 = 1 / t37 / t36;
		    v2rhoa2[i__] = -.4135669939329333 / t29 - t60 * 
			    .01306666666666667 * sigmaaa * t13 - t20 * .0112 *
			     t42 - t6 * .0084 * t67 * t68 + t6 * .0042 * t24 *
			     (t7 * .0784 * t60 * t9 + sigmaaa * .168 / t29 / 
			    t75 * t38 - t81 * .0448 / t2 / t75 / t28 * t87);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t53;
		    t97 = d__1 * d__1;
		    v2sigmaaa2[i__] = t5 * .0084 * t24 * t53 - t6 * .0084 * 
			    t67 * t97 + t6 * .0042 * t24 * (-.0063 / t7 / 
			    sigmaaa * t5 * t9 + .0063 / sigmaaa * t34 * t38 - 
			    .0063 / t2 / t75 / rhoa * t87);
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
		    t7 = 1 / t5;
		    t8 = t7 * sigmaaa;
		    t9 = sqrt(sigmaaa);
		    t10 = t9 * t7;
/* Computing 2nd power */
		    d__1 = t10;
		    t11 = log(t10 + sqrt(d__1 * d__1 + 1));
		    t14 = t10 * .0252 * t11 + 1.;
		    t15 = 1 / t14;
		    t18 = pow_dd(&rhob, &c_b2);
		    t19 = t18 * rhob;
		    t21 = 1 / t19;
		    t22 = t21 * sigmabb;
		    t23 = sqrt(sigmabb);
		    t24 = t23 * t21;
/* Computing 2nd power */
		    d__1 = t24;
		    t25 = log(t24 + sqrt(d__1 * d__1 + 1));
		    t28 = t24 * .0252 * t25 + 1.;
		    t29 = 1 / t28;
		    zk[i__] = t5 * -.9305257363491 - t8 * .0042 * t15 - t19 * 
			    .9305257363491 - t22 * .0042 * t29;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t33 = d__1 * d__1;
		    t35 = 1 / t4 / t33;
		    t36 = t35 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t14;
		    t39 = d__1 * d__1;
		    t40 = 1 / t39;
		    t44 = t33 * rhoa;
/* Computing 2nd power */
		    d__1 = t4;
		    t45 = d__1 * d__1;
		    t47 = 1 / t45 / t44;
		    t50 = 1 / t45 / t33;
		    t52 = sigmaaa * t50 + 1.;
		    t53 = sqrt(t52);
		    t54 = 1 / t53;
		    t57 = t9 * -.0336 * t35 * t11 - sigmaaa * .0336 * t47 * 
			    t54;
		    t58 = t40 * t57;
		    vrhoa[i__] = t4 * -1.2407009817988 + t36 * .0056 * t15 + 
			    t8 * .0042 * t58;
/* Computing 2nd power */
		    d__1 = rhob;
		    t62 = d__1 * d__1;
		    t64 = 1 / t18 / t62;
		    t65 = t64 * sigmabb;
/* Computing 2nd power */
		    d__1 = t28;
		    t68 = d__1 * d__1;
		    t69 = 1 / t68;
		    t73 = t62 * rhob;
/* Computing 2nd power */
		    d__1 = t18;
		    t74 = d__1 * d__1;
		    t76 = 1 / t74 / t73;
		    t79 = 1 / t74 / t62;
		    t81 = sigmabb * t79 + 1.;
		    t82 = sqrt(t81);
		    t83 = 1 / t82;
		    t86 = t23 * -.0336 * t64 * t25 - sigmabb * .0336 * t76 * 
			    t83;
		    t87 = t69 * t86;
		    vrhob[i__] = t18 * -1.2407009817988 + t65 * .0056 * t29 + 
			    t22 * .0042 * t87;
		    t92 = 1 / t9;
		    t98 = t92 * .0126 * t7 * t11 + t50 * .0126 * t54;
		    t99 = t40 * t98;
		    vsigmaaa[i__] = t7 * -.0042 * t15 + t8 * .0042 * t99;
		    vsigmaab[i__] = 0.;
		    t104 = 1 / t23;
		    t110 = t104 * .0126 * t21 * t25 + t79 * .0126 * t83;
		    t111 = t69 * t110;
		    vsigmabb[i__] = t21 * -.0042 * t29 + t22 * .0042 * t111;
		    t117 = 1 / t4 / t44;
		    t124 = 1 / t39 / t14;
/* Computing 2nd power */
		    d__1 = t57;
		    t125 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t33;
		    t132 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t138 = d__1 * d__1;
		    t144 = 1 / t53 / t52;
		    v2rhoa2[i__] = -.4135669939329333 / t45 - t117 * 
			    .01306666666666667 * sigmaaa * t15 - t36 * .0112 *
			     t58 - t8 * .0084 * t124 * t125 + t8 * .0042 * 
			    t40 * (t9 * .0784 * t117 * t11 + sigmaaa * .168 / 
			    t45 / t132 * t54 - t138 * .0448 / t4 / t132 / t44 
			    * t144);
		    t154 = 1 / t18 / t73;
		    t161 = 1 / t68 / t28;
/* Computing 2nd power */
		    d__1 = t86;
		    t162 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t62;
		    t169 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t175 = d__1 * d__1;
		    t181 = 1 / t82 / t81;
		    v2rhob2[i__] = -.4135669939329333 / t74 - t154 * 
			    .01306666666666667 * sigmabb * t29 - t65 * .0112 *
			     t87 - t22 * .0084 * t161 * t162 + t22 * .0042 * 
			    t69 * (t23 * .0784 * t154 * t25 + sigmabb * .168 /
			     t74 / t169 * t83 - t175 * .0448 / t18 / t169 / 
			    t73 * t181);
		    v2rhoab[i__] = 0.;
		    t192 = t7 * t40;
		    v2rhoasigmaaa[i__] = t35 * .0056 * t15 - t36 * .0056 * 
			    t99 + t192 * .0042 * t57 - t8 * .0084 * t124 * 
			    t57 * t98 + t8 * .0042 * t40 * (t92 * -.0168 * 
			    t35 * t11 - t47 * .0504 * t54 + sigmaaa * .0168 / 
			    t4 / t132 / t33 * t144);
		    v2rhoasigmaab[i__] = 0.;
		    v2rhoasigmabb[i__] = 0.;
		    v2rhobsigmaaa[i__] = 0.;
		    v2rhobsigmaab[i__] = 0.;
		    t218 = t21 * t69;
		    v2rhobsigmabb[i__] = t64 * .0056 * t29 - t65 * .0056 * 
			    t111 + t218 * .0042 * t86 - t22 * .0084 * t161 * 
			    t86 * t110 + t22 * .0042 * t69 * (t104 * -.0168 * 
			    t64 * t25 - t76 * .0504 * t83 + sigmabb * .0168 / 
			    t18 / t169 / t62 * t181);
/* Computing 2nd power */
		    d__1 = t98;
		    t242 = d__1 * d__1;
		    v2sigmaaa2[i__] = t192 * .0084 * t98 - t8 * .0084 * t124 *
			     t242 + t8 * .0042 * t40 * (-.0063 / t9 / sigmaaa 
			    * t7 * t11 + .0063 / sigmaaa * t50 * t54 - .0063 /
			     t4 / t132 / rhoa * t144);
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t110;
		    t266 = d__1 * d__1;
		    v2sigmabb2[i__] = t218 * .0084 * t110 - t22 * .0084 * 
			    t161 * t266 + t22 * .0042 * t69 * (-.0063 / t23 / 
			    sigmabb * t21 * t25 + .0063 / sigmabb * t79 * t83 
			    - .0063 / t18 / t169 / rhob * t181);
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
} /* uks_x_b88__ */

/* Subroutine */ int rks_x_b88__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t5, t6, t7, t8, t10, t20, t30, t13, t14, t24, 
	    t25, t35, t18, t40, t21, t39, t29, t32, t38, t43, t44, t49, t55, 
	    t56, t62, t69, t70, t77, t83, t89, t100, t124, rho, sigma;


/*     A.D. Becke */
/*     Density-functional exchange-energy approximation with correct */
/*     asymptotic behaviour */
/*     Phys. Rev. A38 (1988) 3098-3100 */


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
		t5 = 1 / t3;
		t7 = sqrt(sigma);
		t8 = t7 * t5;
/* Computing 2nd power */
		d__1 = t8;
		t10 = log(t8 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		zk[i__] = t3 * -.7385587663820224 - t5 * .005291668409558467 *
			 sigma / (t8 * .0317500104573508 * t10 + 1.);
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
		t5 = 1 / t3;
		t6 = t5 * sigma;
		t7 = sqrt(sigma);
		t8 = t7 * t5;
/* Computing 2nd power */
		d__1 = t8;
		t10 = log(t8 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		t13 = t8 * .0317500104573508 * t10 + 1.;
		t14 = 1 / t13;
		zk[i__] = t3 * -.7385587663820224 - t6 * .005291668409558467 *
			 t14;
/* Computing 2nd power */
		d__1 = rho;
		t18 = d__1 * d__1;
		t20 = 1 / t2 / t18;
/* Computing 2nd power */
		d__1 = t13;
		t24 = d__1 * d__1;
		t25 = 1 / t24;
/* Computing 2nd power */
		d__1 = t2;
		t30 = d__1 * d__1;
		t35 = 1 / t30 / t18;
		t39 = sqrt(sigma * 1.587401051968199 * t35 + 1.);
		t40 = 1 / t39;
		vrhoa[i__] = t2 * -.9847450218426965 + t20 * 
			.00705555787941129 * sigma * t14 + t6 * 
			.002645834204779234 * t25 * (t7 * -.08466669455293548 
			* t20 * t10 - sigma * .106673350692263 / t30 / t18 / 
			rho * t40);
		vsigmaaa[i__] = t5 * -.02116667363823387 * t14 + t6 * 
			.005291668409558467 * t25 * (.06350002091470161 / t7 *
			 t5 * t10 + t35 * .08000501301919725 * t40);
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
		t5 = 1 / t3;
		t6 = t5 * sigma;
		t7 = sqrt(sigma);
		t8 = t7 * t5;
/* Computing 2nd power */
		d__1 = t8;
		t10 = log(t8 * 1.259921049894873 + sqrt(d__1 * d__1 * 
			1.587401051968199 + 1));
		t13 = t8 * .0317500104573508 * t10 + 1.;
		t14 = 1 / t13;
		zk[i__] = t3 * -.7385587663820224 - t6 * .005291668409558467 *
			 t14;
/* Computing 2nd power */
		d__1 = rho;
		t18 = d__1 * d__1;
		t20 = 1 / t2 / t18;
		t21 = t20 * sigma;
/* Computing 2nd power */
		d__1 = t13;
		t24 = d__1 * d__1;
		t25 = 1 / t24;
		t29 = t18 * rho;
/* Computing 2nd power */
		d__1 = t2;
		t30 = d__1 * d__1;
		t32 = 1 / t30 / t29;
		t35 = 1 / t30 / t18;
		t38 = sigma * 1.587401051968199 * t35 + 1.;
		t39 = sqrt(t38);
		t40 = 1 / t39;
		t43 = t7 * -.08466669455293548 * t20 * t10 - sigma * 
			.106673350692263 * t32 * t40;
		t44 = t25 * t43;
		vrhoa[i__] = t2 * -.9847450218426965 + t21 * 
			.00705555787941129 * t14 + t6 * .002645834204779234 * 
			t44;
		t49 = 1 / t7;
		t55 = t49 * .06350002091470161 * t5 * t10 + t35 * 
			.08000501301919725 * t40;
		t56 = t25 * t55;
		vsigmaaa[i__] = t5 * -.02116667363823387 * t14 + t6 * 
			.005291668409558467 * t56;
		t62 = 1 / t2 / t29;
		t69 = 1 / t24 / t13;
/* Computing 2nd power */
		d__1 = t43;
		t70 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t18;
		t77 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = sigma;
		t83 = d__1 * d__1;
		t89 = 1 / t39 / t38;
		v2rhoa2[i__] = -.6564966812284644 / t30 - t62 * 
			.03292593677058602 * sigma * t14 - t21 * 
			.01411111575882258 * t44 - t6 * .005291668409558467 * 
			t69 * t70 + t6 * .002645834204779234 * t25 * (t7 * 
			.3951112412470322 * t62 * t10 + sigma * 
			1.06673350692263 / t30 / t77 * t40 - t83 * 
			.4515557042823225 / t2 / t77 / t29 * t89);
		t100 = t5 * t25;
		v2rhoasigmaaa[i__] = t20 * .02822223151764516 * t14 - t21 * 
			.00705555787941129 * t56 + t100 * .01058333681911693 *
			 t43 - t6 * .005291668409558467 * t69 * t43 * t55 + 
			t6 * .002645834204779234 * t25 * (t49 * 
			-.169333389105871 * t20 * t10 - t32 * 
			.640040104153578 * t40 + sigma * .3386667782117419 / 
			t2 / t77 / t18 * t89);
/* Computing 2nd power */
		d__1 = t55;
		t124 = d__1 * d__1;
		v2sigmaaa2[i__] = t100 * .04233334727646774 * t55 - t6 * 
			.01058333681911693 * t69 * t124 + t6 * 
			.005291668409558467 * t25 * (-.1270000418294032 / t7 /
			 sigma * t5 * t10 + .1600100260383945 / sigma * t35 * 
			t40 - .2540000836588064 / t2 / t77 / rho * t89);
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
} /* rks_x_b88__ */


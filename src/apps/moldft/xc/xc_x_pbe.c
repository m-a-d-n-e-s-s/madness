/* xc_x_pbe.f -- translated by f2c (version 20050501).
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

/* :X_PBEsubrstart */
/*    Generated: Wed Sep  3 12:48:32 GMT 2003 */
/* Subroutine */ int uks_x_pbe__(integer *ideriv, integer *npt, doublereal *
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
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t6, t7, t10, t20, t21, t13, t12, t15, 
	    t41, t26, t18, t19, t29, t36, t37, t45, t46, t43, t38, t47, t65, 
	    t69, t71, t83, t87, t89, rho, rhoa, rhob, sigma, sigmaaa, sigmaab,
	     sigmabb;


/*     J.P. Perdew, K. Burke, and M. Ernzerhof */
/*     Generalized gradient approximation made simple */
/*     Phys. Rev. Lett. 77 (1996) 3865-3868 */


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
/* Computing 2nd power */
		    d__1 = rhob;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhob * (1.804 - .804 / (
			    sigmabb * .00449276922095889 / t5 / t4 + 1.));
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhoa * (1.804 - .804 / (
			    sigmaaa * .00449276922095889 / t5 / t4 + 1.));
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t6 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t7 = d__1 * d__1;
		    t18 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t21 = d__1 * d__1;
		    zk[i__] = t4 * -.9305257363491 * rhoa * (1.804 - .804 / (
			    sigmaaa * .00449276922095889 / t7 / t6 + 1.)) - 
			    t18 * .9305257363491 * rhob * (1.804 - .804 / (
			    sigmabb * .00449276922095889 / t21 / t20 + 1.));
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
/* Computing 2nd power */
		    d__1 = rhob;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    t10 = sigmabb * .00449276922095889 / t5 / t4 + 1.;
		    t13 = 1.804 - .804 / t10;
		    zk[i__] = t3 * -.9305257363491 * t13;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    vrhob[i__] = t2 * -1.2407009817988 * t13 + 
			    .008963286558970112 / t2 / t4 * t21 * sigmabb;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = -.003361232459613792 / t3 * t21;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    t10 = sigmaaa * .00449276922095889 / t5 / t4 + 1.;
		    t13 = 1.804 - .804 / t10;
		    zk[i__] = t3 * -.9305257363491 * t13;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    vrhoa[i__] = t2 * -1.2407009817988 * t13 + 
			    .008963286558970112 / t2 / t4 * t21 * sigmaaa;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = -.003361232459613792 / t3 * t21;
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t6 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t7 = d__1 * d__1;
		    t12 = sigmaaa * .00449276922095889 / t7 / t6 + 1.;
		    t15 = 1.804 - .804 / t12;
		    t18 = pow_dd(&rhob, &c_b2);
		    t19 = t18 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t21 = d__1 * d__1;
		    t26 = sigmabb * .00449276922095889 / t21 / t20 + 1.;
		    t29 = 1.804 - .804 / t26;
		    zk[i__] = t5 * -.9305257363491 * t15 - t19 * 
			    .9305257363491 * t29;
/* Computing 2nd power */
		    d__1 = t12;
		    t36 = d__1 * d__1;
		    t37 = 1 / t36;
		    vrhoa[i__] = t4 * -1.2407009817988 * t15 + 
			    .008963286558970112 / t4 / t6 * t37 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t26;
		    t45 = d__1 * d__1;
		    t46 = 1 / t45;
		    vrhob[i__] = t18 * -1.2407009817988 * t29 + 
			    .008963286558970112 / t18 / t20 * t46 * sigmabb;
		    vsigmaaa[i__] = -.003361232459613792 / t5 * t37;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = -.003361232459613792 / t19 * t46;
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
/* Computing 2nd power */
		    d__1 = rhob;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    t10 = sigmabb * .00449276922095889 / t5 / t4 + 1.;
		    t13 = 1.804 - .804 / t10;
		    zk[i__] = t3 * -.9305257363491 * t13;
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    vrhob[i__] = t2 * -1.2407009817988 * t13 + 
			    .008963286558970112 / t2 / t4 * t21 * sigmabb;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = -.003361232459613792 / t3 * t21;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t4;
		    t37 = d__1 * d__1;
		    t41 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t43 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t5 * t13 - 
			    .008963286558970112 / t2 / t4 / rhob * t21 * 
			    sigmabb + 2.147732158441357e-4 / t37 / t4 * t41 * 
			    t43;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    v2sigmabb2[i__] = 3.020248347808158e-5 / t37 * t41;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = pow_dd(&rhoa, &c_b2);
		    t3 = t2 * rhoa;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t4 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t2;
		    t5 = d__1 * d__1;
		    t10 = sigmaaa * .00449276922095889 / t5 / t4 + 1.;
		    t13 = 1.804 - .804 / t10;
		    zk[i__] = t3 * -.9305257363491 * t13;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    vrhoa[i__] = t2 * -1.2407009817988 * t13 + 
			    .008963286558970112 / t2 / t4 * t21 * sigmaaa;
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = -.003361232459613792 / t3 * t21;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t4;
		    t37 = d__1 * d__1;
		    t41 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t43 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t5 * t13 - 
			    .008963286558970112 / t2 / t4 / rhoa * t21 * 
			    sigmaaa + 2.147732158441357e-4 / t37 / t4 * t41 * 
			    t43;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    v2sigmaaa2[i__] = 3.020248347808158e-5 / t37 * t41;
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
/* Computing 2nd power */
		    d__1 = rhoa;
		    t6 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t7 = d__1 * d__1;
		    t12 = sigmaaa * .00449276922095889 / t7 / t6 + 1.;
		    t15 = 1.804 - .804 / t12;
		    t18 = pow_dd(&rhob, &c_b2);
		    t19 = t18 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t21 = d__1 * d__1;
		    t26 = sigmabb * .00449276922095889 / t21 / t20 + 1.;
		    t29 = 1.804 - .804 / t26;
		    zk[i__] = t5 * -.9305257363491 * t15 - t19 * 
			    .9305257363491 * t29;
/* Computing 2nd power */
		    d__1 = t12;
		    t36 = d__1 * d__1;
		    t37 = 1 / t36;
		    t38 = 1 / t4 / t6 * t37;
		    vrhoa[i__] = t4 * -1.2407009817988 * t15 + t38 * 
			    .008963286558970112 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t26;
		    t45 = d__1 * d__1;
		    t46 = 1 / t45;
		    t47 = 1 / t18 / t20 * t46;
		    vrhob[i__] = t18 * -1.2407009817988 * t29 + t47 * 
			    .008963286558970112 * sigmabb;
		    vsigmaaa[i__] = -.003361232459613792 / t5 * t37;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = -.003361232459613792 / t19 * t46;
/* Computing 2nd power */
		    d__1 = t6;
		    t65 = d__1 * d__1;
		    t69 = 1 / t36 / t12;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t71 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t7 * t15 - 
			    .008963286558970112 / t4 / t6 / rhoa * t37 * 
			    sigmaaa + 2.147732158441357e-4 / t65 / t6 * t69 * 
			    t71;
/* Computing 2nd power */
		    d__1 = t20;
		    t83 = d__1 * d__1;
		    t87 = 1 / t45 / t26;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t89 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t21 * t29 - 
			    .008963286558970112 / t18 / t20 / rhob * t46 * 
			    sigmabb + 2.147732158441357e-4 / t83 / t20 * t87 *
			     t89;
		    v2rhoab[i__] = 0.;
		    v2rhoasigmaaa[i__] = t38 * .004481643279485056 - 
			    8.053995594155087e-5 / t65 / rhoa * t69 * sigmaaa;
		    v2rhoasigmaab[i__] = 0.;
		    v2rhoasigmabb[i__] = 0.;
		    v2rhobsigmaaa[i__] = 0.;
		    v2rhobsigmaab[i__] = 0.;
		    v2rhobsigmabb[i__] = t47 * .004481643279485056 - 
			    8.053995594155087e-5 / t83 / rhob * t87 * sigmabb;
		    v2sigmaaa2[i__] = 3.020248347808158e-5 / t65 * t69;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    v2sigmabb2[i__] = 3.020248347808158e-5 / t83 * t87;
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
} /* uks_x_pbe__ */

/* Subroutine */ int rks_x_pbe__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal t2, t3, t4, t5, t10, t20, t21, t13, t22, t41, t43, t37, 
	    rho, sigma;


/*     J.P. Perdew, K. Burke, and M. Ernzerhof */
/*     Generalized gradient approximation made simple */
/*     Phys. Rev. Lett. 77 (1996) 3865-3868 */


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
/* Computing 2nd power */
		d__1 = rho;
		t4 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t5 = d__1 * d__1;
		zk[i__] = t2 * -.7385587663820224 * rho * (1.804 - .804 / (
			sigma * .007131826587600489 / t5 / t4 + 1.));
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
/* Computing 2nd power */
		d__1 = rho;
		t4 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t5 = d__1 * d__1;
		t10 = sigma * .007131826587600489 / t5 / t4 + 1.;
		t13 = 1.804 - .804 / t10;
		zk[i__] = t3 * -.7385587663820224 * t13;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		vrhoa[i__] = t2 * -.9847450218426965 * t13 + 
			.01129303341188623 / t2 / t4 * t21 * sigma;
		vsigmaaa[i__] = -.01693955011782934 / t3 * t21;
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
/* Computing 2nd power */
		d__1 = rho;
		t4 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t5 = d__1 * d__1;
		t10 = sigma * .007131826587600489 / t5 / t4 + 1.;
		t13 = 1.804 - .804 / t10;
		zk[i__] = t3 * -.7385587663820224 * t13;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		t22 = 1 / t2 / t4 * t21;
		vrhoa[i__] = t2 * -.9847450218426965 * t13 + t22 * 
			.01129303341188623 * sigma;
		vsigmaaa[i__] = -.01693955011782934 / t3 * t21;
/* Computing 2nd power */
		d__1 = t4;
		t37 = d__1 * d__1;
		t41 = 1 / t20 / t10;
/* Computing 2nd power */
		d__1 = sigma;
		t43 = d__1 * d__1;
		v2rhoa2[i__] = -.6564966812284644 / t5 * t13 - 
			.02258606682377246 / t2 / t4 / rho * t21 * sigma + 
			8.590928633765426e-4 / t37 / t4 * t41 * t43;
		v2rhoasigmaaa[i__] = t22 * .02258606682377246 - 
			6.44319647532407e-4 / t37 / rho * t41 * sigma;
		v2sigmaaa2[i__] = 9.664794712986105e-4 / t37 * t41;
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
} /* rks_x_pbe__ */


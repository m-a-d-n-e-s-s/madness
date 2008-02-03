/* xc_x_lda.f -- translated by f2c (version 20050501).
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

/* :X_LDAsubrstart */
/*    Generated: Wed Jan 29 09:08:46 GMT 2003 */
/* Subroutine */ int uks_x_lda__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal t1, t4, t5, t9, t12, rho, rhoa, rhob;


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
		    t1 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhob;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = pow_dd(&rhoa, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = pow_dd(&rhoa, &c_b2);
		    t4 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa - t4 * 
			    .9305257363491 * rhob;
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
		    t1 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhob;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = t1 * -1.2407009817988;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = pow_dd(&rhoa, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa;
		    vrhoa[i__] = t1 * -1.2407009817988;
		    vrhob[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = pow_dd(&rhoa, &c_b2);
		    t4 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa - t4 * 
			    .9305257363491 * rhob;
		    vrhoa[i__] = t1 * -1.2407009817988;
		    vrhob[i__] = t4 * -1.2407009817988;
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
		    t1 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhob;
		    vrhoa[i__] = 0.;
		    vrhob[i__] = t1 * -1.2407009817988;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t1;
		    t5 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t5;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = pow_dd(&rhoa, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa;
		    vrhoa[i__] = t1 * -1.2407009817988;
		    vrhob[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t1;
		    t5 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t5;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = pow_dd(&rhoa, &c_b2);
		    t4 = pow_dd(&rhob, &c_b2);
		    zk[i__] = t1 * -.9305257363491 * rhoa - t4 * 
			    .9305257363491 * rhob;
		    vrhoa[i__] = t1 * -1.2407009817988;
		    vrhob[i__] = t4 * -1.2407009817988;
/* Computing 2nd power */
		    d__1 = t1;
		    t9 = d__1 * d__1;
		    v2rhoa2[i__] = -.4135669939329333 / t9;
/* Computing 2nd power */
		    d__1 = t4;
		    t12 = d__1 * d__1;
		    v2rhob2[i__] = -.4135669939329333 / t12;
		    v2rhoab[i__] = 0.;
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
} /* uks_x_lda__ */

/* Subroutine */ int rks_x_lda__(integer *ideriv, integer *npt, doublereal *
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


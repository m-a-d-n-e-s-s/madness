/* xc_c_pw92.f -- translated by f2c (version 20050501).
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

/* :C_PW92subrstart */
/*    Generated: Mon Oct 18 16:28:06 BST 2004 */
/* Subroutine */ int uks_c_pw92__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal t1, t2, t3, t4, t5, t6, t9, t11, t30, t40, t13, t33, 
	    t34, t17, t35, t19, t36, t39, t42, t43, t44, t45, t46, t48, t64, 
	    t16, t22, t26, t28, t29, t32, t21, t31, t47, t50, t53, t55, t60, 
	    t63, t67, t68, t70, t73, t74, t75, t76, t77, t80, t81, t84, t88, 
	    t91, t96, t97, t101, t102, t121, t130, t131, t115, t116, t117, 
	    t125, t134, t135, t142, t163, t159, t166, t175, t23, t27, t38, 
	    t41, t57, t65, t72, t99, t51, t78, t79, t82, t83, t87, t90, t93, 
	    t94, t98, t100, t103, t104, t109, t110, t111, t112, t113, t114, 
	    t118, t122, t127, t128, t129, t133, t136, t138, t139, t143, t144, 
	    t149, t150, t156, t157, t158, t160, t161, t162, t165, t169, t172, 
	    t177, t178, rho, t179, t181, t182, t183, t184, t185, t189, t190, 
	    t191, t192, t195, t198, t199, t200, t201, t204, t206, t208, t209, 
	    t212, t213, t214, t220, t226, t228, t229, t230, t231, t234, t236, 
	    t238, t239, t240, t242, t245, t246, t249, t250, t255, t256, t259, 
	    t261, t264, t268, t270, t274, t276, t280, t281, t282, t283, t284, 
	    t300, t316, t319, t324, t326, t330, t335, t338, t345, t347, t348, 
	    t350, t351, t353, t354, t357, t361, t366, t369, t374, t375, t376, 
	    t380, t385, t388, t391, t404, t406, t408, t415, t418, t421, t427, 
	    t439, t441, t443, t450, t453, t456, t460, t463, rhoa, rhob, t465, 
	    t484, t502, t513;


/*     J.P. Perdew, Y. Wang */
/*     Accurate and simple analytic representation of */
/*     the electron-gas correlation energy */
/*     Phys. Rev. B 45 (1992) 13244-13249 */


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
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t17 = log(32.1646831778707 / (t6 * 11.12037486309468 + t2 
			    * 3.844746237447211 + t9 * 1.644733775567609 + 
			    t11 * .2405871291288192) + 1.);
		    zk[i__] = rhob * -.03109 * (t2 * .1274696188700087 + 1.) *
			     t17;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t17 = log(32.1646831778707 / (t6 * 11.12037486309468 + t2 
			    * 3.844746237447211 + t9 * 1.644733775567609 + 
			    t11 * .2405871291288192) + 1.);
		    zk[i__] = rhoa * -.03109 * (t2 * .1274696188700087 + 1.) *
			     t17;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = rhoa + rhob;
		    t2 = 1 / t1;
		    t3 = pow_dd(&t2, &c_b2);
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t17 = log(16.0818243221511 / (t6 * 5.98255043577108 + t3 *
			     2.225569421150687 + t9 * .8004286349993634 + t11 
			    * .1897004325747559) + 1.);
		    t19 = (t3 * .1325688999052018 + 1.) * .062182 * t17;
		    t30 = log(29.60857464321668 / (t6 * 8.157414703487641 + 
			    t3 * 2.247591863577616 + t9 * .4300972471276643 + 
			    t11 * .1911512595127338) + 1.);
		    t33 = rhoa - rhob * 1.;
		    t34 = t33 * t2;
		    t35 = t34 + 1.;
		    t36 = pow_dd(&t35, &c_b2);
		    t39 = 1. - t34 * 1.;
		    t40 = pow_dd(&t39, &c_b2);
		    t42 = t36 * t35 + t40 * t39 - 2.;
/* Computing 2nd power */
		    d__1 = t33;
		    t43 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t43;
		    t44 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t1;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
		    t48 = t44 / t46;
		    t64 = log(32.1646831778707 / (t6 * 11.12037486309468 + t3 
			    * 3.844746237447211 + t9 * 1.644733775567609 + 
			    t11 * .2405871291288192) + 1.);
		    zk[i__] = t1 * (-t19 + (t3 * .06901399211255825 + 1.) * 
			    .03799574853701528 * t30 * t42 * (1. - t48 * 1.) 
			    + ((t3 * .1274696188700087 + 1.) * -.03109 * t64 
			    + t19) * 1.923661050931536 * t42 * t48);
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
		    t4 = t2 * .1274696188700087 + 1.;
		    t5 = rhob * t4;
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t2 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.1646831778707 / t13 + 1.;
		    t17 = log(t16);
		    zk[i__] = t5 * -.03109 * t17;
		    vrhoa[i__] = 0.;
		    t22 = 1 / t11;
/* Computing 2nd power */
		    d__1 = t13;
		    t26 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t28 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t28;
		    t29 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t32 = d__1 * d__1;
		    t33 = 1 / t32;
		    vrhob[i__] = t4 * -.03109 * t17 + t1 * 
			    .001321010150222857 * t22 * t17 + t5 * 1. / t26 * 
			    (-1.853395810515781 / t29 / t6 * t33 - t22 * 
			    1.28158207914907 * t33 - .8223668877838045 / t9 * 
			    t33 - .1603914194192128 / t2 * t33) / t16;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t4 = t2 * .1274696188700087 + 1.;
		    t5 = rhoa * t4;
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t2 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.1646831778707 / t13 + 1.;
		    t17 = log(t16);
		    zk[i__] = t5 * -.03109 * t17;
		    t22 = 1 / t11;
/* Computing 2nd power */
		    d__1 = t13;
		    t26 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t28 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t28;
		    t29 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t32 = d__1 * d__1;
		    t33 = 1 / t32;
		    vrhoa[i__] = t4 * -.03109 * t17 + t1 * 
			    .001321010150222857 * t22 * t17 + t5 * 1. / t26 * 
			    (-1.853395810515781 / t29 / t6 * t33 - t22 * 
			    1.28158207914907 * t33 - .8223668877838045 / t9 * 
			    t33 - .1603914194192128 / t2 * t33) / t16;
		    vrhob[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = rhoa + rhob;
		    t2 = 1 / t1;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1325688999052018 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 5.98255043577108 + t3 * 2.225569421150687 + t9 
			    * .8004286349993634 + t11 * .1897004325747559;
		    t16 = 16.0818243221511 / t13 + 1.;
		    t17 = log(t16);
		    t19 = t5 * .062182 * t17;
		    t21 = t3 * .06901399211255825 + 1.;
		    t26 = t6 * 8.157414703487641 + t3 * 2.247591863577616 + 
			    t9 * .4300972471276643 + t11 * .1911512595127338;
		    t29 = 29.60857464321668 / t26 + 1.;
		    t30 = log(t29);
		    t31 = t21 * t30;
		    t33 = rhoa - rhob * 1.;
		    t34 = t33 * t2;
		    t35 = t34 + 1.;
		    t36 = pow_dd(&t35, &c_b2);
		    t39 = 1. - t34 * 1.;
		    t40 = pow_dd(&t39, &c_b2);
		    t42 = t36 * t35 + t40 * t39 - 2.;
/* Computing 2nd power */
		    d__1 = t33;
		    t43 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t43;
		    t44 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t1;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
		    t47 = 1 / t46;
		    t48 = t44 * t47;
		    t50 = 1. - t48 * 1.;
		    t53 = t31 * .03799574853701528 * t42 * t50;
		    t55 = t3 * .1274696188700087 + 1.;
		    t60 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t63 = 32.1646831778707 / t60 + 1.;
		    t64 = log(t63);
		    t67 = t55 * -.03109 * t64 + t19;
		    t68 = t67 * t42;
		    t70 = t68 * 1.923661050931536 * t48;
		    zk[i__] = t1 * (-t19 + t53 + t70);
		    t73 = 1 / t45;
		    t74 = 1 / t11 * t73;
		    t75 = t74 * t17;
		    t76 = t75 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t13;
		    t77 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t80 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t80;
		    t81 = d__1 * d__1;
		    t84 = 1 / t81 / t6 * t73;
		    t88 = 1 / t9 * t73;
		    t91 = 1 / t3 * t73;
		    t96 = t5 / t77 * (t84 * -.99709173929518 - t74 * 
			    .7418564737168958 - t88 * .4002143174996817 - t91 
			    * .1264669550498372) / t16;
		    t97 = t96 * 1.;
		    t101 = t74 * 8.740794299481065e-4 * t30 * t42 * t50;
/* Computing 2nd power */
		    d__1 = t26;
		    t102 = d__1 * d__1;
		    t115 = t21 * 1.124999956683108 / t102 * (t84 * 
			    -1.35956911724794 - t74 * .7491972878592054 - t88 
			    * .2150486235638321 - t91 * .1274341730084892) / 
			    t29 * t42 * t50;
		    t116 = t33 * t73;
		    t117 = t116 * 1.;
		    t121 = t2 * 1.;
		    t125 = t36 * 1.333333333333333 * (t2 - t117) + t40 * 
			    1.333333333333333 * (-t121 + t116);
		    t130 = t43 * t33 * t47;
		    t131 = t130 * 4.;
		    t134 = t44 / t46 / t1;
		    t135 = t134 * 4.;
/* Computing 2nd power */
		    d__1 = t60;
		    t142 = d__1 * d__1;
		    t159 = (t74 * .001321010150222857 * t64 + t55 * 1. / t142 
			    * (t84 * -1.853395810515781 - t74 * 
			    1.28158207914907 - t88 * .8223668877838045 - t91 *
			     .1603914194192128) / t63 - t75 * 
			    .002747799777968419 - t96 * 1.) * 
			    1.923661050931536 * t42 * t48;
		    t163 = t68 * t130;
		    t166 = t68 * 7.694644203726145 * t134;
		    vrhoa[i__] = -t19 + t53 + t70 + t1 * (t76 + t97 - t101 - 
			    t115 + t31 * .03799574853701528 * t125 * t50 + 
			    t31 * .03799574853701528 * t42 * (-t131 + t135) + 
			    t159 + t67 * 1.923661050931536 * t125 * t48 + 
			    t163 * 7.694644203726145 - t166);
		    t175 = t36 * 1.333333333333333 * (-t121 - t117) + t40 * 
			    1.333333333333333 * (t2 + t116);
		    vrhob[i__] = -t19 + t53 + t70 + t1 * (t76 + t97 - t101 - 
			    t115 + t31 * .03799574853701528 * t175 * t50 + 
			    t31 * .03799574853701528 * t42 * (t131 + t135) + 
			    t159 + t67 * 1.923661050931536 * t175 * t48 - 
			    t163 * 7.694644203726145 - t166);
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
		    t4 = t2 * .1274696188700087 + 1.;
		    t5 = rhob * t4;
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t2 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.1646831778707 / t13 + 1.;
		    t17 = log(t16);
		    zk[i__] = t5 * -.03109 * t17;
		    vrhoa[i__] = 0.;
		    t22 = 1 / t11;
		    t23 = t1 * t22;
/* Computing 2nd power */
		    d__1 = t13;
		    t26 = d__1 * d__1;
		    t27 = 1 / t26;
/* Computing 2nd power */
		    d__1 = t6;
		    t28 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t28;
		    t29 = d__1 * d__1;
		    t30 = t29 * t6;
		    t31 = 1 / t30;
/* Computing 2nd power */
		    d__1 = rhob;
		    t32 = d__1 * d__1;
		    t33 = 1 / t32;
		    t38 = 1 / t9;
		    t41 = 1 / t2;
		    t44 = t31 * -1.853395810515781 * t33 - t22 * 
			    1.28158207914907 * t33 - t38 * .8223668877838045 *
			     t33 - t41 * .1603914194192128 * t33;
		    t46 = 1 / t16;
		    t47 = t27 * t44 * t46;
		    vrhob[i__] = t4 * -.03109 * t17 + t23 * 
			    .001321010150222857 * t17 + t5 * 1. * t47;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t55 = 1 / t32 / rhob;
		    t57 = 1 / t11 / t1;
/* Computing 2nd power */
		    d__1 = t44;
		    t65 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t72 = d__1 * d__1;
		    t73 = 1 / t72;
/* Computing 2nd power */
		    d__1 = t26;
		    t99 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t16;
		    t102 = d__1 * d__1;
		    v2rhob2[i__] = t4 * 2. * t27 * t44 * t46 + t55 * 
			    8.806734334819047e-4 * t57 * t17 - t23 * 
			    .08497974591333914 * t47 - t5 * 2. / t26 / t13 * 
			    t65 * t46 + t5 * 1. * t27 * (-1.544496508763151 / 
			    t30 / t1 * t73 + t31 * 3.706791621031562 * t55 - 
			    t57 * .854388052766047 * t73 + t22 * 
			    2.563164158298141 * t55 - .4111834438919023 / t9 /
			     t1 * t73 + t38 * 1.644733775567609 * t55 - 
			    .05346380647307093 / t2 / t1 * t73 + t41 * 
			    .3207828388384256 * t55) * t46 + t5 * 
			    32.1646831778707 / t99 * t65 / t102;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
		    t1 = 1 / rhoa;
		    t2 = pow_dd(&t1, &c_b2);
		    t4 = t2 * .1274696188700087 + 1.;
		    t5 = rhoa * t4;
		    t6 = pow_dd(&t1, &c_b3);
		    t9 = sqrt(t1);
/* Computing 2nd power */
		    d__1 = t2;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t2 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.1646831778707 / t13 + 1.;
		    t17 = log(t16);
		    zk[i__] = t5 * -.03109 * t17;
		    t22 = 1 / t11;
		    t23 = t1 * t22;
/* Computing 2nd power */
		    d__1 = t13;
		    t26 = d__1 * d__1;
		    t27 = 1 / t26;
/* Computing 2nd power */
		    d__1 = t6;
		    t28 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t28;
		    t29 = d__1 * d__1;
		    t30 = t29 * t6;
		    t31 = 1 / t30;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t32 = d__1 * d__1;
		    t33 = 1 / t32;
		    t38 = 1 / t9;
		    t41 = 1 / t2;
		    t44 = t31 * -1.853395810515781 * t33 - t22 * 
			    1.28158207914907 * t33 - t38 * .8223668877838045 *
			     t33 - t41 * .1603914194192128 * t33;
		    t46 = 1 / t16;
		    t47 = t27 * t44 * t46;
		    vrhoa[i__] = t4 * -.03109 * t17 + t23 * 
			    .001321010150222857 * t17 + t5 * 1. * t47;
		    vrhob[i__] = 0.;
		    t55 = 1 / t32 / rhoa;
		    t57 = 1 / t11 / t1;
/* Computing 2nd power */
		    d__1 = t44;
		    t65 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t72 = d__1 * d__1;
		    t73 = 1 / t72;
/* Computing 2nd power */
		    d__1 = t26;
		    t99 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t16;
		    t102 = d__1 * d__1;
		    v2rhoa2[i__] = t4 * 2. * t27 * t44 * t46 + t55 * 
			    8.806734334819047e-4 * t57 * t17 - t23 * 
			    .08497974591333914 * t47 - t5 * 2. / t26 / t13 * 
			    t65 * t46 + t5 * 1. * t27 * (-1.544496508763151 / 
			    t30 / t1 * t73 + t31 * 3.706791621031562 * t55 - 
			    t57 * .854388052766047 * t73 + t22 * 
			    2.563164158298141 * t55 - .4111834438919023 / t9 /
			     t1 * t73 + t38 * 1.644733775567609 * t55 - 
			    .05346380647307093 / t2 / t1 * t73 + t41 * 
			    .3207828388384256 * t55) * t46 + t5 * 
			    32.1646831778707 / t99 * t65 / t102;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		} else {
/* (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol)) */
		    t1 = rhoa + rhob;
		    t2 = 1 / t1;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1325688999052018 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 5.98255043577108 + t3 * 2.225569421150687 + t9 
			    * .8004286349993634 + t11 * .1897004325747559;
		    t16 = 16.0818243221511 / t13 + 1.;
		    t17 = log(t16);
		    t19 = t5 * .062182 * t17;
		    t21 = t3 * .06901399211255825 + 1.;
		    t26 = t6 * 8.157414703487641 + t3 * 2.247591863577616 + 
			    t9 * .4300972471276643 + t11 * .1911512595127338;
		    t29 = 29.60857464321668 / t26 + 1.;
		    t30 = log(t29);
		    t31 = t21 * t30;
		    t33 = rhoa - rhob * 1.;
		    t34 = t33 * t2;
		    t35 = t34 + 1.;
		    t36 = pow_dd(&t35, &c_b2);
		    t39 = 1. - t34 * 1.;
		    t40 = pow_dd(&t39, &c_b2);
		    t42 = t36 * t35 + t40 * t39 - 2.;
/* Computing 2nd power */
		    d__1 = t33;
		    t43 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t43;
		    t44 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t1;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
		    t47 = 1 / t46;
		    t48 = t44 * t47;
		    t50 = 1. - t48 * 1.;
		    t51 = t42 * t50;
		    t53 = t31 * .03799574853701528 * t51;
		    t55 = t3 * .1274696188700087 + 1.;
		    t60 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t63 = 32.1646831778707 / t60 + 1.;
		    t64 = log(t63);
		    t67 = t55 * -.03109 * t64 + t19;
		    t68 = t67 * t42;
		    t70 = t68 * 1.923661050931536 * t48;
		    zk[i__] = t1 * (-t19 + t53 + t70);
		    t72 = 1 / t11;
		    t73 = 1 / t45;
		    t74 = t72 * t73;
		    t75 = t74 * t17;
		    t76 = t75 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t13;
		    t77 = d__1 * d__1;
		    t78 = 1 / t77;
		    t79 = t5 * t78;
/* Computing 2nd power */
		    d__1 = t6;
		    t80 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t80;
		    t81 = d__1 * d__1;
		    t82 = t81 * t6;
		    t83 = 1 / t82;
		    t84 = t83 * t73;
		    t87 = 1 / t9;
		    t88 = t87 * t73;
		    t90 = 1 / t3;
		    t91 = t90 * t73;
		    t93 = t84 * -.99709173929518 - t74 * .7418564737168958 - 
			    t88 * .4002143174996817 - t91 * .1264669550498372;
		    t94 = 1 / t16;
		    t96 = t79 * t93 * t94;
		    t97 = t96 * 1.;
		    t98 = t30 * t42;
		    t99 = t98 * t50;
		    t100 = t74 * t99;
		    t101 = t100 * 8.740794299481065e-4;
/* Computing 2nd power */
		    d__1 = t26;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
		    t104 = t21 * t103;
		    t109 = t84 * -1.35956911724794 - t74 * .7491972878592054 
			    - t88 * .2150486235638321 - t91 * 
			    .1274341730084892;
		    t110 = t104 * t109;
		    t111 = 1 / t29;
		    t112 = t111 * t42;
		    t113 = t112 * t50;
		    t114 = t110 * t113;
		    t115 = t114 * 1.124999956683108;
		    t116 = t33 * t73;
		    t117 = t116 * 1.;
		    t118 = t2 - t117;
		    t121 = t2 * 1.;
		    t122 = -t121 + t116;
		    t125 = t36 * 1.333333333333333 * t118 + t40 * 
			    1.333333333333333 * t122;
		    t127 = t31 * t125 * t50;
		    t128 = t127 * .03799574853701528;
		    t129 = t43 * t33;
		    t130 = t129 * t47;
		    t131 = t130 * 4.;
		    t133 = 1 / t46 / t1;
		    t134 = t44 * t133;
		    t135 = t134 * 4.;
		    t136 = -t131 + t135;
		    t138 = t31 * t42 * t136;
		    t139 = t138 * .03799574853701528;
/* Computing 2nd power */
		    d__1 = t60;
		    t142 = d__1 * d__1;
		    t143 = 1 / t142;
		    t144 = t55 * t143;
		    t149 = t84 * -1.853395810515781 - t74 * 1.28158207914907 
			    - t88 * .8223668877838045 - t91 * 
			    .1603914194192128;
		    t150 = 1 / t63;
		    t156 = t74 * .001321010150222857 * t64 + t144 * 1. * t149 
			    * t150 - t75 * .002747799777968419 - t96 * 1.;
		    t157 = t156 * t42;
		    t158 = t157 * t48;
		    t159 = t158 * 1.923661050931536;
		    t160 = t67 * t125;
		    t161 = t160 * t48;
		    t162 = t161 * 1.923661050931536;
		    t163 = t68 * t130;
		    t165 = t68 * t134;
		    t166 = t165 * 7.694644203726145;
		    vrhoa[i__] = -t19 + t53 + t70 + t1 * (t76 + t97 - t101 - 
			    t115 + t128 + t139 + t159 + t162 + t163 * 
			    7.694644203726145 - t166);
		    t169 = -t121 - t117;
		    t172 = t2 + t116;
		    t175 = t36 * 1.333333333333333 * t169 + t40 * 
			    1.333333333333333 * t172;
		    t177 = t31 * t175 * t50;
		    t178 = t177 * .03799574853701528;
		    t179 = t131 + t135;
		    t181 = t31 * t42 * t179;
		    t182 = t181 * .03799574853701528;
		    t183 = t67 * t175;
		    t184 = t183 * t48;
		    t185 = t184 * 1.923661050931536;
		    vrhob[i__] = -t19 + t53 + t70 + t1 * (t76 + t97 - t101 - 
			    t115 + t178 + t182 + t159 + t185 - t163 * 
			    7.694644203726145 - t166);
		    t189 = t75 * .005495599555936838;
		    t190 = t96 * 2.;
		    t191 = t100 * .001748158859896213;
		    t192 = t114 * 2.249999913366216;
		    t195 = t158 * 3.847322101863073;
		    t198 = t165 * 15.38928840745229;
/* Computing 2nd power */
		    d__1 = t36;
		    t199 = d__1 * d__1;
		    t200 = 1 / t199;
/* Computing 2nd power */
		    d__1 = t118;
		    t201 = d__1 * d__1;
		    t204 = t73 * 2.;
		    t206 = 1 / t45 / t1;
		    t208 = t33 * 2. * t206;
		    t209 = -t204 + t208;
/* Computing 2nd power */
		    d__1 = t40;
		    t212 = d__1 * d__1;
		    t213 = 1 / t212;
/* Computing 2nd power */
		    d__1 = t122;
		    t214 = d__1 * d__1;
		    t220 = t200 * .4444444444444444 * t201 + t36 * 
			    1.333333333333333 * t209 + t213 * 
			    .4444444444444444 * t214 - t40 * 
			    1.333333333333333 * t209;
		    t226 = 1 / t11 / t2 * t47;
		    t228 = t226 * 5.827196199654043e-4 * t99;
		    t229 = t43 * t47;
		    t230 = t68 * t229;
		    t231 = t230 * 23.08393261117844;
		    t234 = t44 / t46 / t45;
		    t236 = t68 * 38.47322101863073 * t234;
		    t238 = t157 * 15.38928840745229 * t134;
		    t239 = t129 * t133;
		    t240 = t68 * t239;
/* Computing 2nd power */
		    d__1 = t77;
		    t242 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t93;
		    t245 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t16;
		    t246 = d__1 * d__1;
		    t249 = t5 / t242 * t245 / t246;
		    t250 = t249 * 16.0818243221511;
		    t255 = t5 / t77 / t13 * t245 * t94;
		    t256 = t255 * 2.;
		    t259 = 1 / t82 / t2 * t47;
		    t261 = t83 * t206;
		    t264 = t72 * t206;
		    t268 = 1 / t9 / t2 * t47;
		    t270 = t87 * t206;
		    t274 = 1 / t3 / t2 * t47;
		    t276 = t90 * t206;
		    t280 = t79 * (t259 * -.8309097827459833 + t261 * 
			    1.99418347859036 - t226 * .4945709824779306 + 
			    t264 * 1.483712947433792 - t268 * 
			    .2001071587498409 + t270 * .8004286349993634 - 
			    t274 * .04215565168327908 + t276 * 
			    .2529339100996745) * t94;
		    t281 = t280 * 1.;
		    t282 = t229 * 12.;
		    t283 = t239 * 32.;
		    t284 = t234 * 20.;
/* Computing 2nd power */
		    d__1 = t149;
		    t300 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t142;
		    t316 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t63;
		    t319 = d__1 * d__1;
		    t324 = t226 * t17;
		    t326 = t264 * t17;
		    t330 = t74 * t78 * t93 * t94;
		    t335 = t226 * 8.806734334819047e-4 * t64 - t264 * 
			    .002642020300445714 * t64 - t74 * 
			    .08497974591333914 * t143 * t149 * t150 - t55 * 
			    2. / t142 / t60 * t300 * t150 + t144 * 1. * (t259 
			    * -1.544496508763151 + t261 * 3.706791621031562 - 
			    t226 * .854388052766047 + t264 * 
			    2.563164158298141 - t268 * .4111834438919023 + 
			    t270 * 1.644733775567609 - t274 * 
			    .05346380647307093 + t276 * .3207828388384256) * 
			    t150 + t55 * 32.1646831778707 / t316 * t300 / 
			    t319 - t324 * .001831866518645613 + t326 * 
			    .005495599555936838 + t330 * .08837926660346786 + 
			    t255 * 2. - t280 * 1. - t249 * 16.0818243221511;
		    t338 = t335 * 1.923661050931536 * t42 * t48;
		    t345 = t157 * t130;
		    t347 = t31 * .03799574853701528 * t220 * t50 - t228 + 
			    t231 + t236 - t238 - t240 * 61.55715362980916 + 
			    t250 - t256 + t281 + t31 * .03799574853701528 * 
			    t42 * (-t282 + t283 - t284) + t338 + t31 * 
			    .07599149707403056 * t125 * t136 + t67 * 
			    1.923661050931536 * t220 * t48 + t345 * 
			    15.38928840745229;
		    t348 = t160 * t130;
		    t350 = t324 * .001831866518645613;
		    t351 = t326 * .005495599555936838;
		    t353 = t264 * .001748158859896213 * t99;
		    t354 = t160 * t134;
		    t357 = t74 * t98 * t136;
		    t361 = t74 * t30 * t125 * t50;
/* Computing 2nd power */
		    d__1 = t109;
		    t366 = d__1 * d__1;
		    t369 = t21 * 2.249999913366216 / t102 / t26 * t366 * t113;
		    t374 = t74 * .05176049209143758 * t103 * t109 * t111 * 
			    t51;
		    t375 = t330 * .08837926660346786;
/* Computing 2nd power */
		    d__1 = t102;
		    t376 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t29;
		    t380 = d__1 * d__1;
		    t385 = t21 * 33.30964519106732 / t376 * t366 / t380 * t42 
			    * t50;
		    t388 = t110 * t111 * t125 * t50;
		    t391 = t110 * t112 * t136;
		    t404 = t104 * 1.124999956683108 * (t259 * 
			    -1.132974264373283 + t261 * 2.71913823449588 - 
			    t226 * .4994648585728036 + t264 * 
			    1.498394575718411 - t268 * .1075243117819161 + 
			    t270 * .4300972471276643 - t274 * 
			    .04247805766949639 + t276 * .2548683460169784) * 
			    t113;
		    t406 = t156 * t125 * t48;
		    t408 = t348 * 15.38928840745229 + t350 - t351 + t353 - 
			    t354 * 15.38928840745229 - t357 * 
			    .001748158859896213 - t361 * .001748158859896213 
			    + t369 + t374 - t375 - t385 - t388 * 
			    2.249999913366216 - t391 * 2.249999913366216 - 
			    t404 + t406 * 3.847322101863073;
		    v2rhoa2[i__] = t189 + t190 - t191 - t192 + t127 * 
			    .07599149707403056 + t138 * .07599149707403056 + 
			    t195 + t161 * 3.847322101863073 + t163 * 
			    15.38928840745229 - t198 + t1 * (t347 + t408);
/* Computing 2nd power */
		    d__1 = t169;
		    t415 = d__1 * d__1;
		    t418 = t204 + t208;
/* Computing 2nd power */
		    d__1 = t172;
		    t421 = d__1 * d__1;
		    t427 = t200 * .4444444444444444 * t415 + t36 * 
			    1.333333333333333 * t418 + t213 * 
			    .4444444444444444 * t421 - t40 * 
			    1.333333333333333 * t418;
		    t439 = t156 * t175 * t48;
		    t441 = t183 * t134;
		    t443 = t183 * t130;
		    t450 = t110 * t111 * t175 * t50;
		    t453 = t31 * .03799574853701528 * t427 * t50 + t31 * 
			    .07599149707403056 * t175 * t179 + t31 * 
			    .03799574853701528 * t42 * (-t282 - t283 - t284) 
			    + t439 * 3.847322101863073 - t441 * 
			    15.38928840745229 - t443 * 15.38928840745229 + 
			    t67 * 1.923661050931536 * t427 * t48 - t450 * 
			    2.249999913366216 - t228 + t231 + t236 - t238 + 
			    t240 * 61.55715362980916 + t250;
		    t456 = t74 * t98 * t179;
		    t460 = t74 * t30 * t175 * t50;
		    t463 = t110 * t112 * t179;
		    t465 = -t256 + t281 + t338 - t345 * 15.38928840745229 + 
			    t350 - t351 + t353 + t369 + t374 - t375 - t385 - 
			    t404 - t456 * .001748158859896213 - t460 * 
			    .001748158859896213 - t463 * 2.249999913366216;
		    v2rhob2[i__] = t189 + t190 - t191 - t192 + t177 * 
			    .07599149707403056 + t181 * .07599149707403056 + 
			    t195 + t184 * 3.847322101863073 - t163 * 
			    15.38928840745229 - t198 + t1 * (t453 + t465);
		    t484 = t200 * .4444444444444444 * t118 * t169 + t36 * 
			    2.666666666666667 * t33 * t206 + t213 * 
			    .4444444444444444 * t122 * t172 - t40 * 
			    2.666666666666667 * t33 * t206;
		    t502 = t439 * 1.923661050931536 - t441 * 
			    7.694644203726145 + t443 * 7.694644203726145 - 
			    t450 * 1.124999956683108 + t31 * 
			    .03799574853701528 * t484 * t50 + t31 * 
			    .03799574853701528 * t125 * t179 + t31 * 
			    .03799574853701528 * t175 * t136 + t31 * 
			    .03799574853701528 * t42 * (t282 - t284) + t67 * 
			    1.923661050931536 * t484 * t48 - t228 - t230 * 
			    23.08393261117844 + t236 - t238 + t250 - t256 + 
			    t281 + t338;
		    t513 = t348 * -7.694644203726145 + t350 - t351 + t353 - 
			    t354 * 7.694644203726145 - t357 * 
			    8.740794299481065e-4 - t361 * 
			    8.740794299481065e-4 + t369 + t374 - t375 - t385 
			    - t388 * 1.124999956683108 - t391 * 
			    1.124999956683108 - t404 + t406 * 
			    1.923661050931536 - t456 * 8.740794299481065e-4 - 
			    t460 * 8.740794299481065e-4 - t463 * 
			    1.124999956683108;
		    v2rhoab[i__] = t189 + t190 - t191 - t192 + t178 + t182 + 
			    t195 + t185 - t198 + t128 + t139 + t162 + t1 * (
			    t502 + t513);
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
} /* uks_c_pw92__ */

/* Subroutine */ int rks_c_pw92__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal t1, t2, t4, t6, t9, t11, t30, t13, t23, t24, t16, t17, 
	    t25, t28, t31, t32, t22, t26, t29, t33, t34, t38, t41, t44, t45, 
	    t47, t55, t56, t57, t59, t61, t62, t64, t75, t77, t82, t83, t86, 
	    t87, t91, t121, t115, rho;


/*     J.P. Perdew, Y. Wang */
/*     Accurate and simple analytic representation of */
/*     the electron-gas correlation energy */
/*     Phys. Rev. B 45 (1992) 13244-13249 */


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
		t6 = pow_dd(&t1, &c_b3);
		t9 = sqrt(t1);
/* Computing 2nd power */
		d__1 = t2;
		t11 = d__1 * d__1;
		t17 = log(16.0818243221511 / (t6 * 5.98255043577108 + t2 * 
			2.225569421150687 + t9 * .8004286349993634 + t11 * 
			.1897004325747559) + 1.);
		zk[i__] = rho * -.062182 * (t2 * .1325688999052018 + 1.) * 
			t17;
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
		t4 = t2 * .1325688999052018 + 1.;
		t6 = pow_dd(&t1, &c_b3);
		t9 = sqrt(t1);
/* Computing 2nd power */
		d__1 = t2;
		t11 = d__1 * d__1;
		t13 = t6 * 5.98255043577108 + t2 * 2.225569421150687 + t9 * 
			.8004286349993634 + t11 * .1897004325747559;
		t16 = 16.0818243221511 / t13 + 1.;
		t17 = log(t16);
		zk[i__] = rho * -.062182 * t4 * t17;
/* Computing 2nd power */
		d__1 = rho;
		t23 = d__1 * d__1;
		t24 = 1 / t23;
		t25 = 1 / t11 * t24;
/* Computing 2nd power */
		d__1 = t13;
		t28 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t6;
		t31 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t31;
		t32 = d__1 * d__1;
		vrhoa[i__] = t4 * -.062182 * t17 + rho * (t25 * 
			.002747799777968419 * t17 + t4 * 1. / t28 * (
			-.99709173929518 / t32 / t6 * t24 - t25 * 
			.7418564737168958 - .4002143174996817 / t9 * t24 - 
			.1264669550498372 / t2 * t24) / t16);
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
		t4 = t2 * .1325688999052018 + 1.;
		t6 = pow_dd(&t1, &c_b3);
		t9 = sqrt(t1);
/* Computing 2nd power */
		d__1 = t2;
		t11 = d__1 * d__1;
		t13 = t6 * 5.98255043577108 + t2 * 2.225569421150687 + t9 * 
			.8004286349993634 + t11 * .1897004325747559;
		t16 = 16.0818243221511 / t13 + 1.;
		t17 = log(t16);
		zk[i__] = rho * -.062182 * t4 * t17;
		t22 = 1 / t11;
/* Computing 2nd power */
		d__1 = rho;
		t23 = d__1 * d__1;
		t24 = 1 / t23;
		t25 = t22 * t24;
		t26 = t25 * t17;
/* Computing 2nd power */
		d__1 = t13;
		t28 = d__1 * d__1;
		t29 = 1 / t28;
		t30 = t4 * t29;
/* Computing 2nd power */
		d__1 = t6;
		t31 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t31;
		t32 = d__1 * d__1;
		t33 = t32 * t6;
		t34 = 1 / t33;
		t38 = 1 / t9;
		t41 = 1 / t2;
		t44 = t34 * -.99709173929518 * t24 - t25 * .7418564737168958 
			- t38 * .4002143174996817 * t24 - t41 * 
			.1264669550498372 * t24;
		t45 = 1 / t16;
		t47 = t30 * t44 * t45;
		vrhoa[i__] = t4 * -.062182 * t17 + rho * (t26 * 
			.002747799777968419 + t47 * 1.);
/* Computing 2nd power */
		d__1 = t23;
		t55 = d__1 * d__1;
		t56 = 1 / t55;
		t57 = 1 / t11 / t1 * t56;
		t59 = t57 * .001831866518645613 * t17;
		t61 = 1 / t23 / rho;
		t62 = t22 * t61;
		t64 = t62 * .005495599555936838 * t17;
		t75 = log(29.60857464321668 / (t6 * 8.157414703487641 + t2 * 
			2.247591863577616 + t9 * .4300972471276643 + t11 * 
			.1911512595127338) + 1.);
		t77 = (t2 * .06901399211255825 + 1.) * t75 * t24;
		t82 = t25 * .08837926660346786 * t29 * t44 * t45;
/* Computing 2nd power */
		d__1 = t28;
		t83 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t44;
		t86 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t16;
		t87 = d__1 * d__1;
		t91 = t4 * 16.0818243221511 / t83 * t86 / t87;
		t115 = t30 * 1. * (-.8309097827459833 / t33 / t1 * t56 + t34 *
			 1.99418347859036 * t61 - t57 * .4945709824779306 + 
			t62 * 1.483712947433792 - .2001071587498409 / t9 / t1 
			* t56 + t38 * .8004286349993634 * t61 - 
			.04215565168327908 / t2 / t1 * t56 + t41 * 
			.2529339100996745 * t61) * t45;
		t121 = t4 * 2. / t28 / t13 * t86 * t45;
		v2rhoa2[i__] = t47 * 4. + t26 * .01099119911187368 + rho * (
			t59 - t64 + t77 * .03377399869956914 - t82 + t91 + 
			t115 - t121) + rho * (t59 - t64 - t77 * 
			.03377399869956914 - t82 + t91 + t115 - t121);
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
} /* rks_c_pw92__ */


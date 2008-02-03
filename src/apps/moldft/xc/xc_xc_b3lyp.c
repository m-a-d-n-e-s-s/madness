/* xc_xc_b3lyp.f -- translated by f2c (version 20050501).
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
static doublereal c_b4 = .16666666666666666;

/* :XC_B3LYPsubrstart */
/*    Generated: Thu Feb 13 08:52:41 GMT 2003 */
/* Subroutine */ int uks_xc_b3lyp__(integer *ideriv, integer *npt, doublereal 
	*rhoa1, doublereal *rhob1, doublereal *sigmaaa1, doublereal *sigmabb1,
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
	    doublereal), atan(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t4, t5, t6, t7, t8, t9, s3, s4, s5, t10,
	     t11, t21, t22, t31, t23, t16, t17, t25, t19, t35, t37, t18, t24, 
	    t32, t33, t36, t38, t39, t42, t43, t45, t47, t51, t52, t55, t56, 
	    t60, t82, t93, t94, t95, t98, t12, t13, t28, t34, t49, t50, t63, 
	    t64, t74, t75, t77, t80, t81, t85, t86, t101, t111, t103, t113, 
	    t120, t130, t107, t117, t109, t125, t128, t138, t139, t140, t143, 
	    t144, t88, t89, t91, t99, t100, t14, t15, t29, t44, t53, t57, t62,
	     t67, t71, t78, t83, t97, t104, t110, t116, t122, t127, t135, 
	    t137, t146, t147, t154, t158, t159, t167, t170, t171, t184, t195, 
	    t200, t203, t204, t210, t211, t212, t214, t215, t217, t221, t225, 
	    t230, t232, t234, t238, t273, t274, rho, t275, t276, t277, t278, 
	    t279, t283, t284, t287, t288, t289, t291, t292, t294, t298, t300, 
	    t302, t303, t309, t320, t325, t326, t329, t336, t337, t373, t376, 
	    t380, t381, t389, t392, t393, t430, t431, t445, t446, t46, t54, 
	    t68, t76, t87, t96, t106, t108, t114, t115, t145, t152, t164, 
	    t175, t177, t178, t188, t190, t198, t199, t201, t272, t58, t69, 
	    t155, t163, t165, t169, t174, t182, t187, t191, t205, t207, t216, 
	    t218, t219, t222, t227, t240, t251, t259, t262, t270, t280, t285, 
	    t290, t297, t299, t304, t306, t307, t311, t314, t317, t318, t319, 
	    t322, rhoa, rhob, t327, t332, t333, t338, t340, t341, t345, t348, 
	    t351, t352, t353, t359, t360, t362, t365, t368, t369, t377, t385, 
	    t387, t391, t396, t397, t403, t407, t411, t415, t420, t426, t427, 
	    t432, t435, t437, t443, t450, t456, t457, t460, t463, t465, t469, 
	    t470, t471, t472, t473, t478, t483, t484, t485, t488, t490, t491, 
	    t494, t497, t500, t501, t502, t503, t504, t506, t510, t513, t514, 
	    t549, t550, t553, t555, t559, t560, t565, t569, t575, t581, t588, 
	    t589, t590, t593, t595, t596, t605, t607, t611, t613, t615, t616, 
	    t618, t620, t625, t628, t631, t633, t638, t639, t642, t645, t647, 
	    sigma, t651, t654, t672, t679, t681, t684, t693, t695, t700, t716,
	     t725, t729, t771, t776, t777, t778, t784, t785, t794, t795, t796,
	     t814, t818, t823, t835, t837, t840, t843, t846, t849, t856, t859,
	     t861, t862, t868, t870, t872, t874, t876, t878, t879, t884, t885,
	     t890, t894, t900, t906, t923, t926, t929, t940, t942, t950, t954,
	     t956, t958, t960, t962, t964, t974, t1003, t1009, t1033, t1040, 
	    t1043, t1046, t1047, t1048, t1049, t1057, t1061, t1063, t1064, 
	    t1066, t1067, t1069, t1070, t1072, t1073, t1078, t1079, t1081, 
	    t1082, t1083, t1084, t1093, t1096, t1099, t1110, t1119, t1120, 
	    t1126, t1155, t1179, sigmaaa, sigmaab, sigmabb;


/*     P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch */
/*     Ab initio calculation of vibrational absorption and circular */
/*     dichroism spectra using density functional force fields */
/*     J. Phys. Chem. 98 (1994) 11623-11627 */


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
		    t16 = 1 / rhob;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t22 = 1 / (t17 * .6203504908994 + t19 * 15.84942278842832 
			    + 101.578);
		    t25 = log(t17 * .6203504908994 * t22);
		    t31 = atan(1.171685277708993 / (t19 * 1.575246635799487 + 
			    20.1231));
/* Computing 2nd power */
		    d__1 = t19 * .7876233178997433 + .743294;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t5 * .003024 * sigmabb /
			     (t8 * .0252 * t9 + 1.) + rhob * .19 * (t25 * 
			    .01554535 + t31 * .6188180297906063 + t37 * 
			    .002667310007273315);
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
		    t16 = 1 / rhoa;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t22 = 1 / (t17 * .6203504908994 + t19 * 15.84942278842832 
			    + 101.578);
		    t25 = log(t17 * .6203504908994 * t22);
		    t31 = atan(1.171685277708993 / (t19 * 1.575246635799487 + 
			    20.1231));
/* Computing 2nd power */
		    d__1 = t19 * .7876233178997433 + .743294;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t5 * .003024 * sigmaaa /
			     (t8 * .0252 * t9 + 1.) + rhoa * .19 * (t25 * 
			    .01554535 + t31 * .6188180297906063 + t37 * 
			    .002667310007273315);
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
		    t32 = pow_dd(&rho, &c_b2);
		    t33 = 1 / t32;
		    t36 = 1 / (t33 * .349 + 1.);
		    t38 = 1 / rho;
		    t39 = rhob * t38;
		    t42 = t33 * .2533;
		    t43 = exp(-t42);
/* Computing 2nd power */
		    d__1 = rho;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t51 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t52 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rhob;
		    t55 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t56 = d__1 * d__1;
		    t60 = t33 * t36;
		    t82 = t45 * .6666666666666667;
		    t93 = pow_dd(&t38, &c_b2);
		    t94 = t93 * .6203504908994;
		    t95 = pow_dd(&t38, &c_b4);
		    t98 = 1 / (t94 + t95 * 10.29581201158544 + 42.7198);
		    t101 = log(t93 * .6203504908994 * t98);
		    t103 = t95 * 1.575246635799487;
		    t107 = atan(.0448998886412873 / (t103 + 13.072));
		    t109 = t95 * .7876233178997433;
/* Computing 2nd power */
		    d__1 = t109 + .409286;
		    t111 = d__1 * d__1;
		    t113 = log(t111 * t98);
		    t117 = 1 / (t94 + t95 * 15.84942278842832 + 101.578);
		    t120 = log(t93 * .6203504908994 * t117);
		    t125 = atan(1.171685277708993 / (t103 + 20.1231));
/* Computing 2nd power */
		    d__1 = t109 + .743294;
		    t128 = d__1 * d__1;
		    t130 = log(t128 * t117);
		    t138 = (rhoa - rhob * 1.) * t38;
		    t139 = t138 + 1.;
		    t140 = pow_dd(&t139, &c_b2);
		    t143 = 1. - t138 * 1.;
		    t144 = pow_dd(&t143, &c_b2);
		    s1 = t5 * -.74442058907928 - t7 * .003024 * sigmaaa / (
			    t10 * .0252 * t11 + 1.) - t19 * .74442058907928;
		    s2 = s1 - t21 * .003024 * sigmabb / (t24 * .0252 * t25 + 
			    1.);
		    zk[i__] = s2 - t36 * .1593432 * rhoa * t39 - t43 * 
			    .0052583256 * t36 / t47 / t45 / rho * (rhoa * 
			    rhob * (t52 * 36.46239897876478 * t51 + t56 * 
			    36.46239897876478 * t55 + (2.611111111111111 - 
			    t33 * .09850555555555556 - t60 * 
			    .1357222222222222) * sigma - (2.5 - t33 * 
			    .01407222222222222 - t60 * .01938888888888889) * 
			    1. * (sigmaaa + sigmabb) - (t42 + t60 * .349 - 
			    11.) * .1111111111111111 * (rhoa * t38 * sigmaaa 
			    + t39 * sigmabb)) - t45 * .6666666666666667 * 
			    sigma + (t82 - t51 * 1.) * sigmabb + (t82 - t55 * 
			    1.) * sigmaaa) + rho * .19 * (t101 * .0310907 + 
			    t107 * 20.5219729378375 + t113 * 
			    .004431373767749538 + (t120 * .01554535 + t125 * 
			    .6188180297906063 + t130 * .002667310007273315 - 
			    t101 * .0310907 - t107 * 20.5219729378375 - t113 *
			     .004431373767749538) * 1.923661050931536 * (t140 
			    * t139 + t144 * t143 - 2.));
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
		    t16 = 1 / rhob;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t21 = t17 * .6203504908994 + t19 * 15.84942278842832 + 
			    101.578;
		    t22 = 1 / t21;
		    t25 = log(t17 * .6203504908994 * t22);
		    t28 = t19 * 1.575246635799487 + 20.1231;
		    t31 = atan(1.171685277708993 / t28);
		    t34 = t19 * .7876233178997433 + .743294;
/* Computing 2nd power */
		    d__1 = t34;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t6 * .003024 * t13 + 
			    rhob * .19 * (t25 * .01554535 + t31 * 
			    .6188180297906063 + t37 * .002667310007273315);
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = rhob;
		    t43 = d__1 * d__1;
		    t45 = 1 / t2 / t43;
/* Computing 2nd power */
		    d__1 = t12;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
/* Computing 2nd power */
		    d__1 = t2;
		    t55 = d__1 * d__1;
		    t60 = 1 / t55 / t43;
		    t63 = sqrt(sigmabb * t60 + 1.);
		    t64 = 1 / t63;
/* Computing 2nd power */
		    d__1 = t17;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t77 = 1 / t43;
/* Computing 2nd power */
		    d__1 = t21;
		    t80 = d__1 * d__1;
		    t81 = 1 / t80;
/* Computing 2nd power */
		    d__1 = t19;
		    t85 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t85;
		    t86 = d__1 * d__1;
		    t88 = 1 / t86 / t19;
		    t89 = t88 * t77;
		    t91 = t75 * -.2067834969664667 * t77 - t89 * 
			    2.641570464738054;
/* Computing 2nd power */
		    d__1 = t28;
		    t99 = d__1 * d__1;
		    t100 = 1 / t99;
		    vrhob[i__] = t2 * -.99256078543904 + t45 * .004032 * 
			    sigmabb * t13 + t6 * .003024 * t50 * (t7 * -.0336 
			    * t45 * t9 - sigmabb * .0336 / t55 / t43 / rhob * 
			    t64) + t25 * .0029536165 + t31 * 
			    .1175754256602152 + t37 * 5.067889013819299e-4 + 
			    rhob * .19 * ((t75 * -.2067834969664667 * t22 * 
			    t77 - t17 * .6203504908994 * t81 * t91) * 
			    .02505897912236993 / t17 * t21 + t100 * 
			    .1903580477513215 * t88 * t77 / (t100 * 
			    1.37284639 + 1.) + (t34 * -.2625411059665811 * 
			    t22 * t89 - t35 * 1. * t81 * t91) * 
			    .002667310007273315 / t35 * t21);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = t5 * -.003024 * t13 + t6 * .003024 * t50 *
			     (.0126 / t7 * t5 * t9 + t60 * .0126 * t64);
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
		    t16 = 1 / rhoa;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t21 = t17 * .6203504908994 + t19 * 15.84942278842832 + 
			    101.578;
		    t22 = 1 / t21;
		    t25 = log(t17 * .6203504908994 * t22);
		    t28 = t19 * 1.575246635799487 + 20.1231;
		    t31 = atan(1.171685277708993 / t28);
		    t34 = t19 * .7876233178997433 + .743294;
/* Computing 2nd power */
		    d__1 = t34;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t6 * .003024 * t13 + 
			    rhoa * .19 * (t25 * .01554535 + t31 * 
			    .6188180297906063 + t37 * .002667310007273315);
/* Computing 2nd power */
		    d__1 = rhoa;
		    t43 = d__1 * d__1;
		    t45 = 1 / t2 / t43;
/* Computing 2nd power */
		    d__1 = t12;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
/* Computing 2nd power */
		    d__1 = t2;
		    t55 = d__1 * d__1;
		    t60 = 1 / t55 / t43;
		    t63 = sqrt(sigmaaa * t60 + 1.);
		    t64 = 1 / t63;
/* Computing 2nd power */
		    d__1 = t17;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t77 = 1 / t43;
/* Computing 2nd power */
		    d__1 = t21;
		    t80 = d__1 * d__1;
		    t81 = 1 / t80;
/* Computing 2nd power */
		    d__1 = t19;
		    t85 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t85;
		    t86 = d__1 * d__1;
		    t88 = 1 / t86 / t19;
		    t89 = t88 * t77;
		    t91 = t75 * -.2067834969664667 * t77 - t89 * 
			    2.641570464738054;
/* Computing 2nd power */
		    d__1 = t28;
		    t99 = d__1 * d__1;
		    t100 = 1 / t99;
		    vrhoa[i__] = t2 * -.99256078543904 + t45 * .004032 * 
			    sigmaaa * t13 + t6 * .003024 * t50 * (t7 * -.0336 
			    * t45 * t9 - sigmaaa * .0336 / t55 / t43 / rhoa * 
			    t64) + t25 * .0029536165 + t31 * 
			    .1175754256602152 + t37 * 5.067889013819299e-4 + 
			    rhoa * .19 * ((t75 * -.2067834969664667 * t22 * 
			    t77 - t17 * .6203504908994 * t81 * t91) * 
			    .02505897912236993 / t17 * t21 + t100 * 
			    .1903580477513215 * t88 * t77 / (t100 * 
			    1.37284639 + 1.) + (t34 * -.2625411059665811 * 
			    t22 * t89 - t35 * 1. * t81 * t91) * 
			    .002667310007273315 / t35 * t21);
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = t5 * -.003024 * t13 + t6 * .003024 * t50 *
			     (.0126 / t7 * t5 * t9 + t60 * .0126 * t64);
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
		    t32 = pow_dd(&rho, &c_b2);
		    t33 = 1 / t32;
		    t35 = t33 * .349 + 1.;
		    t36 = 1 / t35;
		    t37 = t36 * rhoa;
		    t38 = 1 / rho;
		    t39 = rhob * t38;
		    t42 = t33 * .2533;
		    t43 = exp(-t42);
		    t44 = t43 * t36;
/* Computing 2nd power */
		    d__1 = rho;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t47 = d__1 * d__1;
		    t49 = 1 / t47 / t45 / rho;
		    t50 = rhoa * rhob;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t51 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t52 = d__1 * d__1;
		    t53 = t52 * t51;
/* Computing 2nd power */
		    d__1 = rhob;
		    t55 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t56 = d__1 * d__1;
		    t57 = t56 * t55;
		    t60 = t33 * t36;
		    t62 = 2.611111111111111 - t33 * .09850555555555556 - t60 *
			     .1357222222222222;
		    t67 = sigmaaa + sigmabb;
		    t71 = t42 + t60 * .349 - 11.;
		    t75 = rhoa * t38 * sigmaaa + t39 * sigmabb;
		    t78 = t53 * 36.46239897876478 + t57 * 36.46239897876478 + 
			    t62 * sigma - (2.5 - t33 * .01407222222222222 - 
			    t60 * .01938888888888889) * 1. * t67 - t71 * 
			    .1111111111111111 * t75;
		    t82 = t45 * .6666666666666667;
		    t83 = t51 * 1.;
		    t86 = t55 * 1.;
		    t89 = t50 * t78 - t45 * .6666666666666667 * sigma + (t82 
			    - t83) * sigmabb + (t82 - t86) * sigmaaa;
		    t93 = pow_dd(&t38, &c_b2);
		    t94 = t93 * .6203504908994;
		    t95 = pow_dd(&t38, &c_b4);
		    t97 = t94 + t95 * 10.29581201158544 + 42.7198;
		    t98 = 1 / t97;
		    t101 = log(t93 * .6203504908994 * t98);
		    t103 = t95 * 1.575246635799487;
		    t104 = t103 + 13.072;
		    t107 = atan(.0448998886412873 / t104);
		    t109 = t95 * .7876233178997433;
		    t110 = t109 + .409286;
/* Computing 2nd power */
		    d__1 = t110;
		    t111 = d__1 * d__1;
		    t113 = log(t111 * t98);
		    t116 = t94 + t95 * 15.84942278842832 + 101.578;
		    t117 = 1 / t116;
		    t120 = log(t93 * .6203504908994 * t117);
		    t122 = t103 + 20.1231;
		    t125 = atan(1.171685277708993 / t122);
		    t127 = t109 + .743294;
/* Computing 2nd power */
		    d__1 = t127;
		    t128 = d__1 * d__1;
		    t130 = log(t128 * t117);
		    t135 = t120 * .01554535 + t125 * .6188180297906063 + t130 
			    * .002667310007273315 - t101 * .0310907 - t107 * 
			    20.5219729378375 - t113 * .004431373767749538;
		    t137 = rhoa - rhob * 1.;
		    t138 = t137 * t38;
		    t139 = t138 + 1.;
		    t140 = pow_dd(&t139, &c_b2);
		    t143 = 1. - t138 * 1.;
		    t144 = pow_dd(&t143, &c_b2);
		    t146 = t140 * t139 + t144 * t143 - 2.;
		    t147 = t135 * t146;
		    zk[i__] = t5 * -.74442058907928 - t8 * .003024 * t15 - 
			    t19 * .74442058907928 - t22 * .003024 * t29 - t37 
			    * .1593432 * t39 - t44 * .0052583256 * t49 * t89 
			    + rho * .19 * (t101 * .0310907 + t107 * 
			    20.5219729378375 + t113 * .004431373767749538 + 
			    t147 * 1.923661050931536);
		    t154 = 1 / t4 / t51;
/* Computing 2nd power */
		    d__1 = t14;
		    t158 = d__1 * d__1;
		    t159 = 1 / t158;
		    t167 = 1 / t53;
		    t170 = sqrt(sigmaaa * t167 + 1.);
		    t171 = 1 / t170;
		    t184 = t71 * t38;
		    t195 = rho * t135;
		    t200 = t140 * 1.333333333333333 * t38 - t144 * 
			    1.333333333333333 * t38;
/* Computing 2nd power */
		    d__1 = t35;
		    t203 = d__1 * d__1;
		    t204 = 1 / t203;
		    t210 = t204 * .0185369256 * rhoa * rhob / t32 / t45;
		    t211 = 1 / t45;
		    t212 = rhob * t211;
		    t214 = t37 * .1593432 * t212;
/* Computing 2nd power */
		    d__1 = t45;
		    t215 = d__1 * d__1;
		    t217 = 1 / t215 / rho;
		    t221 = t217 * 4.4397795816e-4 * t43 * t36 * t89;
		    t225 = t43 * 6.117185448e-4 * t204 * t217 * t89;
		    t230 = t44 * .0192805272 / t47 / t215 * t89;
		    t232 = 1 / t32 / rho;
		    t234 = t232 * t36;
		    t238 = 1 / t47 / rho * t204;
		    t273 = t44 * .0052583256 * t49 * (t50 * ((t232 * 
			    .03283518518518519 + t234 * .04524074074074074 - 
			    t238 * .01578901851851852) * sigma - (t232 * 
			    .004690740740740741 + t234 * .006462962962962963 
			    - t238 * .002255574074074074) * 1. * t67 - (t232 *
			     -.08443333333333333 - t234 * .1163333333333333 + 
			    t238 * .04060033333333333) * .1111111111111111 * 
			    t75 - t71 * .1111111111111111 * (rhoa * -1. * 
			    t211 * sigmaaa - t212 * 1. * sigmabb)) - rho * 
			    1.333333333333333 * sigma + rho * 
			    1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa);
		    t274 = t101 * .005907233;
		    t275 = t107 * 3.899174858189126;
		    t276 = t113 * 8.419610158724123e-4;
		    t277 = t147 * .3654955996769919;
/* Computing 2nd power */
		    d__1 = t93;
		    t278 = d__1 * d__1;
		    t279 = 1 / t278;
/* Computing 2nd power */
		    d__1 = t97;
		    t283 = d__1 * d__1;
		    t284 = 1 / t283;
		    t287 = t279 * .2067834969664667 * t211;
/* Computing 2nd power */
		    d__1 = t95;
		    t288 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t288;
		    t289 = d__1 * d__1;
		    t291 = 1 / t289 / t95;
		    t292 = t291 * t211;
		    t294 = -t287 - t292 * 1.715968668597574;
		    t298 = 1 / t93;
		    t300 = (t279 * -.2067834969664667 * t98 * t211 - t93 * 
			    .6203504908994 * t284 * t294) * t298 * t97;
/* Computing 2nd power */
		    d__1 = t104;
		    t302 = d__1 * d__1;
		    t303 = 1 / t302;
		    t309 = t303 * t291 * t211 / (t303 * .002016 + 1.);
		    t320 = (t110 * -.2625411059665811 * t98 * t292 - t111 * 
			    1. * t284 * t294) / t111 * t97;
/* Computing 2nd power */
		    d__1 = t116;
		    t325 = d__1 * d__1;
		    t326 = 1 / t325;
		    t329 = -t287 - t292 * 2.641570464738054;
/* Computing 2nd power */
		    d__1 = t122;
		    t336 = d__1 * d__1;
		    t337 = 1 / t336;
		    t373 = rho * .19 * (t300 * .05011795824473985 + t309 * 
			    .2419143800947354 + t320 * .004431373767749538 + (
			    (t279 * -.2067834969664667 * t117 * t211 - t93 * 
			    .6203504908994 * t326 * t329) * 
			    .02505897912236993 * t298 * t116 + t337 * 
			    .1903580477513215 * t291 * t211 / (t337 * 
			    1.37284639 + 1.) + (t127 * -.2625411059665811 * 
			    t117 * t292 - t128 * 1. * t326 * t329) * 
			    .002667310007273315 / t128 * t116 - t300 * 
			    .05011795824473985 - t309 * .2419143800947354 - 
			    t320 * .004431373767749538) * 1.923661050931536 * 
			    t146 + t135 * 1.923661050931536 * (t140 * 
			    -1.333333333333333 * t137 * t211 + t144 * 
			    1.333333333333333 * t137 * t211));
		    vrhoa[i__] = t4 * -.99256078543904 + t154 * .004032 * 
			    sigmaaa * t15 + t8 * .003024 * t159 * (t9 * 
			    -.0336 * t154 * t11 - sigmaaa * .0336 / t52 / t51 
			    / rhoa * t171) - t36 * .1593432 * rhob * t38 - 
			    t44 * .0052583256 * t49 * (rhob * t78 + t50 * (
			    t52 * 97.23306394337274 * rhoa - t184 * 
			    .1111111111111111 * sigmaaa) - rhoa * 2. * 
			    sigmabb) + t195 * .3654955996769919 * t200 - t210 
			    + t214 - t221 - t225 + t230 - t273 + t274 + t275 
			    + t276 + t277 + t373;
		    t376 = 1 / t18 / t55;
/* Computing 2nd power */
		    d__1 = t28;
		    t380 = d__1 * d__1;
		    t381 = 1 / t380;
		    t389 = 1 / t57;
		    t392 = sqrt(sigmabb * t389 + 1.);
		    t393 = 1 / t392;
		    vrhob[i__] = t18 * -.99256078543904 + t376 * .004032 * 
			    sigmabb * t29 + t22 * .003024 * t381 * (t23 * 
			    -.0336 * t376 * t25 - sigmabb * .0336 / t56 / t55 
			    / rhob * t393) - t37 * .1593432 * t38 - t44 * 
			    .0052583256 * t49 * (rhoa * t78 + t50 * (t56 * 
			    97.23306394337274 * rhob - t184 * 
			    .1111111111111111 * sigmabb) - rhob * 2. * 
			    sigmaaa) - t195 * .3654955996769919 * t200 - t210 
			    + t214 - t221 - t225 + t230 - t273 + t274 + t275 
			    + t276 + t277 + t373;
		    t430 = t33 * .01407222222222222;
		    t431 = t60 * .01938888888888889;
		    t445 = t44 * t49 * (t50 * t62 - t45 * .6666666666666667);
		    t446 = t445 * .0052583256;
		    vsigmaaa[i__] = t7 * -.003024 * t15 + t8 * .003024 * t159 
			    * (.0126 / t9 * t7 * t11 + t167 * .0126 * t171) - 
			    t44 * .0052583256 * t49 * (t50 * (t430 - 2.5 + 
			    t431 - t71 * .1111111111111111 * rhoa * t38) + 
			    t82 - t86) - t446;
		    vsigmaab[i__] = t445 * -.0105166512;
		    vsigmabb[i__] = t21 * -.003024 * t29 + t22 * .003024 * 
			    t381 * (.0126 / t23 * t21 * t25 + t389 * .0126 * 
			    t393) - t44 * .0052583256 * t49 * (t50 * (t430 - 
			    2.5 + t431 - t71 * .1111111111111111 * rhob * t38)
			     + t82 - t83) - t446;
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
		    t16 = 1 / rhob;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t21 = t17 * .6203504908994 + t19 * 15.84942278842832 + 
			    101.578;
		    t22 = 1 / t21;
		    t25 = log(t17 * .6203504908994 * t22);
		    t28 = t19 * 1.575246635799487 + 20.1231;
		    t31 = atan(1.171685277708993 / t28);
		    t34 = t19 * .7876233178997433 + .743294;
/* Computing 2nd power */
		    d__1 = t34;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t6 * .003024 * t13 + 
			    rhob * .19 * (t25 * .01554535 + t31 * 
			    .6188180297906063 + t37 * .002667310007273315);
		    vrhoa[i__] = 0.;
/* Computing 2nd power */
		    d__1 = rhob;
		    t43 = d__1 * d__1;
		    t45 = 1 / t2 / t43;
		    t46 = t45 * sigmabb;
/* Computing 2nd power */
		    d__1 = t12;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t54 = t43 * rhob;
/* Computing 2nd power */
		    d__1 = t2;
		    t55 = d__1 * d__1;
		    t60 = 1 / t55 / t43;
		    t62 = sigmabb * t60 + 1.;
		    t63 = sqrt(t62);
		    t64 = 1 / t63;
		    t67 = t7 * -.0336 * t45 * t9 - sigmabb * .0336 / t55 / 
			    t54 * t64;
		    t68 = t50 * t67;
/* Computing 2nd power */
		    d__1 = t17;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t76 = t75 * t22;
		    t77 = 1 / t43;
/* Computing 2nd power */
		    d__1 = t21;
		    t80 = d__1 * d__1;
		    t81 = 1 / t80;
		    t82 = t17 * t81;
/* Computing 2nd power */
		    d__1 = t19;
		    t85 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t85;
		    t86 = d__1 * d__1;
		    t87 = t86 * t19;
		    t88 = 1 / t87;
		    t89 = t88 * t77;
		    t91 = t75 * -.2067834969664667 * t77 - t89 * 
			    2.641570464738054;
		    t94 = t76 * -.2067834969664667 * t77 - t82 * 
			    .6203504908994 * t91;
		    t95 = 1 / t17;
		    t96 = t94 * t95;
		    t97 = t96 * t21;
/* Computing 2nd power */
		    d__1 = t28;
		    t99 = d__1 * d__1;
		    t100 = 1 / t99;
		    t101 = t100 * t88;
		    t103 = t100 * 1.37284639 + 1.;
		    t104 = 1 / t103;
		    t106 = t101 * t77 * t104;
		    t108 = t34 * t22;
		    t111 = t35 * t81;
		    t114 = t108 * -.2625411059665811 * t89 - t111 * 1. * t91;
		    t115 = 1 / t35;
		    t116 = t114 * t115;
		    t117 = t116 * t21;
		    vrhob[i__] = t2 * -.99256078543904 + t46 * .004032 * t13 
			    + t6 * .003024 * t68 + t25 * .0029536165 + t31 * 
			    .1175754256602152 + t37 * 5.067889013819299e-4 + 
			    rhob * .19 * (t97 * .02505897912236993 + t106 * 
			    .1903580477513215 + t117 * .002667310007273315);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t130 = .0126 / t7 * t5 * t9 + t60 * .0126 * t64;
		    vsigmabb[i__] = t5 * -.003024 * t13 + t6 * .003024 * t50 *
			     t130;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t137 = 1 / t2 / t54;
		    t144 = 1 / t49 / t12;
/* Computing 2nd power */
		    d__1 = t67;
		    t145 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t43;
		    t152 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t158 = d__1 * d__1;
		    t164 = 1 / t63 / t62;
		    t175 = 1 / t74 / t16;
		    t177 = 1 / t152;
		    t178 = t175 * t22 * t177;
		    t184 = 1 / t54;
		    t188 = 1 / t80 / t21;
/* Computing 2nd power */
		    d__1 = t91;
		    t190 = d__1 * d__1;
		    t198 = 1 / t87 / t16;
		    t199 = t198 * t177;
		    t201 = t88 * t184;
		    t203 = t175 * -.1378556646443111 * t177 + t75 * 
			    .4135669939329333 * t184 - t199 * 
			    2.201308720615045 + t201 * 5.283140929476108;
		    t221 = t177 * t104;
/* Computing 2nd power */
		    d__1 = t99;
		    t230 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t103;
		    t234 = d__1 * d__1;
		    s1 = -.3308535951463467 / t55 - t137 * .009408 * sigmabb *
			     t13 - t46 * .008064 * t68 - t6 * .006048 * t144 *
			     t145;
		    s2 = s1 + t6 * .003024 * t50 * (t7 * .0784 * t137 * t9 + 
			    sigmabb * .168 / t55 / t152 * t64 - t158 * .0448 /
			     t2 / t152 / t54 * t164) + t97 * 
			    .009522412066500572;
		    s3 = s2 + t106 * .07233605814550218;
		    s4 = s3;
		    s5 = t117 * .00101357780276386 + rhob * .19 * ((t178 * 
			    -.1378556646443111 + t75 * .4135669939329333 * 
			    t81 * t77 * t91 + t76 * .4135669939329333 * t184 
			    + t17 * 1.2407009817988 * t188 * t190 - t82 * 
			    .6203504908994 * t203) * .02505897912236993 * t95 
			    * t21 + t94 * .008352993040789975 / t17 / t16 * 
			    t21 * t77 + t96 * .02505897912236993 * t91 + 
			    .09995362477254242 / t99 / t28 * t175 * t221 + 
			    t100 * .1586317064594346 * t198 * t221 - t101 * 
			    .3807160955026431 * t184 * t104 - 
			    .1372209729363994 / t230 / t28 * t175 * t177 / 
			    t234 + (t178 * .03446391616107778 + t34 * 
			    .5250822119331622 * t81 * t89 * t91 - t108 * 
			    .2187842549721509 * t199 + t108 * 
			    .5250822119331622 * t201 + t35 * 2. * t188 * t190 
			    - t111 * 1. * t203) * .002667310007273315 * t115 *
			     t21 + t114 * 7.002785192652656e-4 / t35 / t34 * 
			    t21 * t88 * t77 + t116 * .002667310007273315 * 
			    t91);
		    v2rhob2[i__] = s4 + s5;
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t130;
		    t272 = d__1 * d__1;
		    v2sigmabb2[i__] = t5 * .006048 * t50 * t130 - t6 * 
			    .006048 * t144 * t272 + t6 * .003024 * t50 * (
			    -.0063 / t7 / sigmabb * t5 * t9 + .0063 / sigmabb 
			    * t60 * t64 - .0063 / t2 / t152 / rhob * t164);
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
		    t16 = 1 / rhoa;
		    t17 = pow_dd(&t16, &c_b2);
		    t19 = pow_dd(&t16, &c_b4);
		    t21 = t17 * .6203504908994 + t19 * 15.84942278842832 + 
			    101.578;
		    t22 = 1 / t21;
		    t25 = log(t17 * .6203504908994 * t22);
		    t28 = t19 * 1.575246635799487 + 20.1231;
		    t31 = atan(1.171685277708993 / t28);
		    t34 = t19 * .7876233178997433 + .743294;
/* Computing 2nd power */
		    d__1 = t34;
		    t35 = d__1 * d__1;
		    t37 = log(t35 * t22);
		    zk[i__] = t3 * -.74442058907928 - t6 * .003024 * t13 + 
			    rhoa * .19 * (t25 * .01554535 + t31 * 
			    .6188180297906063 + t37 * .002667310007273315);
/* Computing 2nd power */
		    d__1 = rhoa;
		    t43 = d__1 * d__1;
		    t45 = 1 / t2 / t43;
		    t46 = t45 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t12;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t54 = t43 * rhoa;
/* Computing 2nd power */
		    d__1 = t2;
		    t55 = d__1 * d__1;
		    t60 = 1 / t55 / t43;
		    t62 = sigmaaa * t60 + 1.;
		    t63 = sqrt(t62);
		    t64 = 1 / t63;
		    t67 = t7 * -.0336 * t45 * t9 - sigmaaa * .0336 / t55 / 
			    t54 * t64;
		    t68 = t50 * t67;
/* Computing 2nd power */
		    d__1 = t17;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t76 = t75 * t22;
		    t77 = 1 / t43;
/* Computing 2nd power */
		    d__1 = t21;
		    t80 = d__1 * d__1;
		    t81 = 1 / t80;
		    t82 = t17 * t81;
/* Computing 2nd power */
		    d__1 = t19;
		    t85 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t85;
		    t86 = d__1 * d__1;
		    t87 = t86 * t19;
		    t88 = 1 / t87;
		    t89 = t88 * t77;
		    t91 = t75 * -.2067834969664667 * t77 - t89 * 
			    2.641570464738054;
		    t94 = t76 * -.2067834969664667 * t77 - t82 * 
			    .6203504908994 * t91;
		    t95 = 1 / t17;
		    t96 = t94 * t95;
		    t97 = t96 * t21;
/* Computing 2nd power */
		    d__1 = t28;
		    t99 = d__1 * d__1;
		    t100 = 1 / t99;
		    t101 = t100 * t88;
		    t103 = t100 * 1.37284639 + 1.;
		    t104 = 1 / t103;
		    t106 = t101 * t77 * t104;
		    t108 = t34 * t22;
		    t111 = t35 * t81;
		    t114 = t108 * -.2625411059665811 * t89 - t111 * 1. * t91;
		    t115 = 1 / t35;
		    t116 = t114 * t115;
		    t117 = t116 * t21;
		    vrhoa[i__] = t2 * -.99256078543904 + t46 * .004032 * t13 
			    + t6 * .003024 * t68 + t25 * .0029536165 + t31 * 
			    .1175754256602152 + t37 * 5.067889013819299e-4 + 
			    rhoa * .19 * (t97 * .02505897912236993 + t106 * 
			    .1903580477513215 + t117 * .002667310007273315);
		    vrhob[i__] = 0.;
		    t130 = .0126 / t7 * t5 * t9 + t60 * .0126 * t64;
		    vsigmaaa[i__] = t5 * -.003024 * t13 + t6 * .003024 * t50 *
			     t130;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t137 = 1 / t2 / t54;
		    t144 = 1 / t49 / t12;
/* Computing 2nd power */
		    d__1 = t67;
		    t145 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t43;
		    t152 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t158 = d__1 * d__1;
		    t164 = 1 / t63 / t62;
		    t175 = 1 / t74 / t16;
		    t177 = 1 / t152;
		    t178 = t175 * t22 * t177;
		    t184 = 1 / t54;
		    t188 = 1 / t80 / t21;
/* Computing 2nd power */
		    d__1 = t91;
		    t190 = d__1 * d__1;
		    t198 = 1 / t87 / t16;
		    t199 = t198 * t177;
		    t201 = t88 * t184;
		    t203 = t175 * -.1378556646443111 * t177 + t75 * 
			    .4135669939329333 * t184 - t199 * 
			    2.201308720615045 + t201 * 5.283140929476108;
		    t221 = t177 * t104;
/* Computing 2nd power */
		    d__1 = t99;
		    t230 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t103;
		    t234 = d__1 * d__1;
		    s1 = -.3308535951463467 / t55 - t137 * .009408 * sigmaaa *
			     t13 - t46 * .008064 * t68 - t6 * .006048 * t144 *
			     t145;
		    s2 = s1 + t6 * .003024 * t50 * (t7 * .0784 * t137 * t9 + 
			    sigmaaa * .168 / t55 / t152 * t64 - t158 * .0448 /
			     t2 / t152 / t54 * t164) + t97 * 
			    .009522412066500572;
		    s3 = s2 + t106 * .07233605814550218;
		    s4 = s3;
		    s5 = t117 * .00101357780276386 + rhoa * .19 * ((t178 * 
			    -.1378556646443111 + t75 * .4135669939329333 * 
			    t81 * t77 * t91 + t76 * .4135669939329333 * t184 
			    + t17 * 1.2407009817988 * t188 * t190 - t82 * 
			    .6203504908994 * t203) * .02505897912236993 * t95 
			    * t21 + t94 * .008352993040789975 / t17 / t16 * 
			    t21 * t77 + t96 * .02505897912236993 * t91 + 
			    .09995362477254242 / t99 / t28 * t175 * t221 + 
			    t100 * .1586317064594346 * t198 * t221 - t101 * 
			    .3807160955026431 * t184 * t104 - 
			    .1372209729363994 / t230 / t28 * t175 * t177 / 
			    t234 + (t178 * .03446391616107778 + t34 * 
			    .5250822119331622 * t81 * t89 * t91 - t108 * 
			    .2187842549721509 * t199 + t108 * 
			    .5250822119331622 * t201 + t35 * 2. * t188 * t190 
			    - t111 * 1. * t203) * .002667310007273315 * t115 *
			     t21 + t114 * 7.002785192652656e-4 / t35 / t34 * 
			    t21 * t88 * t77 + t116 * .002667310007273315 * 
			    t91);
		    v2rhoa2[i__] = s4 + s5;
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t130;
		    t272 = d__1 * d__1;
		    v2sigmaaa2[i__] = t5 * .006048 * t50 * t130 - t6 * 
			    .006048 * t144 * t272 + t6 * .003024 * t50 * (
			    -.0063 / t7 / sigmaaa * t5 * t9 + .0063 / sigmaaa 
			    * t60 * t64 - .0063 / t2 / t152 / rhoa * t164);
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
		    t32 = pow_dd(&rho, &c_b2);
		    t33 = 1 / t32;
		    t35 = t33 * .349 + 1.;
		    t36 = 1 / t35;
		    t37 = t36 * rhoa;
		    t38 = 1 / rho;
		    t39 = rhob * t38;
		    t42 = t33 * .2533;
		    t43 = exp(-t42);
		    t44 = t43 * t36;
/* Computing 2nd power */
		    d__1 = rho;
		    t45 = d__1 * d__1;
		    t46 = t45 * rho;
/* Computing 2nd power */
		    d__1 = t32;
		    t47 = d__1 * d__1;
		    t49 = 1 / t47 / t46;
		    t50 = rhoa * rhob;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t51 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t52 = d__1 * d__1;
		    t53 = t52 * t51;
		    t54 = t53 * 36.46239897876478;
/* Computing 2nd power */
		    d__1 = rhob;
		    t55 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t56 = d__1 * d__1;
		    t57 = t56 * t55;
		    t58 = t57 * 36.46239897876478;
		    t60 = t33 * t36;
		    t62 = 2.611111111111111 - t33 * .09850555555555556 - t60 *
			     .1357222222222222;
		    t63 = t62 * sigma;
		    t67 = sigmaaa + sigmabb;
		    t69 = (2.5 - t33 * .01407222222222222 - t60 * 
			    .01938888888888889) * 1. * t67;
		    t71 = t42 + t60 * .349 - 11.;
		    t75 = rhoa * t38 * sigmaaa + t39 * sigmabb;
		    t77 = t71 * .1111111111111111 * t75;
		    t78 = t54 + t58 + t63 - t69 - t77;
		    t82 = t45 * .6666666666666667;
		    t83 = t51 * 1.;
		    t86 = t55 * 1.;
		    t89 = t50 * t78 - t45 * .6666666666666667 * sigma + (t82 
			    - t83) * sigmabb + (t82 - t86) * sigmaaa;
		    t93 = pow_dd(&t38, &c_b2);
		    t94 = t93 * .6203504908994;
		    t95 = pow_dd(&t38, &c_b4);
		    t97 = t94 + t95 * 10.29581201158544 + 42.7198;
		    t98 = 1 / t97;
		    t101 = log(t93 * .6203504908994 * t98);
		    t103 = t95 * 1.575246635799487;
		    t104 = t103 + 13.072;
		    t107 = atan(.0448998886412873 / t104);
		    t109 = t95 * .7876233178997433;
		    t110 = t109 + .409286;
/* Computing 2nd power */
		    d__1 = t110;
		    t111 = d__1 * d__1;
		    t113 = log(t111 * t98);
		    t116 = t94 + t95 * 15.84942278842832 + 101.578;
		    t117 = 1 / t116;
		    t120 = log(t93 * .6203504908994 * t117);
		    t122 = t103 + 20.1231;
		    t125 = atan(1.171685277708993 / t122);
		    t127 = t109 + .743294;
/* Computing 2nd power */
		    d__1 = t127;
		    t128 = d__1 * d__1;
		    t130 = log(t128 * t117);
		    t135 = t120 * .01554535 + t125 * .6188180297906063 + t130 
			    * .002667310007273315 - t101 * .0310907 - t107 * 
			    20.5219729378375 - t113 * .004431373767749538;
		    t137 = rhoa - rhob * 1.;
		    t138 = t137 * t38;
		    t139 = t138 + 1.;
		    t140 = pow_dd(&t139, &c_b2);
		    t143 = 1. - t138 * 1.;
		    t144 = pow_dd(&t143, &c_b2);
		    t146 = t140 * t139 + t144 * t143 - 2.;
		    t147 = t135 * t146;
		    zk[i__] = t5 * -.74442058907928 - t8 * .003024 * t15 - 
			    t19 * .74442058907928 - t22 * .003024 * t29 - t37 
			    * .1593432 * t39 - t44 * .0052583256 * t49 * t89 
			    + rho * .19 * (t101 * .0310907 + t107 * 
			    20.5219729378375 + t113 * .004431373767749538 + 
			    t147 * 1.923661050931536);
		    t154 = 1 / t4 / t51;
		    t155 = t154 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t14;
		    t158 = d__1 * d__1;
		    t159 = 1 / t158;
		    t163 = t51 * rhoa;
		    t165 = 1 / t52 / t163;
		    t167 = 1 / t53;
		    t169 = sigmaaa * t167 + 1.;
		    t170 = sqrt(t169);
		    t171 = 1 / t170;
		    t174 = t9 * -.0336 * t154 * t11 - sigmaaa * .0336 * t165 *
			     t171;
		    t175 = t159 * t174;
		    t178 = t36 * rhob;
		    t182 = t52 * rhoa;
		    t184 = t71 * t38;
		    t187 = t182 * 97.23306394337274 - t184 * 
			    .1111111111111111 * sigmaaa;
		    t191 = rhob * t78 + t50 * t187 - rhoa * 2. * sigmabb;
		    t195 = rho * t135;
		    t200 = t140 * 1.333333333333333 * t38 - t144 * 
			    1.333333333333333 * t38;
/* Computing 2nd power */
		    d__1 = t35;
		    t203 = d__1 * d__1;
		    t204 = 1 / t203;
		    t205 = t204 * rhoa;
		    t207 = 1 / t32 / t45;
		    t210 = t205 * .0185369256 * rhob * t207;
		    t211 = 1 / t45;
		    t212 = rhob * t211;
		    t214 = t37 * .1593432 * t212;
/* Computing 2nd power */
		    d__1 = t45;
		    t215 = d__1 * d__1;
		    t216 = t215 * rho;
		    t217 = 1 / t216;
		    t218 = t217 * t43;
		    t219 = t36 * t89;
		    t221 = t218 * 4.4397795816e-4 * t219;
		    t222 = t43 * t204;
		    t225 = t222 * 6.117185448e-4 * t217 * t89;
		    t227 = 1 / t47 / t215;
		    t230 = t44 * .0192805272 * t227 * t89;
		    t232 = 1 / t32 / rho;
		    t234 = t232 * t36;
		    t238 = 1 / t47 / rho * t204;
		    t240 = t232 * .03283518518518519 + t234 * 
			    .04524074074074074 - t238 * .01578901851851852;
		    t251 = t232 * -.08443333333333333 - t234 * 
			    .1163333333333333 + t238 * .04060033333333333;
		    t259 = rhoa * -1. * t211 * sigmaaa - t212 * 1. * sigmabb;
		    t262 = t240 * sigma - (t232 * .004690740740740741 + t234 *
			     .006462962962962963 - t238 * .002255574074074074)
			     * 1. * t67 - t251 * .1111111111111111 * t75 - 
			    t71 * .1111111111111111 * t259;
		    t270 = t50 * t262 - rho * 1.333333333333333 * sigma + rho 
			    * 1.333333333333333 * sigmabb + rho * 
			    1.333333333333333 * sigmaaa;
		    t273 = t44 * .0052583256 * t49 * t270;
		    t274 = t101 * .005907233;
		    t275 = t107 * 3.899174858189126;
		    t276 = t113 * 8.419610158724123e-4;
		    t277 = t147 * .3654955996769919;
/* Computing 2nd power */
		    d__1 = t93;
		    t278 = d__1 * d__1;
		    t279 = 1 / t278;
		    t280 = t279 * t98;
/* Computing 2nd power */
		    d__1 = t97;
		    t283 = d__1 * d__1;
		    t284 = 1 / t283;
		    t285 = t93 * t284;
		    t287 = t279 * .2067834969664667 * t211;
/* Computing 2nd power */
		    d__1 = t95;
		    t288 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t288;
		    t289 = d__1 * d__1;
		    t290 = t289 * t95;
		    t291 = 1 / t290;
		    t292 = t291 * t211;
		    t294 = -t287 - t292 * 1.715968668597574;
		    t297 = t280 * -.2067834969664667 * t211 - t285 * 
			    .6203504908994 * t294;
		    t298 = 1 / t93;
		    t299 = t297 * t298;
		    t300 = t299 * t97;
/* Computing 2nd power */
		    d__1 = t104;
		    t302 = d__1 * d__1;
		    t303 = 1 / t302;
		    t304 = t303 * t291;
		    t306 = t303 * .002016 + 1.;
		    t307 = 1 / t306;
		    t309 = t304 * t211 * t307;
		    t311 = t110 * t98;
		    t314 = t111 * t284;
		    t317 = t311 * -.2625411059665811 * t292 - t314 * 1. * 
			    t294;
		    t318 = 1 / t111;
		    t319 = t317 * t318;
		    t320 = t319 * t97;
		    t322 = t279 * t117;
/* Computing 2nd power */
		    d__1 = t116;
		    t325 = d__1 * d__1;
		    t326 = 1 / t325;
		    t327 = t93 * t326;
		    t329 = -t287 - t292 * 2.641570464738054;
		    t332 = t322 * -.2067834969664667 * t211 - t327 * 
			    .6203504908994 * t329;
		    t333 = t332 * t298;
/* Computing 2nd power */
		    d__1 = t122;
		    t336 = d__1 * d__1;
		    t337 = 1 / t336;
		    t338 = t337 * t291;
		    t340 = t337 * 1.37284639 + 1.;
		    t341 = 1 / t340;
		    t345 = t127 * t117;
		    t348 = t128 * t326;
		    t351 = t345 * -.2625411059665811 * t292 - t348 * 1. * 
			    t329;
		    t352 = 1 / t128;
		    t353 = t351 * t352;
		    t359 = t333 * .02505897912236993 * t116 + t338 * 
			    .1903580477513215 * t211 * t341 + t353 * 
			    .002667310007273315 * t116 - t300 * 
			    .05011795824473985 - t309 * .2419143800947354 - 
			    t320 * .004431373767749538;
		    t360 = t359 * t146;
		    t362 = t140 * t137;
		    t365 = t144 * t137;
		    t368 = t362 * -1.333333333333333 * t211 + t365 * 
			    1.333333333333333 * t211;
		    t369 = t135 * t368;
		    t373 = rho * .19 * (t300 * .05011795824473985 + t309 * 
			    .2419143800947354 + t320 * .004431373767749538 + 
			    t360 * 1.923661050931536 + t369 * 
			    1.923661050931536);
		    vrhoa[i__] = t4 * -.99256078543904 + t155 * .004032 * t15 
			    + t8 * .003024 * t175 - t178 * .1593432 * t38 - 
			    t44 * .0052583256 * t49 * t191 + t195 * 
			    .3654955996769919 * t200 - t210 + t214 - t221 - 
			    t225 + t230 - t273 + t274 + t275 + t276 + t277 + 
			    t373;
		    t376 = 1 / t18 / t55;
		    t377 = t376 * sigmabb;
/* Computing 2nd power */
		    d__1 = t28;
		    t380 = d__1 * d__1;
		    t381 = 1 / t380;
		    t385 = t55 * rhob;
		    t387 = 1 / t56 / t385;
		    t389 = 1 / t57;
		    t391 = sigmabb * t389 + 1.;
		    t392 = sqrt(t391);
		    t393 = 1 / t392;
		    t396 = t23 * -.0336 * t376 * t25 - sigmabb * .0336 * t387 
			    * t393;
		    t397 = t381 * t396;
		    t403 = t56 * rhob;
		    t407 = t403 * 97.23306394337274 - t184 * 
			    .1111111111111111 * sigmabb;
		    t411 = rhoa * t78 + t50 * t407 - rhob * 2. * sigmaaa;
		    t415 = -t200;
		    vrhob[i__] = t18 * -.99256078543904 + t377 * .004032 * 
			    t29 + t22 * .003024 * t397 - t37 * .1593432 * t38 
			    - t44 * .0052583256 * t49 * t411 + t195 * 
			    .3654955996769919 * t415 - t210 + t214 - t221 - 
			    t225 + t230 - t273 + t274 + t275 + t276 + t277 + 
			    t373;
		    t420 = 1 / t9;
		    t426 = t420 * .0126 * t7 * t11 + t167 * .0126 * t171;
		    t427 = t159 * t426;
		    t430 = t33 * .01407222222222222;
		    t431 = t60 * .01938888888888889;
		    t432 = t71 * rhoa;
		    t435 = t430 - 2.5 + t431 - t432 * .1111111111111111 * t38;
		    t437 = t50 * t435 + t82 - t86;
		    t443 = t50 * t62 - t45 * .6666666666666667;
		    t445 = t44 * t49 * t443;
		    t446 = t445 * .0052583256;
		    vsigmaaa[i__] = t7 * -.003024 * t15 + t8 * .003024 * t427 
			    - t44 * .0052583256 * t49 * t437 - t446;
		    vsigmaab[i__] = t445 * -.0105166512;
		    t450 = 1 / t23;
		    t456 = t450 * .0126 * t21 * t25 + t389 * .0126 * t393;
		    t457 = t381 * t456;
		    t460 = t71 * rhob;
		    t463 = t430 - 2.5 + t431 - t460 * .1111111111111111 * t38;
		    t465 = t50 * t463 + t82 - t83;
		    vsigmabb[i__] = t21 * -.003024 * t29 + t22 * .003024 * 
			    t457 - t44 * .0052583256 * t49 * t465 - t446;
		    t469 = t309 * .09192746443599945;
		    t470 = t300 * .01904482413300114;
		    t471 = t320 * .001683922031744825;
		    t472 = t360 * .7309911993539838;
		    t473 = t369 * .7309911993539838;
		    t478 = t44 * .0899757936 / t47 / t216 * t89;
		    t483 = t205 * .061789752 * rhob / t32 / t46;
		    t484 = t215 * t45;
		    t485 = 1 / t484;
		    t488 = t222 * .0053015607216 * t485 * t89;
		    t490 = 1 / t32 / t484;
		    t491 = t490 * t43;
		    t494 = t491 * 1.0329887159856e-4 * t204 * t89;
		    t497 = t44 * .0385610544 * t227 * t270;
		    t500 = t222 * .0012234370896 * t217 * t270;
		    t501 = 1 / t46;
		    t502 = rhob * t501;
		    t503 = t37 * t502;
		    t504 = t503 * .3186864;
		    t506 = t207 * t36;
		    t510 = 1 / t47 / t45 * t204;
		    t513 = 1 / t203 / t35;
		    t514 = t501 * t513;
		    t549 = t44 * t49 * (t50 * ((t207 * -.04378024691358025 - 
			    t506 * .06032098765432099 + t510 * 
			    .03157803703703704 - t514 * .003673578308641975) *
			     sigma - (t207 * -.006254320987654321 - t506 * 
			    .008617283950617284 + t510 * .004511148148148148 
			    - t514 * 5.247969012345679e-4) * 1. * t67 - (t207 
			    * .1125777777777778 + t506 * .1551111111111111 - 
			    t510 * .08120066666666667 + t514 * 
			    .009446344222222222) * .1111111111111111 * t75 - 
			    t251 * .2222222222222222 * t259 - t71 * 
			    .1111111111111111 * (rhoa * 2. * t501 * sigmaaa + 
			    t502 * 2. * sigmabb)) - sigma * 1.333333333333333 
			    + sigmabb * 1.333333333333333 + sigmaaa * 
			    1.333333333333333);
		    t550 = t549 * .0052583256;
		    t553 = t218 * 8.8795591632e-4 * t36 * t270;
		    t555 = t491 * 3.7486538933976e-5 * t219;
		    t559 = 1 / t158 / t14;
/* Computing 2nd power */
		    d__1 = t174;
		    t560 = d__1 * d__1;
		    t565 = 1 / t4 / t163;
/* Computing 2nd power */
		    d__1 = t51;
		    t569 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t575 = d__1 * d__1;
		    t581 = 1 / t170 / t169;
		    t588 = t469 + t470 + t471 + t472 + t473 - t478 + t483 + 
			    t488 - t494 + t497 - t500 - t504 - t550 - t553 - 
			    t555 - t155 * .008064 * t175 - t8 * .006048 * 
			    t559 * t560 + t8 * .003024 * t159 * (t9 * .0784 * 
			    t565 * t11 + sigmaaa * .168 / t52 / t569 * t171 - 
			    t575 * .0448 / t4 / t569 / t163 * t581);
		    t589 = rho * t359;
		    t590 = t589 * t200;
		    t593 = 1 / t278 / t38;
		    t595 = 1 / t215;
		    t596 = t593 * t98 * t595;
		    t605 = 1 / t283 / t97;
/* Computing 2nd power */
		    d__1 = t294;
		    t607 = d__1 * d__1;
		    t611 = t593 * .1378556646443111 * t595;
		    t613 = t279 * .4135669939329333 * t501;
		    t615 = 1 / t290 / t38;
		    t616 = t615 * t595;
		    t618 = t291 * t501;
		    t620 = -t611 + t613 - t616 * 1.429973890497978 + t618 * 
			    3.431937337195148;
		    t625 = (t596 * -.1378556646443111 + t279 * 
			    .4135669939329333 * t284 * t211 * t294 + t280 * 
			    .4135669939329333 * t501 + t93 * 1.2407009817988 *
			     t605 * t607 - t285 * .6203504908994 * t620) * 
			    t298 * t97;
		    t628 = 1 / t93 / t38;
		    t631 = t297 * t628 * t97 * t211;
		    t633 = t299 * t294;
		    t638 = t595 * t307;
		    t639 = 1 / t302 / t104 * t593 * t638;
		    t642 = t303 * t615 * t638;
		    t645 = t304 * t501 * t307;
/* Computing 2nd power */
		    d__1 = t302;
		    t647 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t306;
		    t651 = d__1 * d__1;
		    t654 = 1 / t647 / t104 * t593 * t595 / t651;
		    t672 = (t596 * .03446391616107778 + t110 * 
			    .5250822119331622 * t284 * t292 * t294 - t311 * 
			    .2187842549721509 * t616 + t311 * 
			    .5250822119331622 * t618 + t111 * 2. * t605 * 
			    t607 - t314 * 1. * t620) * t318 * t97;
		    t679 = t317 / t111 / t110 * t97 * t291 * t211;
		    t681 = t319 * t294;
		    t684 = t593 * t117 * t595;
		    t693 = 1 / t325 / t116;
/* Computing 2nd power */
		    d__1 = t329;
		    t695 = d__1 * d__1;
		    t700 = -t611 + t613 - t616 * 2.201308720615045 + t618 * 
			    5.283140929476108;
		    t716 = t595 * t341;
/* Computing 2nd power */
		    d__1 = t336;
		    t725 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t340;
		    t729 = d__1 * d__1;
		    s1 = (t684 * -.1378556646443111 + t279 * 
			    .4135669939329333 * t326 * t211 * t329 + t322 * 
			    .4135669939329333 * t501 + t93 * 1.2407009817988 *
			     t693 * t695 - t327 * .6203504908994 * t700) * 
			    .02505897912236993 * t298 * t116 + t332 * 
			    .008352993040789975 * t628 * t116 * t211 + t333 * 
			    .02505897912236993 * t329 + .09995362477254242 / 
			    t336 / t122 * t593 * t716 + t337 * 
			    .1586317064594346 * t615 * t716 - t338 * 
			    .3807160955026431 * t501 * t341 - 
			    .1372209729363994 / t725 / t122 * t593 * t595 / 
			    t729 + (t684 * .03446391616107778 + t127 * 
			    .5250822119331622 * t326 * t292 * t329 - t345 * 
			    .2187842549721509 * t616 + t345 * 
			    .5250822119331622 * t618 + t128 * 2. * t693 * 
			    t695 - t348 * 1. * t700) * .002667310007273315 * 
			    t352 * t116 + t351 * 7.002785192652656e-4 / t128 /
			     t127 * t116 * t291 * t211 + t353 * 
			    .002667310007273315 * t329;
		    t771 = s1 - t625 * .05011795824473985 - t631 * 
			    .01670598608157995 - t633 * .05011795824473985 - 
			    t639 * .1270249377985834 - t642 * 
			    .2015953167456128 + t645 * .4838287601894708 + 
			    t654 * 2.560822746019441e-4 - t672 * 
			    .004431373767749538 - t679 * .001163417769936259 
			    - t681 * .004431373767749538;
/* Computing 2nd power */
		    d__1 = t140;
		    t776 = d__1 * d__1;
		    t777 = 1 / t776;
/* Computing 2nd power */
		    d__1 = t137;
		    t778 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t144;
		    t784 = d__1 * d__1;
		    t785 = 1 / t784;
		    t794 = t625 * .05011795824473985 + t631 * 
			    .01670598608157995 + t633 * .05011795824473985 + 
			    t639 * .1270249377985834 + t642 * 
			    .2015953167456128 - t645 * .4838287601894708 - 
			    t654 * 2.560822746019441e-4 + t672 * 
			    .004431373767749538 + t679 * .001163417769936259 
			    + t681 * .004431373767749538 + t771 * 
			    1.923661050931536 * t146 + t359 * 
			    3.847322101863073 * t368 + t135 * 
			    1.923661050931536 * (t777 * .4444444444444444 * 
			    t778 * t595 + t362 * 2.666666666666667 * t501 + 
			    t785 * .4444444444444444 * t778 * t595 - t365 * 
			    2.666666666666667 * t501);
		    t795 = rho * t794;
		    t796 = t795 * .19;
		    t814 = t777 * -.4444444444444444 * t137 * t501 - t140 * 
			    1.333333333333333 * t211 - t785 * 
			    .4444444444444444 * t137 * t501 + t144 * 
			    1.333333333333333 * t211;
		    t818 = rho * (t359 * 1.923661050931536 * t200 + t135 * 
			    1.923661050931536 * t814);
		    t823 = t43 * 1.423265147568e-4 * t513 * t490 * t89;
		    t835 = t485 * .00384780897072 * t43 * t219;
		    t837 = t218 * t36 * t191;
		    t840 = t222 * t217 * t191;
		    t843 = t44 * t227 * t191;
		    t846 = t251 * t38;
		    t849 = t71 * t211;
		    t856 = t44 * t49 * (rhob * t262 + t50 * (t846 * 
			    -.1111111111111111 * sigmaaa + t849 * 
			    .1111111111111111 * sigmaaa));
		    t859 = rhob * t49;
		    t861 = t513 * .0043129246896 * rhoa * t859;
		    t862 = t135 * t200;
		    t868 = t777 * .4444444444444444 * t211 + t785 * 
			    .4444444444444444 * t211;
		    t870 = t195 * .3654955996769919 * t868;
		    t872 = t204 * rhob * t207;
		    t874 = t178 * t211;
		    t876 = t195 * t814;
		    t878 = t590 * .3654955996769919 + t796 - t565 * .009408 * 
			    sigmaaa * t15 - .3308535951463467 / t52 + t818 * 
			    .19 - t823 - t44 * .0052583256 * t49 * (rhob * 2. 
			    * t187 + t182 * 162.0551065722879 * rhob - 
			    sigmabb * 2.) + t835 - t837 * 8.8795591632e-4 - 
			    t840 * .0012234370896 + t843 * .0385610544 - t856 
			    * .0105166512 - t861 + t862 * .7309911993539838 + 
			    t870 - t872 * .0370738512 + t874 * .3186864 + 
			    t876 * .3654955996769919;
		    v2rhoa2[i__] = t588 + t878;
		    t879 = t135 * t415;
		    t884 = 1 / t380 / t28;
/* Computing 2nd power */
		    d__1 = t396;
		    t885 = d__1 * d__1;
		    t890 = 1 / t18 / t385;
/* Computing 2nd power */
		    d__1 = t55;
		    t894 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t900 = d__1 * d__1;
		    t906 = 1 / t392 / t391;
		    t923 = t218 * t36 * t411;
		    t926 = t222 * t217 * t411;
		    t929 = t44 * t227 * t411;
		    t940 = t44 * t49 * (rhoa * t262 + t50 * (t846 * 
			    -.1111111111111111 * sigmabb + t849 * 
			    .1111111111111111 * sigmabb));
		    t942 = t469 + t470 + t471 + t472 + t473 + t879 * 
			    .7309911993539838 - t377 * .008064 * t397 - t22 * 
			    .006048 * t884 * t885 + t22 * .003024 * t381 * (
			    t23 * .0784 * t890 * t25 + sigmabb * .168 / t56 / 
			    t894 * t393 - t900 * .0448 / t18 / t894 / t385 * 
			    t906) - t44 * .0052583256 * t49 * (rhoa * 2. * 
			    t407 + rhoa * 162.0551065722879 * t403 - sigmaaa *
			     2.) - t923 * 8.8795591632e-4 - t926 * 
			    .0012234370896 + t929 * .0385610544 - t940 * 
			    .0105166512 - t478 + t483 + t488 - t494;
		    t950 = -t814;
		    t954 = rho * (t359 * 1.923661050931536 * t415 + t135 * 
			    1.923661050931536 * t950);
		    t956 = t205 * t207;
		    t958 = t37 * t211;
		    t960 = t195 * t950;
		    t962 = t589 * t415;
		    t964 = t497 - t500 - t504 - t550 - t553 - t555 - t890 * 
			    .009408 * sigmabb * t29 - .3308535951463467 / t56 
			    + t954 * .19 - t956 * .0370738512 + t958 * 
			    .3186864 + t960 * .3654955996769919 + t962 * 
			    .3654955996769919 + t796 - t823 + t835 - t861 + 
			    t870;
		    v2rhob2[i__] = t942 + t964;
		    t974 = -t478 + t483 + t488 - t494 + t497 - t500 - t503 * 
			    .3186864 - t549 * .0052583256 - t553 - t555 + 
			    t590 * .182747799838496;
		    t1003 = t840 * -6.117185448e-4 + t843 * .0192805272 - 
			    t856 * .0052583256 - t861 + t862 * 
			    .3654955996769919 - t872 * .0185369256 + t874 * 
			    .1593432 + t876 * .182747799838496 - t44 * 
			    .0052583256 * t49 * (t54 + t58 + t63 - t69 - t77 
			    + rhob * t407 + rhoa * t187) - t36 * .1593432 * 
			    t38 - t195 * .3654955996769919 * t868;
		    v2rhoab[i__] = t469 + t470 + t471 + t472 + t473 + t879 * 
			    .3654955996769919 - t923 * 4.4397795816e-4 - t926 
			    * 6.117185448e-4 + t929 * .0192805272 - t940 * 
			    .0052583256 + t974 + t954 * .095 - t956 * 
			    .0185369256 + t958 * .1593432 + t960 * 
			    .182747799838496 + t962 * .182747799838496 + t795 
			    * .19 + t818 * .095 - t823 + t835 - t837 * 
			    4.4397795816e-4 + t1003;
		    t1009 = t7 * t159;
		    t1033 = t50 * .1111111111111111 * t184;
		    t1040 = t218 * 4.4397795816e-4 * t36 * t437;
		    t1043 = t222 * 6.117185448e-4 * t217 * t437;
		    t1046 = t44 * .0192805272 * t227 * t437;
		    t1047 = t232 * .004690740740740741;
		    t1048 = t234 * .006462962962962963;
		    t1049 = t238 * .002255574074074074;
		    t1057 = rho * 1.333333333333333;
		    t1061 = t44 * .0052583256 * t49 * (t50 * (-t1047 - t1048 
			    + t1049 - t251 * .1111111111111111 * rhoa * t38 + 
			    t432 * .1111111111111111 * t211) + t1057);
		    t1063 = t44 * t859 * t62;
		    t1064 = t1063 * .0052583256;
		    t1066 = t218 * t36 * t443;
		    t1067 = t1066 * 4.4397795816e-4;
		    t1069 = t222 * t217 * t443;
		    t1070 = t1069 * 6.117185448e-4;
		    t1072 = t44 * t227 * t443;
		    t1073 = t1072 * .0192805272;
		    t1078 = t44 * t49 * (t50 * t240 - rho * 1.333333333333333)
			    ;
		    t1079 = t1078 * .0052583256;
		    v2rhoasigmaaa[i__] = t154 * .004032 * t15 - t155 * 
			    .004032 * t427 + t1009 * .003024 * t174 - t8 * 
			    .006048 * t559 * t174 * t426 + t8 * .003024 * 
			    t159 * (t420 * -.0168 * t154 * t11 - t165 * .0504 
			    * t171 + sigmaaa * .0168 / t4 / t569 / t51 * t581)
			     - t44 * .0052583256 * t49 * (rhob * t435 - t1033)
			     - t1040 - t1043 + t1046 - t1061 - t1064 - t1067 
			    - t1070 + t1073 - t1079;
		    t1081 = t1066 * 8.8795591632e-4;
		    t1082 = t1069 * .0012234370896;
		    t1083 = t1072 * .0385610544;
		    t1084 = t1078 * .0105166512;
		    v2rhoasigmaab[i__] = t1063 * -.0105166512 - t1081 - t1082 
			    + t1083 - t1084;
		    t1093 = t218 * 4.4397795816e-4 * t36 * t465;
		    t1096 = t222 * 6.117185448e-4 * t217 * t465;
		    t1099 = t44 * .0192805272 * t227 * t465;
		    t1110 = t44 * .0052583256 * t49 * (t50 * (-t1047 - t1048 
			    + t1049 - t251 * .1111111111111111 * rhob * t38 + 
			    t460 * .1111111111111111 * t211) + t1057);
		    v2rhoasigmabb[i__] = t44 * -.0052583256 * t49 * (rhob * 
			    t463 - rhoa * 2.) - t1093 - t1096 + t1099 - t1110 
			    - t1064 - t1067 - t1070 + t1073 - t1079;
		    t1119 = t44 * t49 * rhoa * t62;
		    t1120 = t1119 * .0052583256;
		    v2rhobsigmaaa[i__] = t44 * -.0052583256 * t49 * (rhoa * 
			    t435 - rhob * 2.) - t1040 - t1043 + t1046 - t1061 
			    - t1120 - t1067 - t1070 + t1073 - t1079;
		    v2rhobsigmaab[i__] = t1119 * -.0105166512 - t1081 - t1082 
			    + t1083 - t1084;
		    t1126 = t21 * t381;
		    v2rhobsigmabb[i__] = t376 * .004032 * t29 - t377 * 
			    .004032 * t457 + t1126 * .003024 * t396 - t22 * 
			    .006048 * t884 * t396 * t456 + t22 * .003024 * 
			    t381 * (t450 * -.0168 * t376 * t25 - t387 * .0504 
			    * t393 + sigmabb * .0168 / t18 / t894 / t55 * 
			    t906) - t44 * .0052583256 * t49 * (rhoa * t463 - 
			    t1033) - t1093 - t1096 + t1099 - t1110 - t1120 - 
			    t1067 - t1070 + t1073 - t1079;
/* Computing 2nd power */
		    d__1 = t426;
		    t1155 = d__1 * d__1;
		    v2sigmaaa2[i__] = t1009 * .006048 * t426 - t8 * .006048 * 
			    t559 * t1155 + t8 * .003024 * t159 * (-.0063 / t9 
			    / sigmaaa * t7 * t11 + .0063 / sigmaaa * t167 * 
			    t171 - .0063 / t4 / t569 / rhoa * t581);
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t456;
		    t1179 = d__1 * d__1;
		    v2sigmabb2[i__] = t1126 * .006048 * t456 - t22 * .006048 *
			     t884 * t1179 + t22 * .003024 * t381 * (-.0063 / 
			    t23 / sigmabb * t21 * t25 + .0063 / sigmabb * 
			    t389 * t393 - .0063 / t18 / t894 / rhob * t906);
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
} /* uks_xc_b3lyp__ */

/* Subroutine */ int rks_xc_b3lyp__(integer *ideriv, integer *npt, doublereal 
	*rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), exp(doublereal), atan(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, s3, t5, t6, t7, t8, s4, s5, s7, s8, s6, 
	    t10, t20, t30, t13, t23, t24, t34, t17, t26, t28, t56, t57, t62, 
	    t59, t65, t71, t75, t77, t14, t19, t25, t31, t36, t44, t47, t52, 
	    t61, t68, t74, t84, t88, t89, t94, t98, t99, t27, t32, t37, t42, 
	    t46, t85, t97, t200, t201, t112, t102, t123, t124, t103, t117, 
	    t109, t127, t129, t144, t146, t147, t175, t176, t178, t181, t182, 
	    t186, t187, t189, t190, t192, t111, t114, t119, t128, t130, t131, 
	    t134, t139, t149, t160, t164, t168, t177, t183, t188, t195, t196, 
	    t197, t198, t202, t204, t205, t207, t209, t212, t215, t216, t217, 
	    t218, t225, t231, t232, t240, t254, t258, t269, t271, t273, t277, 
	    t287, t289, t291, t292, rho, t311, t314, t342, t346, t363, t369, 
	    t377, t378, t392, t394, t395, t404, t406, t414, t415, t417, t419, 
	    t437, t446, t450, t485, t490, t554, sigma;


/*     P.J. Stephens, F.J. Devlin, C.F. Chabalowski, M.J. Frisch */
/*     Ab initio calculation of vibrational absorption and circular */
/*     dichroism spectra using density functional force fields */
/*     J. Phys. Chem. 98 (1994) 11623-11627 */


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
		t17 = 1 / t2;
		t20 = 1 / (t17 * .349 + 1.);
		t23 = t17 * .2533;
		t24 = exp(-t23);
/* Computing 2nd power */
		d__1 = rho;
		t26 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t28 = d__1 * d__1;
		t34 = t17 * t20;
		t56 = 1 / rho;
		t57 = pow_dd(&t56, &c_b2);
		t59 = pow_dd(&t56, &c_b4);
		t62 = 1 / (t57 * .6203504908994 + t59 * 10.29581201158544 + 
			42.7198);
		t65 = log(t57 * .6203504908994 * t62);
		t71 = atan(.0448998886412873 / (t59 * 1.575246635799487 + 
			13.072));
/* Computing 2nd power */
		d__1 = t59 * .7876233178997433 + .409286;
		t75 = d__1 * d__1;
		t77 = log(t75 * t62);
		zk[i__] = t3 * -.5908470131056179 - t5 * .003810001254882096 *
			 sigma / (t8 * .0317500104573508 * t10 + 1.) - t20 * 
			.0398358 * rho - t24 * .0052583256 * t20 / t28 / t26 /
			 rho * (t26 * .25 * (t28 * 11.48493600075277 * t26 + (
			2.611111111111111 - t17 * .09850555555555556 - t34 * 
			.1357222222222222) * sigma - (2.5 - t17 * 
			.01407222222222222 - t34 * .01938888888888889) * .5 * 
			sigma - (t23 + t34 * .349 - 11.) * .02777777777777778 
			* sigma) - t26 * .4583333333333333 * sigma) + rho * 
			.19 * (t65 * .0310907 + t71 * 20.5219729378375 + t77 *
			 .004431373767749538);
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
		t17 = 1 / t2;
		t19 = t17 * .349 + 1.;
		t20 = 1 / t19;
		t23 = t17 * .2533;
		t24 = exp(-t23);
		t25 = t24 * t20;
/* Computing 2nd power */
		d__1 = rho;
		t26 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t2;
		t28 = d__1 * d__1;
		t30 = 1 / t28 / t26 / rho;
		t31 = t28 * t26;
		t34 = t17 * t20;
		t36 = 2.611111111111111 - t17 * .09850555555555556 - t34 * 
			.1357222222222222;
		t44 = t23 + t34 * .349 - 11.;
		t47 = t31 * 11.48493600075277 + t36 * sigma - (2.5 - t17 * 
			.01407222222222222 - t34 * .01938888888888889) * .5 * 
			sigma - t44 * .02777777777777778 * sigma;
		t52 = t26 * .25 * t47 - t26 * .4583333333333333 * sigma;
		t56 = 1 / rho;
		t57 = pow_dd(&t56, &c_b2);
		t59 = pow_dd(&t56, &c_b4);
		t61 = t57 * .6203504908994 + t59 * 10.29581201158544 + 
			42.7198;
		t62 = 1 / t61;
		t65 = log(t57 * .6203504908994 * t62);
		t68 = t59 * 1.575246635799487 + 13.072;
		t71 = atan(.0448998886412873 / t68);
		t74 = t59 * .7876233178997433 + .409286;
/* Computing 2nd power */
		d__1 = t74;
		t75 = d__1 * d__1;
		t77 = log(t75 * t62);
		zk[i__] = t3 * -.5908470131056179 - t6 * .003810001254882096 *
			 t14 - t20 * .0398358 * rho - t25 * .0052583256 * t30 
			* t52 + rho * .19 * (t65 * .0310907 + t71 * 
			20.5219729378375 + t77 * .004431373767749538);
		t84 = 1 / t2 / t26;
/* Computing 2nd power */
		d__1 = t13;
		t88 = d__1 * d__1;
		t89 = 1 / t88;
		t94 = 1 / t31;
		t98 = sqrt(sigma * 1.587401051968199 * t94 + 1.);
		t99 = 1 / t98;
		t109 = t28 * rho;
		t112 = t44 * t56 * sigma;
		t117 = rho * sigma;
/* Computing 2nd power */
		d__1 = t19;
		t123 = d__1 * d__1;
		t124 = 1 / t123;
/* Computing 2nd power */
		d__1 = t26;
		t127 = d__1 * d__1;
		t129 = 1 / t127 / rho;
		t144 = t5 * t20;
		t146 = 1 / t109;
		t147 = t146 * t124;
/* Computing 2nd power */
		d__1 = t57;
		t175 = d__1 * d__1;
		t176 = 1 / t175;
		t178 = 1 / t26;
/* Computing 2nd power */
		d__1 = t61;
		t181 = d__1 * d__1;
		t182 = 1 / t181;
/* Computing 2nd power */
		d__1 = t59;
		t186 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t186;
		t187 = d__1 * d__1;
		t189 = 1 / t187 / t59;
		t190 = t189 * t178;
		t192 = t176 * -.2067834969664667 * t178 - t190 * 
			1.715968668597574;
/* Computing 2nd power */
		d__1 = t68;
		t200 = d__1 * d__1;
		t201 = 1 / t200;
		s1 = t2 * -.7877960174741572 + t84 * .005080001673176129 * 
			sigma * t14 + t6 * .001905000627441048 * t89 * (t7 * 
			-.08466669455293548 * t84 * t10 - sigma * 
			.106673350692263 * t30 * t99) - t20 * .0398358 - t25 *
			 .0052583256 * t30 * (rho * .5 * t47 + t26 * .25 * (
			t109 * 30.62649600200738 - t112 * .02777777777777778) 
			- t117 * .25) - t124 * .0046342314 * t17 - t129 * 
			4.4397795816e-4 * t24 * t20 * t52;
		s2 = s1 - t24 * 6.117185448e-4 * t124 * t129 * t52 + t25 * 
			.0192805272 / t28 / t127 * t52 - t25 * .0052583256 * 
			t30 * (t26 * .25 * ((t5 * .03283518518518519 + t144 * 
			.04524074074074074 - t147 * .01578901851851852) * 
			sigma - (t5 * .004690740740740741 + t144 * 
			.006462962962962963 - t147 * .002255574074074074) * 
			.5 * sigma - (t5 * -.08443333333333333 - t144 * 
			.1163333333333333 + t147 * .04060033333333333) * 
			.02777777777777778 * sigma + t112 * 
			.02777777777777778) - t117 * .6666666666666667);
		vrhoa[i__] = s2 + t65 * .005907233 + t71 * 3.899174858189126 
			+ t77 * 8.419610158724123e-4 + rho * .19 * ((t176 * 
			-.2067834969664667 * t62 * t178 - t57 * 
			.6203504908994 * t182 * t192) * .05011795824473985 / 
			t57 * t61 + t201 * .2419143800947354 * t189 * t178 / (
			t201 * .002016 + 1.) + (t74 * -.2625411059665811 * 
			t62 * t190 - t75 * 1. * t182 * t192) * 
			.004431373767749538 / t75 * t61);
		vsigmaaa[i__] = t5 * -.01524000501952839 * t14 + t6 * 
			.003810001254882096 * t89 * (.06350002091470161 / t7 *
			 t5 * t10 + t94 * .08000501301919725 * t99) + t25 * 
			5.842584e-4 * t146 - t25 * .0210333024 * t30 * (t26 * 
			.25 * t36 - t26 * .6666666666666667);
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
		t17 = 1 / t2;
		t19 = t17 * .349 + 1.;
		t20 = 1 / t19;
		t23 = t17 * .2533;
		t24 = exp(-t23);
		t25 = t24 * t20;
/* Computing 2nd power */
		d__1 = rho;
		t26 = d__1 * d__1;
		t27 = t26 * rho;
/* Computing 2nd power */
		d__1 = t2;
		t28 = d__1 * d__1;
		t30 = 1 / t28 / t27;
		t31 = t28 * t26;
		t32 = t31 * 11.48493600075277;
		t34 = t17 * t20;
		t36 = 2.611111111111111 - t17 * .09850555555555556 - t34 * 
			.1357222222222222;
		t37 = t36 * sigma;
		t42 = (2.5 - t17 * .01407222222222222 - t34 * 
			.01938888888888889) * .5 * sigma;
		t44 = t23 + t34 * .349 - 11.;
		t46 = t44 * .02777777777777778 * sigma;
		t47 = t32 + t37 - t42 - t46;
		t52 = t26 * .25 * t47 - t26 * .4583333333333333 * sigma;
		t56 = 1 / rho;
		t57 = pow_dd(&t56, &c_b2);
		t59 = pow_dd(&t56, &c_b4);
		t61 = t57 * .6203504908994 + t59 * 10.29581201158544 + 
			42.7198;
		t62 = 1 / t61;
		t65 = log(t57 * .6203504908994 * t62);
		t68 = t59 * 1.575246635799487 + 13.072;
		t71 = atan(.0448998886412873 / t68);
		t74 = t59 * .7876233178997433 + .409286;
/* Computing 2nd power */
		d__1 = t74;
		t75 = d__1 * d__1;
		t77 = log(t75 * t62);
		zk[i__] = t3 * -.5908470131056179 - t6 * .003810001254882096 *
			 t14 - t20 * .0398358 * rho - t25 * .0052583256 * t30 
			* t52 + rho * .19 * (t65 * .0310907 + t71 * 
			20.5219729378375 + t77 * .004431373767749538);
		t84 = 1 / t2 / t26;
		t85 = t84 * sigma;
/* Computing 2nd power */
		d__1 = t13;
		t88 = d__1 * d__1;
		t89 = 1 / t88;
		t94 = 1 / t31;
		t97 = sigma * 1.587401051968199 * t94 + 1.;
		t98 = sqrt(t97);
		t99 = 1 / t98;
		t102 = t7 * -.08466669455293548 * t84 * t10 - sigma * 
			.106673350692263 * t30 * t99;
		t103 = t89 * t102;
		t109 = t28 * rho;
		t111 = t44 * t56;
		t112 = t111 * sigma;
		t114 = t109 * 30.62649600200738 - t112 * .02777777777777778;
		t117 = rho * sigma;
		t119 = rho * .5 * t47 + t26 * .25 * t114 - t117 * .25;
/* Computing 2nd power */
		d__1 = t19;
		t123 = d__1 * d__1;
		t124 = 1 / t123;
/* Computing 2nd power */
		d__1 = t26;
		t127 = d__1 * d__1;
		t128 = t127 * rho;
		t129 = 1 / t128;
		t130 = t129 * t24;
		t131 = t20 * t52;
		t134 = t24 * t124;
		t139 = 1 / t28 / t127;
		t144 = t5 * t20;
		t146 = 1 / t109;
		t147 = t146 * t124;
		t149 = t5 * .03283518518518519 + t144 * .04524074074074074 - 
			t147 * .01578901851851852;
		t160 = t5 * -.08443333333333333 - t144 * .1163333333333333 + 
			t147 * .04060033333333333;
		t164 = t149 * sigma - (t5 * .004690740740740741 + t144 * 
			.006462962962962963 - t147 * .002255574074074074) * 
			.5 * sigma - t160 * .02777777777777778 * sigma + t112 
			* .02777777777777778;
		t168 = t26 * .25 * t164 - t117 * .6666666666666667;
/* Computing 2nd power */
		d__1 = t57;
		t175 = d__1 * d__1;
		t176 = 1 / t175;
		t177 = t176 * t62;
		t178 = 1 / t26;
/* Computing 2nd power */
		d__1 = t61;
		t181 = d__1 * d__1;
		t182 = 1 / t181;
		t183 = t57 * t182;
/* Computing 2nd power */
		d__1 = t59;
		t186 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t186;
		t187 = d__1 * d__1;
		t188 = t187 * t59;
		t189 = 1 / t188;
		t190 = t189 * t178;
		t192 = t176 * -.2067834969664667 * t178 - t190 * 
			1.715968668597574;
		t195 = t177 * -.2067834969664667 * t178 - t183 * 
			.6203504908994 * t192;
		t196 = 1 / t57;
		t197 = t195 * t196;
		t198 = t197 * t61;
/* Computing 2nd power */
		d__1 = t68;
		t200 = d__1 * d__1;
		t201 = 1 / t200;
		t202 = t201 * t189;
		t204 = t201 * .002016 + 1.;
		t205 = 1 / t204;
		t207 = t202 * t178 * t205;
		t209 = t74 * t62;
		t212 = t75 * t182;
		t215 = t209 * -.2625411059665811 * t190 - t212 * 1. * t192;
		t216 = 1 / t75;
		t217 = t215 * t216;
		t218 = t217 * t61;
		vrhoa[i__] = t2 * -.7877960174741572 + t85 * 
			.005080001673176129 * t14 + t6 * .001905000627441048 *
			 t103 - t20 * .0398358 - t25 * .0052583256 * t30 * 
			t119 - t124 * .0046342314 * t17 - t130 * 
			4.4397795816e-4 * t131 - t134 * 6.117185448e-4 * t129 
			* t52 + t25 * .0192805272 * t139 * t52 - t25 * 
			.0052583256 * t30 * t168 + t65 * .005907233 + t71 * 
			3.899174858189126 + t77 * 8.419610158724123e-4 + rho *
			 .19 * (t198 * .05011795824473985 + t207 * 
			.2419143800947354 + t218 * .004431373767749538);
		t225 = 1 / t7;
		t231 = t225 * .06350002091470161 * t5 * t10 + t94 * 
			.08000501301919725 * t99;
		t232 = t89 * t231;
		t240 = t26 * .25 * t36 - t26 * .6666666666666667;
		vsigmaaa[i__] = t5 * -.01524000501952839 * t14 + t6 * 
			.003810001254882096 * t232 + t25 * 5.842584e-4 * t146 
			- t25 * .0210333024 * t30 * t240;
		t254 = 1 / t2 / t27;
		t258 = rho * t114;
		t269 = 1 / t123 / t19;
		t271 = t127 * t26;
		t273 = 1 / t2 / t271;
		t277 = 1 / t271;
		t287 = t84 * t20;
		t289 = t94 * t124;
		t291 = 1 / t27;
		t292 = t291 * t269;
		t311 = t160 * t56 * sigma;
		t314 = t44 * t178 * sigma;
		s1 = t130 * -.00177591183264 * t20 * t119 - t134 * 
			.0024468741792 * t129 * t119 + t25 * .0771221088 * 
			t139 * t119 - t254 * .02370667447482193 * sigma * t14 
			- t25 * .0052583256 * t30 * (t258 + t31 * 
			25.52208000167282 - sigma * .5) + t25 * .0771221088 * 
			t139 * t168;
		s2 = s1 - t24 * 2.846530295136e-4 * t269 * t273 * t52 + t134 *
			 .0106031214432 * t277 * t52 - t25 * .1799515872 / 
			t28 / t128 * t52;
		t342 = s2 - t25 * .0105166512 * t30 * (t26 * .25 * ((t84 * 
			-.04378024691358025 - t287 * .06032098765432099 + 
			t289 * .03157803703703704 - t292 * 
			.003673578308641975) * sigma - (t84 * 
			-.006254320987654321 - t287 * .008617283950617284 + 
			t289 * .004511148148148148 - t292 * 
			5.247969012345679e-4) * .5 * sigma - (t84 * 
			.1125777777777778 + t287 * .1551111111111111 - t289 * 
			.08120066666666667 + t292 * .009446344222222222) * 
			.02777777777777778 * sigma + t311 * 
			.05555555555555556 - t314 * .05555555555555556) - 
			sigma * .6666666666666667) + t277 * .00769561794144 * 
			t24 * t131 - t25 * .0052583256 * t30 * (t32 + t37 - 
			t42 - t46 + t258) - t25 * .0210333024 * t30 * (rho * 
			.5 * t164 + t26 * .25 * (t311 * -.02777777777777778 + 
			t314 * .02777777777777778));
		t346 = t273 * t24;
/* Computing 2nd power */
		d__1 = sigma;
		t363 = d__1 * d__1;
		t369 = 1 / t98 / t97;
		t377 = 1 / t88 / t13;
/* Computing 2nd power */
		d__1 = t102;
		t378 = d__1 * d__1;
		t392 = 1 / t175 / t56;
		t394 = 1 / t127;
		t395 = t392 * t62 * t394;
		t404 = 1 / t181 / t61;
/* Computing 2nd power */
		d__1 = t192;
		t406 = d__1 * d__1;
		t414 = 1 / t188 / t56;
		t415 = t414 * t394;
		t417 = t189 * t291;
		t419 = t392 * -.1378556646443111 * t394 + t176 * 
			.4135669939329333 * t291 - t415 * 1.429973890497978 + 
			t417 * 3.431937337195148;
		t437 = t394 * t205;
/* Computing 2nd power */
		d__1 = t200;
		t446 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t204;
		t450 = d__1 * d__1;
		s1 = t130 * -.00177591183264 * t20 * t168 - t346 * 
			7.4973077867952e-5 * t131 - t346 * 2.0659774319712e-4 
			* t124 * t52 + t207 * .1838549288719989 + t198 * 
			.03808964826600229 + t218 * .003367844063489649 - t85 
			* .01016000334635226 * t103;
		s2 = s1 + t6 * .001905000627441048 * t89 * (t7 * 
			.3951112412470322 * t254 * t10 + sigma * 
			1.06673350692263 * t139 * t99 - t363 * 
			.4515557042823225 / t2 / t127 / t27 * t369) - t6 * 
			.003810001254882096 * t377 * t378 - t134 * 
			.0024468741792 * t129 * t168;
		s3 = s2 - .5251973449827715 / t28;
		s4 = s3 - t124 * .0061789752 * t5;
		s5 = s4;
		s7 = t269 * -.0021564623448 * t146;
		s8 = rho * .38 * ((t395 * -.1378556646443111 + t176 * 
			.4135669939329333 * t182 * t178 * t192 + t177 * 
			.4135669939329333 * t291 + t57 * 1.2407009817988 * 
			t404 * t406 - t183 * .6203504908994 * t419) * 
			.05011795824473985 * t196 * t61 + t195 * 
			.01670598608157995 / t57 / t56 * t61 * t178 + t197 * 
			.05011795824473985 * t192 + .1270249377985834 / t200 /
			 t68 * t392 * t437 + t201 * .2015953167456128 * t414 *
			 t437 - t202 * .4838287601894708 * t291 * t205 - 
			2.560822746019441e-4 / t446 / t68 * t392 * t394 / 
			t450 + (t395 * .03446391616107778 + t74 * 
			.5250822119331622 * t182 * t190 * t192 - t209 * 
			.2187842549721509 * t415 + t209 * .5250822119331622 * 
			t417 + t75 * 2. * t404 * t406 - t212 * 1. * t419) * 
			.004431373767749538 * t216 * t61 + t215 * 
			.001163417769936259 / t75 / t74 * t61 * t189 * t178 + 
			t217 * .004431373767749538 * t192);
		s6 = s7 + s8;
		t485 = s5 + s6;
		v2rhoa2[i__] = t342 + t485;
		t490 = t5 * t89;
		s1 = t84 * .02032000669270451 * t14 - t85 * 
			.005080001673176129 * t232 + t490 * 
			.007620002509764193 * t102 - t6 * .003810001254882096 
			* t377 * t102 * t231 + t6 * .001905000627441048 * t89 
			* (t225 * -.169333389105871 * t84 * t10 - t30 * 
			.640040104153578 * t99 + sigma * .3386667782117419 * 
			t273 * t369) - t25 * .0052583256 * t30 * (rho * 
			-.9444444444444444 - rho * .02777777777777778 * t44) 
			+ t291 * 4.933088424e-5 * t24 * t20;
		v2rhoasigmaaa[i__] = s1 + t134 * 6.79687272e-5 * t291 + t25 * 
			.0080822412 * t94 - t25 * .0105166512 * t30 * (t26 * 
			.25 * (t144 * -3e-23 + t111 * .05555555555555556) + 
			rho * 1.333333333333333) - t25 * .0105166512 * t94 * 
			t36 - t130 * .00177591183264 * t20 * t240 - t134 * 
			.0024468741792 * t129 * t240 + t25 * .0771221088 * 
			t139 * t240 - t25 * .0210333024 * t30 * (t26 * .25 * 
			t149 - rho * 1.333333333333333);
/* Computing 2nd power */
		d__1 = t231;
		t554 = d__1 * d__1;
		v2sigmaaa2[i__] = t490 * .03048001003905677 * t231 - t6 * 
			.007620002509764193 * t377 * t554 + t6 * 
			.003810001254882096 * t89 * (-.1270000418294032 / t7 /
			 sigma * t5 * t10 + .1600100260383945 / sigma * t94 * 
			t99 - .2540000836588064 / t2 / t128 * t369);
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
} /* rks_xc_b3lyp__ */


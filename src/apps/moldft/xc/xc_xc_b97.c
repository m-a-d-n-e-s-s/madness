/* xc_xc_b97.f -- translated by f2c (version 20050501).
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

/* :XC_B97subrstart */
/*    Generated: Tue Jul 22 14:39:42 GMT 2003 */
/* Subroutine */ int uks_xc_b97__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal s1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t20, t12, t21, 
	    t14, t15, t32, t35, t27, t19, t28, t37, t43, t45, t49, t16, t17, 
	    t22, t29, t30, t33, t34, t39, t47, t51, t59, t61, t62, t65, t67, 
	    t71, t72, t76, t77, t84, t85, t88, t89, t92, t94, t11, t18, t24, 
	    t31, t42, t46, t50, t100, t200, t102, t130, t122, t114, t106, 
	    t115, t116, t124, t119, t132, t143, t146, t147, t148, t149, t152, 
	    t153, t155, t156, t157, t158, t159, t161, t177, t193, t196, t201, 
	    t53, t54, t68, t75, t86, t90, t93, t97, t128, t13, t23, t26, t41, 
	    t44, t48, t52, t55, t56, t60, t64, t78, t81, t87, t96, t99, t103, 
	    t107, t110, t111, t118, t126, t129, t134, t139, t142, t144, t160, 
	    t163, t166, rho, t168, t173, t176, t180, t181, t183, t190, t197, 
	    t202, t205, t212, t218, t222, t225, t227, t233, t236, t237, t240, 
	    t241, t243, t244, t247, t258, t259, t269, t277, t278, t279, t280, 
	    t281, t284, t285, t288, t292, t295, t300, t301, t305, t306, t319, 
	    t320, t321, t325, t329, t334, t335, t338, t339, t346, t363, t367, 
	    t370, t384, t389, t399, t405, t409, t412, t414, t420, t423, t424, 
	    t427, t428, t430, t431, t434, t445, t446, t456, t469, t502, t506, 
	    t532, t536, t80, t83, t91, t95, t105, t108, t109, t123, t131, 
	    t164, t169, t174, t175, t186, t188, t207, t251, rhoa, rhob, t252, 
	    t261, t266, t209, t211, t217, t221, t224, t230, t242, t245, t246, 
	    t255, t260, t272, t273, t276, t282, t283, t286, t287, t291, t294, 
	    t297, t298, t302, t303, t304, t307, t308, t313, t314, t315, t316, 
	    t317, t318, t322, t326, t331, t332, t333, t337, t340, t342, t343, 
	    t347, t348, t353, t354, t360, t361, t362, t364, t365, t366, t369, 
	    t377, t380, t388, t392, t396, t398, t404, t408, t411, t417, t429, 
	    t432, t433, t439, t442, t447, t448, t459, t460, t463, t466, t471, 
	    t472, t473, t475, t476, t477, t478, t479, t487, t490, t498, t505, 
	    t509, t518, t519, t528, t535, t539, t548, t549, t558, sigma, t562,
	     t566, t567, t574, t576, t577, t582, t603, t609, t610, t612, t613,
	     t617, t618, t621, t623, t624, t625, t628, t629, t630, t632, t633,
	     t634, t635, t638, t639, t641, t643, t646, t647, t650, t651, t652,
	     t653, t654, t655, t658, t659, t664, t665, t666, t669, t671, t672,
	     t675, t676, t677, t683, t689, t692, t694, t699, t703, t706, t709,
	     t711, t717, t719, t723, t725, t730, t733, t735, t737, t742, t744,
	     t746, t748, t749, t761, t777, t780, t792, t804, t807, t810, t812,
	     t817, t819, t820, t824, t829, t830, t831, t834, t835, t838, t842,
	     t848, t851, t856, t858, t876, t878, t881, t885, t886, t887, t899,
	     t900, t937, t943, t962, t964, t967, t971, t972, t973, t987, t989,
	     t990, t995, t1000, t1021, t1027, t1029, t1034, t1035, t1037, 
	    t1038, t1068, t1070, t1073, t1078, t1081, t1084, t1090, t1099, 
	    t1105, t1109, t1112, t1116, t1118, t1121, t1129, t1130, t1160, 
	    t1170, t1182, t1185, t1195, t1207, t1210, t1216, t1260, t1268, 
	    t1271, t1282, t1285, t1291, t1333, t1338, t1364, t1372, t1377, 
	    sigmaaa, sigmaab, sigmabb;


/*     A.D. Becke */
/*     Density-functional thermochemistry. V. Systematic optimization of */
/*     exchange-correlation functionals */
/*     J. Chem. Phys. 107 (1997) 8554-8560 */


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
		    t8 = sigmabb / t5 / t4;
		    t10 = t8 * .004 + 1.;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t19 = t14 / t2 / t15 / rhob;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t27 = 1 / rhob;
		    t28 = pow_dd(&t27, &c_b2);
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t43 = log(32.1646831778707 / (t32 * 11.12037486309468 + 
			    t28 * 3.844746237447211 + t35 * 1.644733775567609 
			    + t37 * .2405871291288192) + 1.);
		    t45 = t8 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhob * (t8 * .0020292 / 
			    t10 + .8094 + t19 * 1.19696e-5 / t20) - rhob * 
			    .03109 * (t28 * .1274696188700087 + 1.) * t43 * (
			    t8 * .46974 / t45 + .1737 - t19 * .099472 / t49);
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
		    t8 = sigmaaa / t5 / t4;
		    t10 = t8 * .004 + 1.;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t19 = t14 / t2 / t15 / rhoa;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t27 = 1 / rhoa;
		    t28 = pow_dd(&t27, &c_b2);
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t43 = log(32.1646831778707 / (t32 * 11.12037486309468 + 
			    t28 * 3.844746237447211 + t35 * 1.644733775567609 
			    + t37 * .2405871291288192) + 1.);
		    t45 = t8 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhoa * (t8 * .0020292 / 
			    t10 + .8094 + t19 * 1.19696e-5 / t20) - rhoa * 
			    .03109 * (t28 * .1274696188700087 + 1.) * t43 * (
			    t8 * .46974 / t45 + .1737 - t19 * .099472 / t49);
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
		    t10 = sigmaaa / t7 / t6;
		    t12 = t10 * .004 + 1.;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t16 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t17 = d__1 * d__1;
		    t21 = t16 / t4 / t17 / rhoa;
/* Computing 2nd power */
		    d__1 = t12;
		    t22 = d__1 * d__1;
		    t29 = 1 / rhoa;
		    t30 = pow_dd(&t29, &c_b2);
		    t33 = rhoa * (t30 * .1274696188700087 + 1.);
		    t34 = pow_dd(&t29, &c_b4);
		    t37 = sqrt(t29);
/* Computing 2nd power */
		    d__1 = t30;
		    t39 = d__1 * d__1;
		    t45 = log(32.1646831778707 / (t34 * 11.12037486309468 + 
			    t30 * 3.844746237447211 + t37 * 1.644733775567609 
			    + t39 * .2405871291288192) + 1.);
		    t47 = t10 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t47;
		    t51 = d__1 * d__1;
		    t59 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = rhob;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t59;
		    t62 = d__1 * d__1;
		    t65 = sigmabb / t62 / t61;
		    t67 = t65 * .004 + 1.;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t71 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t61;
		    t72 = d__1 * d__1;
		    t76 = t71 / t59 / t72 / rhob;
/* Computing 2nd power */
		    d__1 = t67;
		    t77 = d__1 * d__1;
		    t84 = 1 / rhob;
		    t85 = pow_dd(&t84, &c_b2);
		    t88 = rhob * (t85 * .1274696188700087 + 1.);
		    t89 = pow_dd(&t84, &c_b4);
		    t92 = sqrt(t84);
/* Computing 2nd power */
		    d__1 = t85;
		    t94 = d__1 * d__1;
		    t100 = log(32.1646831778707 / (t89 * 11.12037486309468 + 
			    t85 * 3.844746237447211 + t92 * 1.644733775567609 
			    + t94 * .2405871291288192) + 1.);
		    t102 = t65 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t102;
		    t106 = d__1 * d__1;
		    t114 = rhoa + rhob;
		    t115 = 1 / t114;
		    t116 = pow_dd(&t115, &c_b2);
		    t119 = pow_dd(&t115, &c_b4);
		    t122 = sqrt(t115);
/* Computing 2nd power */
		    d__1 = t116;
		    t124 = d__1 * d__1;
		    t130 = log(16.0818243221511 / (t119 * 5.98255043577108 + 
			    t116 * 2.225569421150687 + t122 * 
			    .8004286349993634 + t124 * .1897004325747559) + 
			    1.);
		    t132 = (t116 * .1325688999052018 + 1.) * .062182 * t130;
		    t143 = log(29.60857464321668 / (t119 * 8.157414703487641 
			    + t116 * 2.247591863577616 + t122 * 
			    .4300972471276643 + t124 * .1911512595127338) + 
			    1.);
		    t146 = rhoa - rhob * 1.;
		    t147 = t146 * t115;
		    t148 = t147 + 1.;
		    t149 = pow_dd(&t148, &c_b2);
		    t152 = 1. - t147 * 1.;
		    t153 = pow_dd(&t152, &c_b2);
		    t155 = t149 * t148 + t153 * t152 - 2.;
/* Computing 2nd power */
		    d__1 = t146;
		    t156 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t156;
		    t157 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t114;
		    t158 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t158;
		    t159 = d__1 * d__1;
		    t161 = t157 / t159;
		    t177 = log(32.1646831778707 / (t119 * 11.12037486309468 + 
			    t116 * 3.844746237447211 + t122 * 
			    1.644733775567609 + t124 * .2405871291288192) + 
			    1.);
		    t193 = t10 * .5 + t65 * .5;
		    t196 = t10 * .003 + 1. + t65 * .003;
/* Computing 2nd power */
		    d__1 = t193;
		    t200 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t196;
		    t201 = d__1 * d__1;
		    zk[i__] = t4 * -.9305257363491 * rhoa * (t10 * .0020292 / 
			    t12 + .8094 + t21 * 1.19696e-5 / t22) - t33 * 
			    .03109 * t45 * (t10 * .46974 / t47 + .1737 - t21 *
			     .099472 / t51) - t59 * .9305257363491 * rhob * (
			    t65 * .0020292 / t67 + .8094 + t76 * 1.19696e-5 / 
			    t77) - t88 * .03109 * t100 * (t65 * .46974 / t102 
			    + .1737 - t76 * .099472 / t106) + (t114 * (-t132 
			    + (t116 * .06901399211255825 + 1.) * 
			    .03799574853701528 * t143 * t155 * (1. - t161 * 
			    1.) + ((t116 * .1274696188700087 + 1.) * -.03109 *
			     t177 + t132) * 1.923661050931536 * t155 * t161) 
			    + t33 * .03109 * t45 + t88 * .03109 * t100) * (
			    t193 * .0044826 / t196 + .9454 - t200 * 
			    1.654596e-4 / t201);
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
		    t7 = 1 / t5 / t4;
		    t8 = sigmabb * t7;
		    t10 = t8 * .004 + 1.;
		    t11 = 1 / t10;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t18 = 1 / t2 / t15 / rhob;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t8 * .0020292 * t11 + .8094 + t19 * 1.19696e-5 * 
			    t21;
		    t27 = 1 / rhob;
		    t28 = pow_dd(&t27, &c_b2);
		    t30 = t28 * .1274696188700087 + 1.;
		    t31 = rhob * t30;
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t39 = t32 * 11.12037486309468 + t28 * 3.844746237447211 + 
			    t35 * 1.644733775567609 + t37 * .2405871291288192;
		    t42 = 32.1646831778707 / t39 + 1.;
		    t43 = log(t42);
		    t45 = t8 * .2 + 1.;
		    t46 = 1 / t45;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t8 * .46974 * t46 + .1737 - t19 * .099472 * t50;
		    t54 = t43 * t53;
		    zk[i__] = t3 * -.9305257363491 * t24 - t31 * .03109 * t54;
		    vrhoa[i__] = 0.;
		    t62 = sigmabb / t5 / t4 / rhob;
		    t68 = t14 / t2 / t15 / t4;
/* Computing 2nd power */
		    d__1 = t15;
		    t72 = d__1 * d__1;
		    t75 = t14 * sigmabb / t72 / rhob;
		    t77 = 1 / t20 / t10;
		    t86 = 1 / t37;
/* Computing 2nd power */
		    d__1 = t39;
		    t90 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t93 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t93;
		    t94 = d__1 * d__1;
		    t97 = 1 / t4;
		    t119 = 1 / t49 / t45;
		    vrhob[i__] = t2 * -1.2407009817988 * t24 - t3 * 
			    .9305257363491 * (t62 * -.0054112 * t11 - t68 * 
			    4.219306666666667e-5 * t21 + t75 * 
			    2.553514666666667e-7 * t77) - t30 * .03109 * t43 *
			     t53 + t27 * .001321010150222857 * t86 * t54 + 
			    t31 * 1. / t90 * (-1.853395810515781 / t94 / t32 *
			     t97 - t86 * 1.28158207914907 * t97 - 
			    .8223668877838045 / t35 * t97 - .1603914194192128 
			    / t28 * t97) / t42 * t53 - t31 * .03109 * t43 * (
			    t62 * -1.25264 * t46 + t68 * .7810453333333333 * 
			    t50 - t75 * .1061034666666667 * t119);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t128 = sigmabb * t18;
		    t132 = t14 / t72;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * .0020292 * 
			    t11 + t128 * 1.58224e-5 * t21 - t132 * 9.57568e-8 
			    * t77) - t31 * .03109 * t43 * (t7 * .46974 * t46 
			    - t128 * .292892 * t50 + t132 * .0397888 * t119);
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
		    t7 = 1 / t5 / t4;
		    t8 = sigmaaa * t7;
		    t10 = t8 * .004 + 1.;
		    t11 = 1 / t10;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t18 = 1 / t2 / t15 / rhoa;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t8 * .0020292 * t11 + .8094 + t19 * 1.19696e-5 * 
			    t21;
		    t27 = 1 / rhoa;
		    t28 = pow_dd(&t27, &c_b2);
		    t30 = t28 * .1274696188700087 + 1.;
		    t31 = rhoa * t30;
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t39 = t32 * 11.12037486309468 + t28 * 3.844746237447211 + 
			    t35 * 1.644733775567609 + t37 * .2405871291288192;
		    t42 = 32.1646831778707 / t39 + 1.;
		    t43 = log(t42);
		    t45 = t8 * .2 + 1.;
		    t46 = 1 / t45;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t8 * .46974 * t46 + .1737 - t19 * .099472 * t50;
		    t54 = t43 * t53;
		    zk[i__] = t3 * -.9305257363491 * t24 - t31 * .03109 * t54;
		    t62 = sigmaaa / t5 / t4 / rhoa;
		    t68 = t14 / t2 / t15 / t4;
/* Computing 2nd power */
		    d__1 = t15;
		    t72 = d__1 * d__1;
		    t75 = t14 * sigmaaa / t72 / rhoa;
		    t77 = 1 / t20 / t10;
		    t86 = 1 / t37;
/* Computing 2nd power */
		    d__1 = t39;
		    t90 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t32;
		    t93 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t93;
		    t94 = d__1 * d__1;
		    t97 = 1 / t4;
		    t119 = 1 / t49 / t45;
		    vrhoa[i__] = t2 * -1.2407009817988 * t24 - t3 * 
			    .9305257363491 * (t62 * -.0054112 * t11 - t68 * 
			    4.219306666666667e-5 * t21 + t75 * 
			    2.553514666666667e-7 * t77) - t30 * .03109 * t43 *
			     t53 + t27 * .001321010150222857 * t86 * t54 + 
			    t31 * 1. / t90 * (-1.853395810515781 / t94 / t32 *
			     t97 - t86 * 1.28158207914907 * t97 - 
			    .8223668877838045 / t35 * t97 - .1603914194192128 
			    / t28 * t97) / t42 * t53 - t31 * .03109 * t43 * (
			    t62 * -1.25264 * t46 + t68 * .7810453333333333 * 
			    t50 - t75 * .1061034666666667 * t119);
		    vrhob[i__] = 0.;
		    t128 = sigmaaa * t18;
		    t132 = t14 / t72;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * .0020292 * 
			    t11 + t128 * 1.58224e-5 * t21 - t132 * 9.57568e-8 
			    * t77) - t31 * .03109 * t43 * (t7 * .46974 * t46 
			    - t128 * .292892 * t50 + t132 * .0397888 * t119);
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
		    t9 = 1 / t7 / t6;
		    t10 = sigmaaa * t9;
		    t12 = t10 * .004 + 1.;
		    t13 = 1 / t12;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t16 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t17 = d__1 * d__1;
		    t20 = 1 / t4 / t17 / rhoa;
		    t21 = t16 * t20;
/* Computing 2nd power */
		    d__1 = t12;
		    t22 = d__1 * d__1;
		    t23 = 1 / t22;
		    t26 = t10 * .0020292 * t13 + .8094 + t21 * 1.19696e-5 * 
			    t23;
		    t29 = 1 / rhoa;
		    t30 = pow_dd(&t29, &c_b2);
		    t32 = t30 * .1274696188700087 + 1.;
		    t33 = rhoa * t32;
		    t34 = pow_dd(&t29, &c_b4);
		    t37 = sqrt(t29);
/* Computing 2nd power */
		    d__1 = t30;
		    t39 = d__1 * d__1;
		    t41 = t34 * 11.12037486309468 + t30 * 3.844746237447211 + 
			    t37 * 1.644733775567609 + t39 * .2405871291288192;
		    t44 = 32.1646831778707 / t41 + 1.;
		    t45 = log(t44);
		    t47 = t10 * .2 + 1.;
		    t48 = 1 / t47;
/* Computing 2nd power */
		    d__1 = t47;
		    t51 = d__1 * d__1;
		    t52 = 1 / t51;
		    t55 = t10 * .46974 * t48 + .1737 - t21 * .099472 * t52;
		    t56 = t45 * t55;
		    t59 = pow_dd(&rhob, &c_b2);
		    t60 = t59 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t59;
		    t62 = d__1 * d__1;
		    t64 = 1 / t62 / t61;
		    t65 = sigmabb * t64;
		    t67 = t65 * .004 + 1.;
		    t68 = 1 / t67;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t71 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t61;
		    t72 = d__1 * d__1;
		    t75 = 1 / t59 / t72 / rhob;
		    t76 = t71 * t75;
/* Computing 2nd power */
		    d__1 = t67;
		    t77 = d__1 * d__1;
		    t78 = 1 / t77;
		    t81 = t65 * .0020292 * t68 + .8094 + t76 * 1.19696e-5 * 
			    t78;
		    t84 = 1 / rhob;
		    t85 = pow_dd(&t84, &c_b2);
		    t87 = t85 * .1274696188700087 + 1.;
		    t88 = rhob * t87;
		    t89 = pow_dd(&t84, &c_b4);
		    t92 = sqrt(t84);
/* Computing 2nd power */
		    d__1 = t85;
		    t94 = d__1 * d__1;
		    t96 = t89 * 11.12037486309468 + t85 * 3.844746237447211 + 
			    t92 * 1.644733775567609 + t94 * .2405871291288192;
		    t99 = 32.1646831778707 / t96 + 1.;
		    t100 = log(t99);
		    t102 = t65 * .2 + 1.;
		    t103 = 1 / t102;
/* Computing 2nd power */
		    d__1 = t102;
		    t106 = d__1 * d__1;
		    t107 = 1 / t106;
		    t110 = t65 * .46974 * t103 + .1737 - t76 * .099472 * t107;
		    t111 = t100 * t110;
		    t114 = rhoa + rhob;
		    t115 = 1 / t114;
		    t116 = pow_dd(&t115, &c_b2);
		    t118 = t116 * .1325688999052018 + 1.;
		    t119 = pow_dd(&t115, &c_b4);
		    t122 = sqrt(t115);
/* Computing 2nd power */
		    d__1 = t116;
		    t124 = d__1 * d__1;
		    t126 = t119 * 5.98255043577108 + t116 * 2.225569421150687 
			    + t122 * .8004286349993634 + t124 * 
			    .1897004325747559;
		    t129 = 16.0818243221511 / t126 + 1.;
		    t130 = log(t129);
		    t132 = t118 * .062182 * t130;
		    t134 = t116 * .06901399211255825 + 1.;
		    t139 = t119 * 8.157414703487641 + t116 * 
			    2.247591863577616 + t122 * .4300972471276643 + 
			    t124 * .1911512595127338;
		    t142 = 29.60857464321668 / t139 + 1.;
		    t143 = log(t142);
		    t144 = t134 * t143;
		    t146 = rhoa - rhob * 1.;
		    t147 = t146 * t115;
		    t148 = t147 + 1.;
		    t149 = pow_dd(&t148, &c_b2);
		    t152 = 1. - t147 * 1.;
		    t153 = pow_dd(&t152, &c_b2);
		    t155 = t149 * t148 + t153 * t152 - 2.;
/* Computing 2nd power */
		    d__1 = t146;
		    t156 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t156;
		    t157 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t114;
		    t158 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t158;
		    t159 = d__1 * d__1;
		    t160 = 1 / t159;
		    t161 = t157 * t160;
		    t163 = 1. - t161 * 1.;
		    t166 = t144 * .03799574853701528 * t155 * t163;
		    t168 = t116 * .1274696188700087 + 1.;
		    t173 = t119 * 11.12037486309468 + t116 * 
			    3.844746237447211 + t122 * 1.644733775567609 + 
			    t124 * .2405871291288192;
		    t176 = 32.1646831778707 / t173 + 1.;
		    t177 = log(t176);
		    t180 = t168 * -.03109 * t177 + t132;
		    t181 = t180 * t155;
		    t183 = t181 * 1.923661050931536 * t161;
		    t190 = t114 * (-t132 + t166 + t183) + t33 * .03109 * t45 
			    + t88 * .03109 * t100;
		    t193 = t10 * .5 + t65 * .5;
		    t196 = t10 * .003 + 1. + t65 * .003;
		    t197 = 1 / t196;
/* Computing 2nd power */
		    d__1 = t193;
		    t200 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t196;
		    t201 = d__1 * d__1;
		    t202 = 1 / t201;
		    t205 = t193 * .0044826 * t197 + .9454 - t200 * 
			    1.654596e-4 * t202;
		    zk[i__] = t5 * -.9305257363491 * t26 - t33 * .03109 * t56 
			    - t60 * .9305257363491 * t81 - t88 * .03109 * 
			    t111 + t190 * t205;
		    t212 = sigmaaa / t7 / t6 / rhoa;
		    t218 = t16 / t4 / t17 / t6;
/* Computing 2nd power */
		    d__1 = t17;
		    t222 = d__1 * d__1;
		    t225 = t16 * sigmaaa / t222 / rhoa;
		    t227 = 1 / t22 / t12;
		    t233 = t32 * t45;
		    t236 = 1 / t39;
		    t237 = t29 * t236;
/* Computing 2nd power */
		    d__1 = t41;
		    t240 = d__1 * d__1;
		    t241 = 1 / t240;
/* Computing 2nd power */
		    d__1 = t34;
		    t243 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t243;
		    t244 = d__1 * d__1;
		    t247 = 1 / t6;
		    t258 = -1.853395810515781 / t244 / t34 * t247 - t236 * 
			    1.28158207914907 * t247 - .8223668877838045 / t37 
			    * t247 - .1603914194192128 / t30 * t247;
		    t259 = 1 / t44;
		    t269 = 1 / t51 / t47;
		    t277 = 1 / t158;
		    t278 = 1 / t124 * t277;
		    t279 = t278 * t130;
		    t280 = t279 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t126;
		    t281 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t119;
		    t284 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t284;
		    t285 = d__1 * d__1;
		    t288 = 1 / t285 / t119 * t277;
		    t292 = 1 / t122 * t277;
		    t295 = 1 / t116 * t277;
		    t300 = t118 / t281 * (t288 * -.99709173929518 - t278 * 
			    .7418564737168958 - t292 * .4002143174996817 - 
			    t295 * .1264669550498372) / t129;
		    t301 = t300 * 1.;
		    t305 = t278 * 8.740794299481065e-4 * t143 * t155 * t163;
/* Computing 2nd power */
		    d__1 = t139;
		    t306 = d__1 * d__1;
		    t319 = t134 * 1.124999956683108 / t306 * (t288 * 
			    -1.35956911724794 - t278 * .7491972878592054 - 
			    t292 * .2150486235638321 - t295 * 
			    .1274341730084892) / t142 * t155 * t163;
		    t320 = t146 * t277;
		    t321 = t320 * 1.;
		    t325 = t115 * 1.;
		    t329 = t149 * 1.333333333333333 * (t115 - t321) + t153 * 
			    1.333333333333333 * (-t325 + t320);
		    t334 = t156 * t146 * t160;
		    t335 = t334 * 4.;
		    t338 = t157 / t159 / t114;
		    t339 = t338 * 4.;
/* Computing 2nd power */
		    d__1 = t173;
		    t346 = d__1 * d__1;
		    t363 = (t278 * .001321010150222857 * t177 + t168 * 1. / 
			    t346 * (t288 * -1.853395810515781 - t278 * 
			    1.28158207914907 - t292 * .8223668877838045 - 
			    t295 * .1603914194192128) / t176 - t279 * 
			    .002747799777968419 - t300 * 1.) * 
			    1.923661050931536 * t155 * t161;
		    t367 = t181 * t334;
		    t370 = t181 * 7.694644203726145 * t338;
		    t384 = t193 * t202;
		    t389 = t200 / t201 / t196;
		    vrhoa[i__] = t4 * -1.2407009817988 * t26 - t5 * 
			    .9305257363491 * (t212 * -.0054112 * t13 - t218 * 
			    4.219306666666667e-5 * t23 + t225 * 
			    2.553514666666667e-7 * t227) - t233 * .03109 * 
			    t55 + t237 * .001321010150222857 * t56 + t33 * 1. 
			    * t241 * t258 * t259 * t55 - t33 * .03109 * t45 * 
			    (t212 * -1.25264 * t48 + t218 * .7810453333333333 
			    * t52 - t225 * .1061034666666667 * t269) + (-t132 
			    + t166 + t183 + t114 * (t280 + t301 - t305 - t319 
			    + t144 * .03799574853701528 * t329 * t163 + t144 *
			     .03799574853701528 * t155 * (-t335 + t339) + 
			    t363 + t180 * 1.923661050931536 * t329 * t161 + 
			    t367 * 7.694644203726145 - t370) + t233 * .03109 
			    - t237 * .001321010150222857 * t45 - t33 * 1. * 
			    t241 * t258 * t259) * t205 + t190 * (t212 * 
			    -.0059768 * t197 + t384 * 4.770864e-4 * t212 - 
			    t389 * 2.6473536e-6 * t212);
		    t399 = sigmabb / t62 / t61 / rhob;
		    t405 = t71 / t59 / t72 / t61;
/* Computing 2nd power */
		    d__1 = t72;
		    t409 = d__1 * d__1;
		    t412 = t71 * sigmabb / t409 / rhob;
		    t414 = 1 / t77 / t67;
		    t420 = t87 * t100;
		    t423 = 1 / t94;
		    t424 = t84 * t423;
/* Computing 2nd power */
		    d__1 = t96;
		    t427 = d__1 * d__1;
		    t428 = 1 / t427;
/* Computing 2nd power */
		    d__1 = t89;
		    t430 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t430;
		    t431 = d__1 * d__1;
		    t434 = 1 / t61;
		    t445 = -1.853395810515781 / t431 / t89 * t434 - t423 * 
			    1.28158207914907 * t434 - .8223668877838045 / t92 
			    * t434 - .1603914194192128 / t85 * t434;
		    t446 = 1 / t99;
		    t456 = 1 / t106 / t102;
		    t469 = t149 * 1.333333333333333 * (-t325 - t321) + t153 * 
			    1.333333333333333 * (t115 + t320);
		    s1 = t59 * -1.2407009817988 * t81 - t60 * .9305257363491 *
			     (t399 * -.0054112 * t68 - t405 * 
			    4.219306666666667e-5 * t78 + t412 * 
			    2.553514666666667e-7 * t414) - t420 * .03109 * 
			    t110 + t424 * .001321010150222857 * t111;
		    vrhob[i__] = s1 + t88 * 1. * t428 * t445 * t446 * t110 - 
			    t88 * .03109 * t100 * (t399 * -1.25264 * t103 + 
			    t405 * .7810453333333333 * t107 - t412 * 
			    .1061034666666667 * t456) + (-t132 + t166 + t183 
			    + t114 * (t280 + t301 - t305 - t319 + t144 * 
			    .03799574853701528 * t469 * t163 + t144 * 
			    .03799574853701528 * t155 * (t335 + t339) + t363 
			    + t180 * 1.923661050931536 * t469 * t161 - t367 * 
			    7.694644203726145 - t370) + t420 * .03109 - t424 *
			     .001321010150222857 * t100 - t88 * 1. * t428 * 
			    t445 * t446) * t205 + t190 * (t399 * -.0059768 * 
			    t197 + t384 * 4.770864e-4 * t399 - t389 * 
			    2.6473536e-6 * t399);
		    t502 = sigmaaa * t20;
		    t506 = t16 / t222;
		    vsigmaaa[i__] = t5 * -.9305257363491 * (t9 * .0020292 * 
			    t13 + t502 * 1.58224e-5 * t23 - t506 * 9.57568e-8 
			    * t227) - t33 * .03109 * t45 * (t9 * .46974 * t48 
			    - t502 * .292892 * t52 + t506 * .0397888 * t269) 
			    + t190 * (t9 * .0022413 * t197 - t384 * 
			    1.789074e-4 * t9 + t389 * 9.927576e-7 * t9);
		    vsigmaab[i__] = 0.;
		    t532 = sigmabb * t75;
		    t536 = t71 / t409;
		    vsigmabb[i__] = t60 * -.9305257363491 * (t64 * .0020292 * 
			    t68 + t532 * 1.58224e-5 * t78 - t536 * 9.57568e-8 
			    * t414) - t88 * .03109 * t100 * (t64 * .46974 * 
			    t103 - t532 * .292892 * t107 + t536 * .0397888 * 
			    t456) + t190 * (t64 * .0022413 * t197 - t384 * 
			    1.789074e-4 * t64 + t389 * 9.927576e-7 * t64);
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
		    t7 = 1 / t5 / t4;
		    t8 = sigmabb * t7;
		    t10 = t8 * .004 + 1.;
		    t11 = 1 / t10;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t18 = 1 / t2 / t15 / rhob;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t8 * .0020292 * t11 + .8094 + t19 * 1.19696e-5 * 
			    t21;
		    t27 = 1 / rhob;
		    t28 = pow_dd(&t27, &c_b2);
		    t30 = t28 * .1274696188700087 + 1.;
		    t31 = rhob * t30;
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t39 = t32 * 11.12037486309468 + t28 * 3.844746237447211 + 
			    t35 * 1.644733775567609 + t37 * .2405871291288192;
		    t42 = 32.1646831778707 / t39 + 1.;
		    t43 = log(t42);
		    t45 = t8 * .2 + 1.;
		    t46 = 1 / t45;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t8 * .46974 * t46 + .1737 - t19 * .099472 * t50;
		    t54 = t43 * t53;
		    zk[i__] = t3 * -.9305257363491 * t24 - t31 * .03109 * t54;
		    vrhoa[i__] = 0.;
		    t59 = t4 * rhob;
		    t62 = sigmabb / t5 / t59;
		    t68 = t14 / t2 / t15 / t4;
		    t71 = t14 * sigmabb;
/* Computing 2nd power */
		    d__1 = t15;
		    t72 = d__1 * d__1;
		    t75 = t71 / t72 / rhob;
		    t77 = 1 / t20 / t10;
		    t80 = t62 * -.0054112 * t11 - t68 * 4.219306666666667e-5 *
			     t21 + t75 * 2.553514666666667e-7 * t77;
		    t83 = t30 * t43;
		    t86 = 1 / t37;
		    t87 = t27 * t86;
/* Computing 2nd power */
		    d__1 = t39;
		    t90 = d__1 * d__1;
		    t91 = 1 / t90;
		    t92 = t31 * t91;
/* Computing 2nd power */
		    d__1 = t32;
		    t93 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t93;
		    t94 = d__1 * d__1;
		    t95 = t94 * t32;
		    t96 = 1 / t95;
		    t97 = 1 / t4;
		    t102 = 1 / t35;
		    t105 = 1 / t28;
		    t108 = t96 * -1.853395810515781 * t97 - t86 * 
			    1.28158207914907 * t97 - t102 * .8223668877838045 
			    * t97 - t105 * .1603914194192128 * t97;
		    t109 = 1 / t42;
		    t110 = t108 * t109;
		    t111 = t110 * t53;
		    t119 = 1 / t49 / t45;
		    t122 = t62 * -1.25264 * t46 + t68 * .7810453333333333 * 
			    t50 - t75 * .1061034666666667 * t119;
		    t123 = t43 * t122;
		    vrhob[i__] = t2 * -1.2407009817988 * t24 - t3 * 
			    .9305257363491 * t80 - t83 * .03109 * t53 + t87 * 
			    .001321010150222857 * t54 + t92 * 1. * t111 - t31 
			    * .03109 * t123;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t128 = sigmabb * t18;
		    t131 = 1 / t72;
		    t132 = t14 * t131;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * .0020292 * 
			    t11 + t128 * 1.58224e-5 * t21 - t132 * 9.57568e-8 
			    * t77) - t31 * .03109 * t43 * (t7 * .46974 * t46 
			    - t128 * .292892 * t50 + t132 * .0397888 * t119);
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t155 = sigmabb / t5 / t15;
		    t161 = t14 / t2 / t15 / t59;
		    t164 = t72 * t4;
		    t166 = t71 / t164;
/* Computing 2nd power */
		    d__1 = t14;
		    t169 = d__1 * d__1;
		    t173 = t169 / t5 / t72 / t15;
/* Computing 2nd power */
		    d__1 = t20;
		    t174 = d__1 * d__1;
		    t175 = 1 / t174;
		    t186 = 1 / t59;
		    t188 = 1 / t37 / t27;
/* Computing 2nd power */
		    d__1 = t108;
		    t200 = d__1 * d__1;
		    t207 = 1 / t15;
/* Computing 2nd power */
		    d__1 = t90;
		    t233 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t42;
		    t236 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t251 = d__1 * d__1;
		    t252 = 1 / t251;
		    s1 = -.4135669939329333 / t5 * t24 - t2 * 2.4814019635976 
			    * t80 - t3 * .9305257363491 * (t155 * 
			    .01984106666666667 * t11 + t161 * 
			    2.095032888888889e-4 * t21 - t166 * 
			    3.198281955555556e-6 * t77 + t173 * 
			    8.171246933333333e-9 * t175) + t30 * 2. * t91 * 
			    t111 - t83 * .06218 * t122 + t186 * 
			    8.806734334819047e-4 * t188 * t54;
		    v2rhob2[i__] = s1 - t87 * .08497974591333914 * t91 * t111 
			    + t87 * .002642020300445714 * t123 - t31 * 2. / 
			    t90 / t39 * t200 * t109 * t53 + t92 * 1. * (
			    -1.544496508763151 / t95 / t27 * t207 + t96 * 
			    3.706791621031562 * t186 - t188 * 
			    .854388052766047 * t207 + t86 * 2.563164158298141 
			    * t186 - .4111834438919023 / t35 / t27 * t207 + 
			    t102 * 1.644733775567609 * t186 - 
			    .05346380647307093 / t28 / t27 * t207 + t105 * 
			    .3207828388384256 * t186) * t109 * t53 + t31 * 
			    32.1646831778707 / t233 * t200 / t236 * t53 + t92 
			    * 2. * t110 * t122 - t31 * .03109 * t43 * (t155 * 
			    4.593013333333333 * t46 - t161 * 
			    5.614695111111111 * t50 + t166 * 
			    1.788046222222222 * t119 - t173 * 
			    .1697655466666667 * t252);
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t261 = sigmabb * t131;
		    t266 = t14 / t5 / t164;
		    v2sigmabb2[i__] = t3 * -.9305257363491 * (t18 * 7.7056e-6 
			    * t21 - t261 * 3.180928e-7 * t77 + t266 * 
			    1.1490816e-9 * t175) - t31 * .03109 * t43 * (t18 *
			     -.38684 * t50 + t261 * .1967344 * t119 - t266 * 
			    .02387328 * t252);
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
		    t7 = 1 / t5 / t4;
		    t8 = sigmaaa * t7;
		    t10 = t8 * .004 + 1.;
		    t11 = 1 / t10;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t14 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t4;
		    t15 = d__1 * d__1;
		    t18 = 1 / t2 / t15 / rhoa;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t8 * .0020292 * t11 + .8094 + t19 * 1.19696e-5 * 
			    t21;
		    t27 = 1 / rhoa;
		    t28 = pow_dd(&t27, &c_b2);
		    t30 = t28 * .1274696188700087 + 1.;
		    t31 = rhoa * t30;
		    t32 = pow_dd(&t27, &c_b4);
		    t35 = sqrt(t27);
/* Computing 2nd power */
		    d__1 = t28;
		    t37 = d__1 * d__1;
		    t39 = t32 * 11.12037486309468 + t28 * 3.844746237447211 + 
			    t35 * 1.644733775567609 + t37 * .2405871291288192;
		    t42 = 32.1646831778707 / t39 + 1.;
		    t43 = log(t42);
		    t45 = t8 * .2 + 1.;
		    t46 = 1 / t45;
/* Computing 2nd power */
		    d__1 = t45;
		    t49 = d__1 * d__1;
		    t50 = 1 / t49;
		    t53 = t8 * .46974 * t46 + .1737 - t19 * .099472 * t50;
		    t54 = t43 * t53;
		    zk[i__] = t3 * -.9305257363491 * t24 - t31 * .03109 * t54;
		    t59 = t4 * rhoa;
		    t62 = sigmaaa / t5 / t59;
		    t68 = t14 / t2 / t15 / t4;
		    t71 = t14 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t15;
		    t72 = d__1 * d__1;
		    t75 = t71 / t72 / rhoa;
		    t77 = 1 / t20 / t10;
		    t80 = t62 * -.0054112 * t11 - t68 * 4.219306666666667e-5 *
			     t21 + t75 * 2.553514666666667e-7 * t77;
		    t83 = t30 * t43;
		    t86 = 1 / t37;
		    t87 = t27 * t86;
/* Computing 2nd power */
		    d__1 = t39;
		    t90 = d__1 * d__1;
		    t91 = 1 / t90;
		    t92 = t31 * t91;
/* Computing 2nd power */
		    d__1 = t32;
		    t93 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t93;
		    t94 = d__1 * d__1;
		    t95 = t94 * t32;
		    t96 = 1 / t95;
		    t97 = 1 / t4;
		    t102 = 1 / t35;
		    t105 = 1 / t28;
		    t108 = t96 * -1.853395810515781 * t97 - t86 * 
			    1.28158207914907 * t97 - t102 * .8223668877838045 
			    * t97 - t105 * .1603914194192128 * t97;
		    t109 = 1 / t42;
		    t110 = t108 * t109;
		    t111 = t110 * t53;
		    t119 = 1 / t49 / t45;
		    t122 = t62 * -1.25264 * t46 + t68 * .7810453333333333 * 
			    t50 - t75 * .1061034666666667 * t119;
		    t123 = t43 * t122;
		    vrhoa[i__] = t2 * -1.2407009817988 * t24 - t3 * 
			    .9305257363491 * t80 - t83 * .03109 * t53 + t87 * 
			    .001321010150222857 * t54 + t92 * 1. * t111 - t31 
			    * .03109 * t123;
		    vrhob[i__] = 0.;
		    t128 = sigmaaa * t18;
		    t131 = 1 / t72;
		    t132 = t14 * t131;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * .0020292 * 
			    t11 + t128 * 1.58224e-5 * t21 - t132 * 9.57568e-8 
			    * t77) - t31 * .03109 * t43 * (t7 * .46974 * t46 
			    - t128 * .292892 * t50 + t132 * .0397888 * t119);
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t155 = sigmaaa / t5 / t15;
		    t161 = t14 / t2 / t15 / t59;
		    t164 = t72 * t4;
		    t166 = t71 / t164;
/* Computing 2nd power */
		    d__1 = t14;
		    t169 = d__1 * d__1;
		    t173 = t169 / t5 / t72 / t15;
/* Computing 2nd power */
		    d__1 = t20;
		    t174 = d__1 * d__1;
		    t175 = 1 / t174;
		    t186 = 1 / t59;
		    t188 = 1 / t37 / t27;
/* Computing 2nd power */
		    d__1 = t108;
		    t200 = d__1 * d__1;
		    t207 = 1 / t15;
/* Computing 2nd power */
		    d__1 = t90;
		    t233 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t42;
		    t236 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t251 = d__1 * d__1;
		    t252 = 1 / t251;
		    s1 = -.4135669939329333 / t5 * t24 - t2 * 2.4814019635976 
			    * t80 - t3 * .9305257363491 * (t155 * 
			    .01984106666666667 * t11 + t161 * 
			    2.095032888888889e-4 * t21 - t166 * 
			    3.198281955555556e-6 * t77 + t173 * 
			    8.171246933333333e-9 * t175) + t30 * 2. * t91 * 
			    t111 - t83 * .06218 * t122 + t186 * 
			    8.806734334819047e-4 * t188 * t54;
		    v2rhoa2[i__] = s1 - t87 * .08497974591333914 * t91 * t111 
			    + t87 * .002642020300445714 * t123 - t31 * 2. / 
			    t90 / t39 * t200 * t109 * t53 + t92 * 1. * (
			    -1.544496508763151 / t95 / t27 * t207 + t96 * 
			    3.706791621031562 * t186 - t188 * 
			    .854388052766047 * t207 + t86 * 2.563164158298141 
			    * t186 - .4111834438919023 / t35 / t27 * t207 + 
			    t102 * 1.644733775567609 * t186 - 
			    .05346380647307093 / t28 / t27 * t207 + t105 * 
			    .3207828388384256 * t186) * t109 * t53 + t31 * 
			    32.1646831778707 / t233 * t200 / t236 * t53 + t92 
			    * 2. * t110 * t122 - t31 * .03109 * t43 * (t155 * 
			    4.593013333333333 * t46 - t161 * 
			    5.614695111111111 * t50 + t166 * 
			    1.788046222222222 * t119 - t173 * 
			    .1697655466666667 * t252);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t261 = sigmaaa * t131;
		    t266 = t14 / t5 / t164;
		    v2sigmaaa2[i__] = t3 * -.9305257363491 * (t18 * 7.7056e-6 
			    * t21 - t261 * 3.180928e-7 * t77 + t266 * 
			    1.1490816e-9 * t175) - t31 * .03109 * t43 * (t18 *
			     -.38684 * t50 + t261 * .1967344 * t119 - t266 * 
			    .02387328 * t252);
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
		    t9 = 1 / t7 / t6;
		    t10 = sigmaaa * t9;
		    t12 = t10 * .004 + 1.;
		    t13 = 1 / t12;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t16 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t17 = d__1 * d__1;
		    t20 = 1 / t4 / t17 / rhoa;
		    t21 = t16 * t20;
/* Computing 2nd power */
		    d__1 = t12;
		    t22 = d__1 * d__1;
		    t23 = 1 / t22;
		    t26 = t10 * .0020292 * t13 + .8094 + t21 * 1.19696e-5 * 
			    t23;
		    t29 = 1 / rhoa;
		    t30 = pow_dd(&t29, &c_b2);
		    t32 = t30 * .1274696188700087 + 1.;
		    t33 = rhoa * t32;
		    t34 = pow_dd(&t29, &c_b4);
		    t37 = sqrt(t29);
/* Computing 2nd power */
		    d__1 = t30;
		    t39 = d__1 * d__1;
		    t41 = t34 * 11.12037486309468 + t30 * 3.844746237447211 + 
			    t37 * 1.644733775567609 + t39 * .2405871291288192;
		    t44 = 32.1646831778707 / t41 + 1.;
		    t45 = log(t44);
		    t47 = t10 * .2 + 1.;
		    t48 = 1 / t47;
/* Computing 2nd power */
		    d__1 = t47;
		    t51 = d__1 * d__1;
		    t52 = 1 / t51;
		    t55 = t10 * .46974 * t48 + .1737 - t21 * .099472 * t52;
		    t56 = t45 * t55;
		    t59 = pow_dd(&rhob, &c_b2);
		    t60 = t59 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t61 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t59;
		    t62 = d__1 * d__1;
		    t64 = 1 / t62 / t61;
		    t65 = sigmabb * t64;
		    t67 = t65 * .004 + 1.;
		    t68 = 1 / t67;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t71 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t61;
		    t72 = d__1 * d__1;
		    t75 = 1 / t59 / t72 / rhob;
		    t76 = t71 * t75;
/* Computing 2nd power */
		    d__1 = t67;
		    t77 = d__1 * d__1;
		    t78 = 1 / t77;
		    t81 = t65 * .0020292 * t68 + .8094 + t76 * 1.19696e-5 * 
			    t78;
		    t84 = 1 / rhob;
		    t85 = pow_dd(&t84, &c_b2);
		    t87 = t85 * .1274696188700087 + 1.;
		    t88 = rhob * t87;
		    t89 = pow_dd(&t84, &c_b4);
		    t92 = sqrt(t84);
/* Computing 2nd power */
		    d__1 = t85;
		    t94 = d__1 * d__1;
		    t96 = t89 * 11.12037486309468 + t85 * 3.844746237447211 + 
			    t92 * 1.644733775567609 + t94 * .2405871291288192;
		    t99 = 32.1646831778707 / t96 + 1.;
		    t100 = log(t99);
		    t102 = t65 * .2 + 1.;
		    t103 = 1 / t102;
/* Computing 2nd power */
		    d__1 = t102;
		    t106 = d__1 * d__1;
		    t107 = 1 / t106;
		    t110 = t65 * .46974 * t103 + .1737 - t76 * .099472 * t107;
		    t111 = t100 * t110;
		    t114 = rhoa + rhob;
		    t115 = 1 / t114;
		    t116 = pow_dd(&t115, &c_b2);
		    t118 = t116 * .1325688999052018 + 1.;
		    t119 = pow_dd(&t115, &c_b4);
		    t122 = sqrt(t115);
/* Computing 2nd power */
		    d__1 = t116;
		    t124 = d__1 * d__1;
		    t126 = t119 * 5.98255043577108 + t116 * 2.225569421150687 
			    + t122 * .8004286349993634 + t124 * 
			    .1897004325747559;
		    t129 = 16.0818243221511 / t126 + 1.;
		    t130 = log(t129);
		    t132 = t118 * .062182 * t130;
		    t134 = t116 * .06901399211255825 + 1.;
		    t139 = t119 * 8.157414703487641 + t116 * 
			    2.247591863577616 + t122 * .4300972471276643 + 
			    t124 * .1911512595127338;
		    t142 = 29.60857464321668 / t139 + 1.;
		    t143 = log(t142);
		    t144 = t134 * t143;
		    t146 = rhoa - rhob * 1.;
		    t147 = t146 * t115;
		    t148 = t147 + 1.;
		    t149 = pow_dd(&t148, &c_b2);
		    t152 = 1. - t147 * 1.;
		    t153 = pow_dd(&t152, &c_b2);
		    t155 = t149 * t148 + t153 * t152 - 2.;
/* Computing 2nd power */
		    d__1 = t146;
		    t156 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t156;
		    t157 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t114;
		    t158 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t158;
		    t159 = d__1 * d__1;
		    t160 = 1 / t159;
		    t161 = t157 * t160;
		    t163 = 1. - t161 * 1.;
		    t164 = t155 * t163;
		    t166 = t144 * .03799574853701528 * t164;
		    t168 = t116 * .1274696188700087 + 1.;
		    t173 = t119 * 11.12037486309468 + t116 * 
			    3.844746237447211 + t122 * 1.644733775567609 + 
			    t124 * .2405871291288192;
		    t176 = 32.1646831778707 / t173 + 1.;
		    t177 = log(t176);
		    t180 = t168 * -.03109 * t177 + t132;
		    t181 = t180 * t155;
		    t183 = t181 * 1.923661050931536 * t161;
		    t190 = t114 * (-t132 + t166 + t183) + t33 * .03109 * t45 
			    + t88 * .03109 * t100;
		    t193 = t10 * .5 + t65 * .5;
		    t196 = t10 * .003 + 1. + t65 * .003;
		    t197 = 1 / t196;
/* Computing 2nd power */
		    d__1 = t193;
		    t200 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t196;
		    t201 = d__1 * d__1;
		    t202 = 1 / t201;
		    t205 = t193 * .0044826 * t197 + .9454 - t200 * 
			    1.654596e-4 * t202;
		    zk[i__] = t5 * -.9305257363491 * t26 - t33 * .03109 * t56 
			    - t60 * .9305257363491 * t81 - t88 * .03109 * 
			    t111 + t190 * t205;
		    t209 = t6 * rhoa;
		    t211 = 1 / t7 / t209;
		    t212 = sigmaaa * t211;
		    t217 = 1 / t4 / t17 / t6;
		    t218 = t16 * t217;
		    t221 = t16 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t17;
		    t222 = d__1 * d__1;
		    t224 = 1 / t222 / rhoa;
		    t225 = t221 * t224;
		    t227 = 1 / t22 / t12;
		    t230 = t212 * -.0054112 * t13 - t218 * 
			    4.219306666666667e-5 * t23 + t225 * 
			    2.553514666666667e-7 * t227;
		    t233 = t32 * t45;
		    t236 = 1 / t39;
		    t237 = t29 * t236;
/* Computing 2nd power */
		    d__1 = t41;
		    t240 = d__1 * d__1;
		    t241 = 1 / t240;
		    t242 = t33 * t241;
/* Computing 2nd power */
		    d__1 = t34;
		    t243 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t243;
		    t244 = d__1 * d__1;
		    t245 = t244 * t34;
		    t246 = 1 / t245;
		    t247 = 1 / t6;
		    t252 = 1 / t37;
		    t255 = 1 / t30;
		    t258 = t246 * -1.853395810515781 * t247 - t236 * 
			    1.28158207914907 * t247 - t252 * 
			    .8223668877838045 * t247 - t255 * 
			    .1603914194192128 * t247;
		    t259 = 1 / t44;
		    t260 = t258 * t259;
		    t261 = t260 * t55;
		    t269 = 1 / t51 / t47;
		    t272 = t212 * -1.25264 * t48 + t218 * .7810453333333333 * 
			    t52 - t225 * .1061034666666667 * t269;
		    t273 = t45 * t272;
		    t276 = 1 / t124;
		    t277 = 1 / t158;
		    t278 = t276 * t277;
		    t279 = t278 * t130;
		    t280 = t279 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t126;
		    t281 = d__1 * d__1;
		    t282 = 1 / t281;
		    t283 = t118 * t282;
/* Computing 2nd power */
		    d__1 = t119;
		    t284 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t284;
		    t285 = d__1 * d__1;
		    t286 = t285 * t119;
		    t287 = 1 / t286;
		    t288 = t287 * t277;
		    t291 = 1 / t122;
		    t292 = t291 * t277;
		    t294 = 1 / t116;
		    t295 = t294 * t277;
		    t297 = t288 * -.99709173929518 - t278 * .7418564737168958 
			    - t292 * .4002143174996817 - t295 * 
			    .1264669550498372;
		    t298 = 1 / t129;
		    t300 = t283 * t297 * t298;
		    t301 = t300 * 1.;
		    t302 = t143 * t155;
		    t303 = t302 * t163;
		    t304 = t278 * t303;
		    t305 = t304 * 8.740794299481065e-4;
/* Computing 2nd power */
		    d__1 = t139;
		    t306 = d__1 * d__1;
		    t307 = 1 / t306;
		    t308 = t134 * t307;
		    t313 = t288 * -1.35956911724794 - t278 * 
			    .7491972878592054 - t292 * .2150486235638321 - 
			    t295 * .1274341730084892;
		    t314 = t308 * t313;
		    t315 = 1 / t142;
		    t316 = t315 * t155;
		    t317 = t316 * t163;
		    t318 = t314 * t317;
		    t319 = t318 * 1.124999956683108;
		    t320 = t146 * t277;
		    t321 = t320 * 1.;
		    t322 = t115 - t321;
		    t325 = t115 * 1.;
		    t326 = -t325 + t320;
		    t329 = t149 * 1.333333333333333 * t322 + t153 * 
			    1.333333333333333 * t326;
		    t331 = t144 * t329 * t163;
		    t332 = t331 * .03799574853701528;
		    t333 = t156 * t146;
		    t334 = t333 * t160;
		    t335 = t334 * 4.;
		    t337 = 1 / t159 / t114;
		    t338 = t157 * t337;
		    t339 = t338 * 4.;
		    t340 = -t335 + t339;
		    t342 = t144 * t155 * t340;
		    t343 = t342 * .03799574853701528;
/* Computing 2nd power */
		    d__1 = t173;
		    t346 = d__1 * d__1;
		    t347 = 1 / t346;
		    t348 = t168 * t347;
		    t353 = t288 * -1.853395810515781 - t278 * 
			    1.28158207914907 - t292 * .8223668877838045 - 
			    t295 * .1603914194192128;
		    t354 = 1 / t176;
		    t360 = t278 * .001321010150222857 * t177 + t348 * 1. * 
			    t353 * t354 - t279 * .002747799777968419 - t300 * 
			    1.;
		    t361 = t360 * t155;
		    t362 = t361 * t161;
		    t363 = t362 * 1.923661050931536;
		    t364 = t180 * t329;
		    t365 = t364 * t161;
		    t366 = t365 * 1.923661050931536;
		    t367 = t181 * t334;
		    t369 = t181 * t338;
		    t370 = t369 * 7.694644203726145;
		    t377 = t241 * t258 * t259;
		    t380 = -t132 + t166 + t183 + t114 * (t280 + t301 - t305 - 
			    t319 + t332 + t343 + t363 + t366 + t367 * 
			    7.694644203726145 - t370) + t233 * .03109 - t237 *
			     .001321010150222857 * t45 - t33 * 1. * t377;
		    t384 = t193 * t202;
		    t388 = 1 / t201 / t196;
		    t389 = t200 * t388;
		    t392 = t212 * -.0059768 * t197 + t384 * 4.770864e-4 * 
			    t212 - t389 * 2.6473536e-6 * t212;
		    vrhoa[i__] = t4 * -1.2407009817988 * t26 - t5 * 
			    .9305257363491 * t230 - t233 * .03109 * t55 + 
			    t237 * .001321010150222857 * t56 + t242 * 1. * 
			    t261 - t33 * .03109 * t273 + t380 * t205 + t190 * 
			    t392;
		    t396 = t61 * rhob;
		    t398 = 1 / t62 / t396;
		    t399 = sigmabb * t398;
		    t404 = 1 / t59 / t72 / t61;
		    t405 = t71 * t404;
		    t408 = t71 * sigmabb;
/* Computing 2nd power */
		    d__1 = t72;
		    t409 = d__1 * d__1;
		    t411 = 1 / t409 / rhob;
		    t412 = t408 * t411;
		    t414 = 1 / t77 / t67;
		    t417 = t399 * -.0054112 * t68 - t405 * 
			    4.219306666666667e-5 * t78 + t412 * 
			    2.553514666666667e-7 * t414;
		    t420 = t87 * t100;
		    t423 = 1 / t94;
		    t424 = t84 * t423;
/* Computing 2nd power */
		    d__1 = t96;
		    t427 = d__1 * d__1;
		    t428 = 1 / t427;
		    t429 = t88 * t428;
/* Computing 2nd power */
		    d__1 = t89;
		    t430 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t430;
		    t431 = d__1 * d__1;
		    t432 = t431 * t89;
		    t433 = 1 / t432;
		    t434 = 1 / t61;
		    t439 = 1 / t92;
		    t442 = 1 / t85;
		    t445 = t433 * -1.853395810515781 * t434 - t423 * 
			    1.28158207914907 * t434 - t439 * 
			    .8223668877838045 * t434 - t442 * 
			    .1603914194192128 * t434;
		    t446 = 1 / t99;
		    t447 = t445 * t446;
		    t448 = t447 * t110;
		    t456 = 1 / t106 / t102;
		    t459 = t399 * -1.25264 * t103 + t405 * .7810453333333333 *
			     t107 - t412 * .1061034666666667 * t456;
		    t460 = t100 * t459;
		    t463 = -t325 - t321;
		    t466 = t115 + t320;
		    t469 = t149 * 1.333333333333333 * t463 + t153 * 
			    1.333333333333333 * t466;
		    t471 = t144 * t469 * t163;
		    t472 = t471 * .03799574853701528;
		    t473 = t335 + t339;
		    t475 = t144 * t155 * t473;
		    t476 = t475 * .03799574853701528;
		    t477 = t180 * t469;
		    t478 = t477 * t161;
		    t479 = t478 * 1.923661050931536;
		    t487 = t428 * t445 * t446;
		    t490 = -t132 + t166 + t183 + t114 * (t280 + t301 - t305 - 
			    t319 + t472 + t476 + t363 + t479 - t367 * 
			    7.694644203726145 - t370) + t420 * .03109 - t424 *
			     .001321010150222857 * t100 - t88 * 1. * t487;
		    t498 = t399 * -.0059768 * t197 + t384 * 4.770864e-4 * 
			    t399 - t389 * 2.6473536e-6 * t399;
		    vrhob[i__] = t59 * -1.2407009817988 * t81 - t60 * 
			    .9305257363491 * t417 - t420 * .03109 * t110 + 
			    t424 * .001321010150222857 * t111 + t429 * 1. * 
			    t448 - t88 * .03109 * t460 + t490 * t205 + t190 * 
			    t498;
		    t502 = sigmaaa * t20;
		    t505 = 1 / t222;
		    t506 = t16 * t505;
		    t509 = t9 * .0020292 * t13 + t502 * 1.58224e-5 * t23 - 
			    t506 * 9.57568e-8 * t227;
		    t518 = t9 * .46974 * t48 - t502 * .292892 * t52 + t506 * 
			    .0397888 * t269;
		    t519 = t45 * t518;
		    t528 = t9 * .0022413 * t197 - t384 * 1.789074e-4 * t9 + 
			    t389 * 9.927576e-7 * t9;
		    vsigmaaa[i__] = t5 * -.9305257363491 * t509 - t33 * 
			    .03109 * t519 + t190 * t528;
		    vsigmaab[i__] = 0.;
		    t532 = sigmabb * t75;
		    t535 = 1 / t409;
		    t536 = t71 * t535;
		    t539 = t64 * .0020292 * t68 + t532 * 1.58224e-5 * t78 - 
			    t536 * 9.57568e-8 * t414;
		    t548 = t64 * .46974 * t103 - t532 * .292892 * t107 + t536 
			    * .0397888 * t456;
		    t549 = t100 * t548;
		    t558 = t64 * .0022413 * t197 - t384 * 1.789074e-4 * t64 + 
			    t389 * 9.927576e-7 * t64;
		    vsigmabb[i__] = t60 * -.9305257363491 * t539 - t88 * 
			    .03109 * t549 + t190 * t558;
		    t562 = t32 * t241;
		    t566 = 1 / t240 / t41;
/* Computing 2nd power */
		    d__1 = t258;
		    t567 = d__1 * d__1;
		    t574 = 1 / t209;
		    t576 = 1 / t39 / t29;
		    t577 = t574 * t576;
		    t582 = 1 / t17;
		    t603 = -1.544496508763151 / t245 / t29 * t582 + t246 * 
			    3.706791621031562 * t574 - t576 * 
			    .854388052766047 * t582 + t236 * 
			    2.563164158298141 * t574 - .4111834438919023 / 
			    t37 / t29 * t582 + t252 * 1.644733775567609 * 
			    t574 - .05346380647307093 / t30 / t29 * t582 + 
			    t255 * .3207828388384256 * t574;
/* Computing 2nd power */
		    d__1 = t240;
		    t609 = d__1 * d__1;
		    t610 = 1 / t609;
/* Computing 2nd power */
		    d__1 = t44;
		    t612 = d__1 * d__1;
		    t613 = 1 / t612;
		    t617 = t279 * .005495599555936838;
		    t618 = t300 * 2.;
		    t621 = t369 * 15.38928840745229;
		    t623 = t304 * .001748158859896213;
		    t624 = t318 * 2.249999913366216;
		    t625 = t362 * 3.847322101863073;
		    t628 = 1 / t124 / t115 * t160;
		    t629 = t628 * t130;
		    t630 = t629 * .001831866518645613;
		    t632 = 1 / t158 / t114;
		    t633 = t276 * t632;
		    t634 = t633 * t130;
		    t635 = t634 * .005495599555936838;
		    t638 = t278 * t282 * t297 * t298;
		    t639 = t638 * .08837926660346786;
		    t641 = t278 * t302 * t340;
/* Computing 2nd power */
		    d__1 = t281;
		    t643 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t297;
		    t646 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t129;
		    t647 = d__1 * d__1;
		    t650 = t118 / t643 * t646 / t647;
		    t651 = t650 * 16.0818243221511;
		    t652 = t156 * t160;
		    t653 = t652 * 12.;
		    t654 = t333 * t337;
		    t655 = t654 * 32.;
		    t658 = t157 / t159 / t158;
		    t659 = t658 * 20.;
/* Computing 2nd power */
		    d__1 = t149;
		    t664 = d__1 * d__1;
		    t665 = 1 / t664;
/* Computing 2nd power */
		    d__1 = t322;
		    t666 = d__1 * d__1;
		    t669 = t277 * 2.;
		    t671 = t146 * 2. * t632;
		    t672 = -t669 + t671;
/* Computing 2nd power */
		    d__1 = t153;
		    t675 = d__1 * d__1;
		    t676 = 1 / t675;
/* Computing 2nd power */
		    d__1 = t326;
		    t677 = d__1 * d__1;
		    t683 = t665 * .4444444444444444 * t666 + t149 * 
			    1.333333333333333 * t672 + t676 * 
			    .4444444444444444 * t677 - t153 * 
			    1.333333333333333 * t672;
		    t689 = t278 * t143 * t329 * t163;
		    t692 = t628 * 5.827196199654043e-4 * t303;
		    t694 = t633 * .001748158859896213 * t303;
		    t699 = t278 * .05176049209143758 * t307 * t313 * t315 * 
			    t164;
/* Computing 2nd power */
		    d__1 = t313;
		    t703 = d__1 * d__1;
		    t706 = t134 * 2.249999913366216 / t306 / t139 * t703 * 
			    t317;
		    t709 = 1 / t286 / t115 * t160;
		    t711 = t287 * t632;
		    t717 = 1 / t122 / t115 * t160;
		    t719 = t291 * t632;
		    t723 = 1 / t116 / t115 * t160;
		    t725 = t294 * t632;
		    t730 = t308 * 1.124999956683108 * (t709 * 
			    -1.132974264373283 + t711 * 2.71913823449588 - 
			    t628 * .4994648585728036 + t633 * 
			    1.498394575718411 - t717 * .1075243117819161 + 
			    t719 * .4300972471276643 - t723 * 
			    .04247805766949639 + t725 * .2548683460169784) * 
			    t317;
		    t733 = t314 * t315 * t329 * t163;
		    t735 = t630 - t635 - t639 - t641 * .001748158859896213 + 
			    t651 + t144 * .03799574853701528 * t155 * (-t653 
			    + t655 - t659) + t144 * .03799574853701528 * t683 
			    * t163 - t689 * .001748158859896213 - t692 + t694 
			    + t699 + t706 - t730 - t733 * 2.249999913366216;
		    t737 = t314 * t316 * t340;
		    t742 = t364 * t338;
		    t744 = t364 * t334;
		    t746 = t181 * t654;
		    t748 = t181 * t652;
		    t749 = t748 * 23.08393261117844;
/* Computing 2nd power */
		    d__1 = t353;
		    t761 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t346;
		    t777 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t176;
		    t780 = d__1 * d__1;
		    t792 = t118 / t281 / t126 * t646 * t298;
		    t804 = t283 * (t709 * -.8309097827459833 + t711 * 
			    1.99418347859036 - t628 * .4945709824779306 + 
			    t633 * 1.483712947433792 - t717 * 
			    .2001071587498409 + t719 * .8004286349993634 - 
			    t723 * .04215565168327908 + t725 * 
			    .2529339100996745) * t298;
		    t807 = t628 * 8.806734334819047e-4 * t177 - t633 * 
			    .002642020300445714 * t177 - t278 * 
			    .08497974591333914 * t347 * t353 * t354 - t168 * 
			    2. / t346 / t173 * t761 * t354 + t348 * 1. * (
			    t709 * -1.544496508763151 + t711 * 
			    3.706791621031562 - t628 * .854388052766047 + 
			    t633 * 2.563164158298141 - t717 * 
			    .4111834438919023 + t719 * 1.644733775567609 - 
			    t723 * .05346380647307093 + t725 * 
			    .3207828388384256) * t354 + t168 * 
			    32.1646831778707 / t777 * t761 / t780 - t629 * 
			    .001831866518645613 + t634 * .005495599555936838 
			    + t638 * .08837926660346786 + t792 * 2. - t804 * 
			    1. - t650 * 16.0818243221511;
		    t810 = t807 * 1.923661050931536 * t155 * t161;
		    t812 = t181 * 38.47322101863073 * t658;
		    t817 = t360 * t329 * t161;
		    t819 = t804 * 1.;
/* Computing 2nd power */
		    d__1 = t306;
		    t820 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t142;
		    t824 = d__1 * d__1;
		    t829 = t134 * 33.30964519106732 / t820 * t703 / t824 * 
			    t155 * t163;
		    t830 = t792 * 2.;
		    t831 = t361 * t334;
		    t834 = t361 * 15.38928840745229 * t338;
		    t835 = t737 * -2.249999913366216 + t180 * 
			    1.923661050931536 * t683 * t161 - t742 * 
			    15.38928840745229 + t744 * 15.38928840745229 - 
			    t746 * 61.55715362980916 + t749 + t810 + t812 + 
			    t144 * .07599149707403056 * t329 * t340 + t817 * 
			    3.847322101863073 + t819 - t829 - t830 + t831 * 
			    15.38928840745229 - t834;
		    t838 = t562 * -2. * t260 + t33 * 2. * t566 * t567 * t259 
			    + t237 * .08497974591333914 * t377 - t577 * 
			    8.806734334819047e-4 * t45 - t33 * 1. * t241 * 
			    t603 * t259 + t331 * .07599149707403056 - t33 * 
			    32.1646831778707 * t610 * t567 * t613 + t617 + 
			    t618 + t365 * 3.847322101863073 + t367 * 
			    15.38928840745229 - t621 + t342 * 
			    .07599149707403056 - t623 - t624 + t625 + t114 * (
			    t735 + t835);
		    t842 = sigmaaa / t7 / t17;
		    t848 = t16 / t4 / t17 / t209;
		    t851 = t193 * t388;
/* Computing 2nd power */
		    d__1 = t201;
		    t856 = d__1 * d__1;
		    t858 = t200 / t856;
		    t876 = t222 * t6;
		    t878 = t221 / t876;
/* Computing 2nd power */
		    d__1 = t16;
		    t881 = d__1 * d__1;
		    t885 = t881 / t7 / t222 / t17;
/* Computing 2nd power */
		    d__1 = t22;
		    t886 = d__1 * d__1;
		    t887 = 1 / t886;
/* Computing 2nd power */
		    d__1 = t51;
		    t899 = d__1 * d__1;
		    t900 = 1 / t899;
		    s1 = t380 * 2. * t392 + t838 * t205 + t190 * (t842 * 
			    .02191493333333333 * t197 - t848 * 6.839296e-4 * 
			    t202 + t851 * 1.4692992e-5 * t848 - t384 * 
			    .0017493168 * t842 - t858 * 6.35364864e-8 * t848 
			    + t389 * 9.7069632e-6 * t842) - t233 * .06218 * 
			    t272 - .4135669939329333 / t7 * t26 - t4 * 
			    2.4814019635976 * t230 - t5 * .9305257363491 * (
			    t842 * .01984106666666667 * t13 + t848 * 
			    2.095032888888889e-4 * t23 - t878 * 
			    3.198281955555556e-6 * t227 + t885 * 
			    8.171246933333333e-9 * t887) - t33 * .03109 * t45 
			    * (t842 * 4.593013333333333 * t48 - t848 * 
			    5.614695111111111 * t52 + t878 * 
			    1.788046222222222 * t269 - t885 * 
			    .1697655466666667 * t900);
		    v2rhoa2[i__] = s1 + t562 * 2. * t261 + t237 * 
			    .002642020300445714 * t273 + t577 * 
			    8.806734334819047e-4 * t56 - t237 * 
			    .08497974591333914 * t241 * t261 - t33 * 2. * 
			    t566 * t567 * t259 * t55 + t242 * 1. * t603 * 
			    t259 * t55 + t242 * 2. * t260 * t272 + t33 * 
			    32.1646831778707 * t610 * t567 * t613 * t55;
		    t937 = sigmabb / t62 / t72;
		    t943 = t71 / t59 / t72 / t396;
		    t962 = t409 * t61;
		    t964 = t408 / t962;
/* Computing 2nd power */
		    d__1 = t71;
		    t967 = d__1 * d__1;
		    t971 = t967 / t62 / t409 / t72;
/* Computing 2nd power */
		    d__1 = t77;
		    t972 = d__1 * d__1;
		    t973 = 1 / t972;
		    t987 = 1 / t396;
		    t989 = 1 / t94 / t84;
		    t990 = t987 * t989;
		    t995 = t87 * t428;
		    t1000 = 1 / t72;
		    t1021 = -1.544496508763151 / t432 / t84 * t1000 + t433 * 
			    3.706791621031562 * t987 - t989 * 
			    .854388052766047 * t1000 + t423 * 
			    2.563164158298141 * t987 - .4111834438919023 / 
			    t92 / t84 * t1000 + t439 * 1.644733775567609 * 
			    t987 - .05346380647307093 / t85 / t84 * t1000 + 
			    t442 * .3207828388384256 * t987;
		    t1027 = 1 / t427 / t96;
/* Computing 2nd power */
		    d__1 = t445;
		    t1029 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t427;
		    t1034 = d__1 * d__1;
		    t1035 = 1 / t1034;
/* Computing 2nd power */
		    d__1 = t99;
		    t1037 = d__1 * d__1;
		    t1038 = 1 / t1037;
		    t1068 = t477 * t338;
		    t1070 = t477 * t334;
		    t1073 = t1068 * -15.38928840745229 - t1070 * 
			    15.38928840745229 + t630 - t635 - t639 + t651 - 
			    t692 + t694 + t699 + t706 - t730 + t746 * 
			    61.55715362980916 + t749 + t810;
/* Computing 2nd power */
		    d__1 = t463;
		    t1078 = d__1 * d__1;
		    t1081 = t669 + t671;
/* Computing 2nd power */
		    d__1 = t466;
		    t1084 = d__1 * d__1;
		    t1090 = t665 * .4444444444444444 * t1078 + t149 * 
			    1.333333333333333 * t1081 + t676 * 
			    .4444444444444444 * t1084 - t153 * 
			    1.333333333333333 * t1081;
		    t1099 = t360 * t469 * t161;
		    t1105 = t278 * t302 * t473;
		    t1109 = t278 * t143 * t469 * t163;
		    t1112 = t314 * t316 * t473;
		    t1116 = t314 * t315 * t469 * t163;
		    t1118 = t812 + t819 - t829 - t830 - t831 * 
			    15.38928840745229 - t834 + t144 * 
			    .07599149707403056 * t469 * t473 + t144 * 
			    .03799574853701528 * t1090 * t163 + t144 * 
			    .03799574853701528 * t155 * (-t653 - t655 - t659) 
			    + t1099 * 3.847322101863073 + t180 * 
			    1.923661050931536 * t1090 * t161 - t1105 * 
			    .001748158859896213 - t1109 * .001748158859896213 
			    - t1112 * 2.249999913366216 - t1116 * 
			    2.249999913366216;
		    t1121 = t88 * -1. * t428 * t1021 * t446 + t88 * 2. * 
			    t1027 * t1029 * t446 - t990 * 
			    8.806734334819047e-4 * t100 + t424 * 
			    .08497974591333914 * t487 - t995 * 2. * t447 - 
			    t88 * 32.1646831778707 * t1035 * t1029 * t1038 + 
			    t617 + t618 - t367 * 15.38928840745229 - t621 - 
			    t623 - t624 + t625 + t478 * 3.847322101863073 + 
			    t471 * .07599149707403056 + t475 * 
			    .07599149707403056 + t114 * (t1073 + t1118);
/* Computing 2nd power */
		    d__1 = t106;
		    t1129 = d__1 * d__1;
		    t1130 = 1 / t1129;
		    s1 = t490 * 2. * t498 + t190 * (t937 * .02191493333333333 
			    * t197 - t943 * 6.839296e-4 * t202 + t851 * 
			    1.4692992e-5 * t943 - t384 * .0017493168 * t937 - 
			    t858 * 6.35364864e-8 * t943 + t389 * 9.7069632e-6 
			    * t937) - t420 * .06218 * t459 - t60 * 
			    .9305257363491 * (t937 * .01984106666666667 * t68 
			    + t943 * 2.095032888888889e-4 * t78 - t964 * 
			    3.198281955555556e-6 * t414 + t971 * 
			    8.171246933333333e-9 * t973) - t59 * 
			    2.4814019635976 * t417 - .4135669939329333 / t62 *
			     t81 - t424 * .08497974591333914 * t428 * t448 + 
			    t990 * 8.806734334819047e-4 * t111;
		    v2rhob2[i__] = s1 + t424 * .002642020300445714 * t460 + 
			    t995 * 2. * t448 + t429 * 1. * t1021 * t446 * 
			    t110 - t88 * 2. * t1027 * t1029 * t446 * t110 + 
			    t88 * 32.1646831778707 * t1035 * t1029 * t1038 * 
			    t110 + t429 * 2. * t447 * t459 + t1121 * t205 - 
			    t88 * .03109 * t100 * (t937 * 4.593013333333333 * 
			    t103 - t943 * 5.614695111111111 * t107 + t964 * 
			    1.788046222222222 * t456 - t971 * 
			    .1697655466666667 * t1130);
		    t1160 = t665 * .4444444444444444 * t322 * t463 + t149 * 
			    2.666666666666667 * t146 * t632 + t676 * 
			    .4444444444444444 * t326 * t466 - t153 * 
			    2.666666666666667 * t146 * t632;
		    t1170 = t1068 * -7.694644203726145 + t1070 * 
			    7.694644203726145 + t630 - t635 - t639 - t641 * 
			    8.740794299481065e-4 + t651 - t689 * 
			    8.740794299481065e-4 + t144 * .03799574853701528 *
			     t469 * t340 + t144 * .03799574853701528 * t155 * 
			    (t653 - t659) + t144 * .03799574853701528 * t1160 
			    * t163 + t144 * .03799574853701528 * t329 * t473 
			    + t180 * 1.923661050931536 * t1160 * t161 - t692 
			    + t694 + t699 + t706;
		    t1182 = -t730 - t733 * 1.124999956683108 - t737 * 
			    1.124999956683108 - t742 * 7.694644203726145 - 
			    t744 * 7.694644203726145 - t748 * 
			    23.08393261117844 + t810 + t812 + t817 * 
			    1.923661050931536 + t819 - t829 - t830 - t834 + 
			    t1099 * 1.923661050931536 - t1105 * 
			    8.740794299481065e-4 - t1109 * 
			    8.740794299481065e-4 - t1112 * 1.124999956683108 
			    - t1116 * 1.124999956683108;
		    t1185 = t617 + t618 - t623 - t624 + t472 + t476 + t625 + 
			    t479 - t621 + t332 + t343 + t366 + t114 * (t1170 
			    + t1182);
		    t1195 = t211 * sigmabb * t398;
		    v2rhoab[i__] = t1185 * t205 + t380 * t498 + t490 * t392 + 
			    t190 * (t212 * -6.839296e-4 * t202 * sigmabb * 
			    t398 + t851 * 1.4692992e-5 * sigmaaa * t1195 - 
			    t858 * 6.35364864e-8 * sigmaaa * t1195);
		    t1207 = sigmaaa * t217;
		    t1210 = t16 * t224;
		    t1216 = t221 / t7 / t222 / t209;
		    v2rhoasigmaaa[i__] = t4 * -1.2407009817988 * t509 - t5 * 
			    .9305257363491 * (t211 * -.0054112 * t13 - t1207 *
			     6.274133333333333e-5 * t23 + t1210 * 
			    1.103598933333333e-6 * t227 - t1216 * 
			    3.0642176e-9 * t887) - t233 * .03109 * t518 + 
			    t237 * .001321010150222857 * t519 + t242 * 1. * 
			    t260 * t518 - t33 * .03109 * t45 * (t211 * 
			    -1.25264 * t48 + t1207 * 1.812618666666667 * t52 
			    - t1210 * .6307285333333333 * t269 + t1216 * 
			    .06366208 * t900) + t380 * t528 + t190 * (t211 * 
			    -.0059768 * t197 + t1207 * 2.564736e-4 * t202 - 
			    t851 * 5.509872e-6 * t1207 + t384 * 4.770864e-4 * 
			    t211 + t858 * 2.38261824e-8 * t1207 - t389 * 
			    2.6473536e-6 * t211);
		    v2rhoasigmaab[i__] = 0.;
		    t1260 = t212 * t64;
		    v2rhoasigmabb[i__] = t380 * t558 + t190 * (t212 * 
			    2.564736e-4 * t202 * t64 - t851 * 5.509872e-6 * 
			    t1260 + t858 * 2.38261824e-8 * t1260);
		    t1268 = t202 * t9;
		    t1271 = t399 * t9;
		    v2rhobsigmaaa[i__] = t490 * t528 + t190 * (t399 * 
			    2.564736e-4 * t1268 - t851 * 5.509872e-6 * t1271 
			    + t858 * 2.38261824e-8 * t1271);
		    v2rhobsigmaab[i__] = 0.;
		    t1282 = sigmabb * t404;
		    t1285 = t71 * t411;
		    t1291 = t408 / t62 / t409 / t396;
		    v2rhobsigmabb[i__] = t59 * -1.2407009817988 * t539 - t60 *
			     .9305257363491 * (t398 * -.0054112 * t68 - t1282 
			    * 6.274133333333333e-5 * t78 + t1285 * 
			    1.103598933333333e-6 * t414 - t1291 * 
			    3.0642176e-9 * t973) - t420 * .03109 * t548 + 
			    t424 * .001321010150222857 * t549 + t429 * 1. * 
			    t447 * t548 - t88 * .03109 * t100 * (t398 * 
			    -1.25264 * t103 + t1282 * 1.812618666666667 * 
			    t107 - t1285 * .6307285333333333 * t456 + t1291 * 
			    .06366208 * t1130) + t490 * t558 + t190 * (t398 * 
			    -.0059768 * t197 + t1282 * 2.564736e-4 * t202 - 
			    t851 * 5.509872e-6 * t1282 + t384 * 4.770864e-4 * 
			    t398 + t858 * 2.38261824e-8 * t1282 - t389 * 
			    2.6473536e-6 * t398);
		    t1333 = sigmaaa * t505;
		    t1338 = t16 / t7 / t876;
		    v2sigmaaa2[i__] = t5 * -.9305257363491 * (t20 * 7.7056e-6 
			    * t23 - t1333 * 3.180928e-7 * t227 + t1338 * 
			    1.1490816e-9 * t887) - t33 * .03109 * t45 * (t20 *
			     -.38684 * t52 + t1333 * .1967344 * t269 - t1338 *
			     .02387328 * t900) + t190 * (t20 * -9.61776e-5 * 
			    t202 + t851 * 2.066202e-6 * t20 - t858 * 
			    8.9348184e-9 * t20);
		    v2sigmaaaab[i__] = 0.;
		    t1364 = t9 * t64;
		    v2sigmaaabb[i__] = t190 * (t1268 * -9.61776e-5 * t64 + 
			    t851 * 2.066202e-6 * t1364 - t858 * 8.9348184e-9 *
			     t1364);
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t1372 = sigmabb * t535;
		    t1377 = t71 / t62 / t962;
		    v2sigmabb2[i__] = t60 * -.9305257363491 * (t75 * 
			    7.7056e-6 * t78 - t1372 * 3.180928e-7 * t414 + 
			    t1377 * 1.1490816e-9 * t973) - t88 * .03109 * 
			    t100 * (t75 * -.38684 * t107 + t1372 * .1967344 * 
			    t456 - t1377 * .02387328 * t1130) + t190 * (t75 * 
			    -9.61776e-5 * t202 + t851 * 2.066202e-6 * t75 - 
			    t858 * 8.9348184e-9 * t75);
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
} /* uks_xc_b97__ */

/* Subroutine */ int rks_xc_b97__(integer *ideriv, integer *npt, doublereal *
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
    static doublereal s1, t2, t3, t4, t5, t7, t8, t10, t20, t11, t31, t14, 
	    t15, t32, t35, t27, t19, t28, t37, t43, t45, t49, t68, t75, t79, 
	    t18, t21, t24, t30, t39, t42, t46, t50, t53, t54, t58, t64, t67, 
	    t73, t76, t80, t83, t90, t96, t87, t89, t95, t99, t100, t200, 
	    t111, t103, t121, t105, t114, t115, t122, t118, t119, t125, t126, 
	    t128, t131, t134, t136, t137, t147, t158, t186, t193, t197, t102, 
	    t108, t120, t123, t124, t130, t133, t138, t139, t150, t151, t156, 
	    t159, t160, t165, t166, t168, t176, t179, t189, t196, t209, t210, 
	    t219, t224, t230, t233, t235, t238, t242, t243, t244, t252, t254, 
	    t255, t262, t264, t271, t280, t281, t290, t291, t293, t295, t297, 
	    rho, t301, t303, t307, t309, t311, t321, t322, t333, t335, t340, 
	    t342, t344, t345, t348, t349, t353, t365, t371, t380, t381, t383, 
	    t384, t410, t412, t414, t415, t417, t433, t436, t442, t471, t473, 
	    t476, t485, t490, sigma;


/*     A.D. Becke */
/*     Density-functional thermochemistry. V. Systematic optimization of */
/*     exchange-correlation functionals */
/*     J. Chem. Phys. 107 (1997) 8554-8560 */


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
		t8 = sigma / t5 / t4;
		t10 = t8 * .006349604207872798 + 1.;
/* Computing 2nd power */
		d__1 = sigma;
		t14 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t4;
		t15 = d__1 * d__1;
		t19 = t14 / t2 / t15 / rho;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t27 = 1 / rho;
		t28 = pow_dd(&t27, &c_b2);
		t31 = rho * (t28 * .1606016560364007 + 1.);
		t32 = pow_dd(&t27, &c_b4);
		t35 = sqrt(t27);
/* Computing 2nd power */
		d__1 = t28;
		t37 = d__1 * d__1;
		t43 = log(32.1646831778707 / (t32 * 12.48219874679732 + t28 * 
			4.844076716063854 + t35 * 2.326004811900819 + t37 * 
			.3819082618690966) + 1.);
		t45 = t8 * .3174802103936399 + 1.;
/* Computing 2nd power */
		d__1 = t45;
		t49 = d__1 * d__1;
		t68 = log(16.0818243221511 / (t32 * 5.98255043577108 + t28 * 
			2.225569421150687 + t35 * .8004286349993634 + t37 * 
			.1897004325747559) + 1.);
		t75 = t8 * .009524406311809197 + 1.;
/* Computing 2nd power */
		d__1 = t75;
		t79 = d__1 * d__1;
		zk[i__] = t2 * -.7385587663820224 * rho * (t8 * 
			.00322115421465387 / t10 + .8094 + t19 * 
			3.016150199764335e-5 / t20) - t31 * .03109 * t43 * (
			t8 * .745665770151542 / t45 + .1737 - t19 * 
			.2506537333502856 / t49) + (rho * -.062182 * (t28 * 
			.1325688999052018 + 1.) * t68 + t31 * .03109 * t43) * 
			(t8 * .007115683955552651 / t75 + .9454 - t19 * 
			4.169320658943715e-4 / t79);
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
		t7 = 1 / t5 / t4;
		t8 = sigma * t7;
		t10 = t8 * .006349604207872798 + 1.;
		t11 = 1 / t10;
/* Computing 2nd power */
		d__1 = sigma;
		t14 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t4;
		t15 = d__1 * d__1;
		t18 = 1 / t2 / t15 / rho;
		t19 = t14 * t18;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		t24 = t8 * .00322115421465387 * t11 + .8094 + t19 * 
			3.016150199764335e-5 * t21;
		t27 = 1 / rho;
		t28 = pow_dd(&t27, &c_b2);
		t30 = t28 * .1606016560364007 + 1.;
		t31 = rho * t30;
		t32 = pow_dd(&t27, &c_b4);
		t35 = sqrt(t27);
/* Computing 2nd power */
		d__1 = t28;
		t37 = d__1 * d__1;
		t39 = t32 * 12.48219874679732 + t28 * 4.844076716063854 + t35 
			* 2.326004811900819 + t37 * .3819082618690966;
		t42 = 32.1646831778707 / t39 + 1.;
		t43 = log(t42);
		t45 = t8 * .3174802103936399 + 1.;
		t46 = 1 / t45;
/* Computing 2nd power */
		d__1 = t45;
		t49 = d__1 * d__1;
		t50 = 1 / t49;
		t53 = t8 * .745665770151542 * t46 + .1737 - t19 * 
			.2506537333502856 * t50;
		t54 = t43 * t53;
		t58 = t28 * .1325688999052018 + 1.;
		t64 = t32 * 5.98255043577108 + t28 * 2.225569421150687 + t35 *
			 .8004286349993634 + t37 * .1897004325747559;
		t67 = 16.0818243221511 / t64 + 1.;
		t68 = log(t67);
		t73 = rho * -.062182 * t58 * t68 + t31 * .03109 * t43;
		t75 = t8 * .009524406311809197 + 1.;
		t76 = 1 / t75;
/* Computing 2nd power */
		d__1 = t75;
		t79 = d__1 * d__1;
		t80 = 1 / t79;
		t83 = t8 * .007115683955552651 * t76 + .9454 - t19 * 
			4.169320658943715e-4 * t80;
		zk[i__] = t3 * -.7385587663820224 * t24 - t31 * .03109 * t54 
			+ t73 * t83;
		t90 = sigma / t5 / t4 / rho;
		t96 = t14 / t2 / t15 / t4;
/* Computing 2nd power */
		d__1 = t15;
		t100 = d__1 * d__1;
		t103 = t14 * sigma / t100 / rho;
		t105 = 1 / t20 / t10;
		t111 = t30 * t43;
		t114 = 1 / t37;
		t115 = t27 * t114;
/* Computing 2nd power */
		d__1 = t39;
		t118 = d__1 * d__1;
		t119 = 1 / t118;
/* Computing 2nd power */
		d__1 = t32;
		t121 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t121;
		t122 = d__1 * d__1;
		t125 = 1 / t4;
		t126 = 1 / t122 / t32 * t125;
		t128 = t114 * t125;
		t131 = 1 / t35 * t125;
		t134 = 1 / t28 * t125;
		t136 = t126 * -4.160732915599108 - t128 * 3.229384477375903 - 
			t131 * 2.326004811900819 - t134 * .5092110158254621;
		t137 = 1 / t42;
		t147 = 1 / t49 / t45;
/* Computing 2nd power */
		d__1 = t64;
		t158 = d__1 * d__1;
		t186 = 1 / t79 / t75;
		vrhoa[i__] = t2 * -.9847450218426965 * t24 - t3 * 
			.3692793831910112 * (t90 * -.01717948914482064 * t11 
			- t96 * 2.126397314118042e-4 * t21 + t103 * 
			2.042811733333333e-6 * t105) - t111 * .03109 * t53 + 
			t115 * .001664368495390566 * t54 + t31 * .5 * t119 * 
			t136 * t137 * t53 - t31 * .015545 * t43 * (t90 * 
			-3.976884107474891 * t46 + t96 * 3.936221825555298 * 
			t50 - t103 * .8488277333333333 * t147) + (t58 * 
			-.062182 * t68 + rho * (t128 * .002747799777968419 * 
			t68 + t58 * 1. / t158 * (t126 * -.99709173929518 - 
			t128 * .7418564737168958 - t131 * .4002143174996817 - 
			t134 * .1264669550498372) / t67) + t111 * .03109 - 
			t115 * .001664368495390566 * t43 - t31 * .5 * t119 * 
			t136 * t137) * t83 + t73 * (t90 * -.01897515721480707 
			* t76 + t96 * .002404364791914262 * t80 - t103 * 
			2.11788288e-5 * t186);
		t193 = sigma * t18;
		t197 = t14 / t100;
		vsigmaaa[i__] = t3 * -.7385587663820224 * (t7 * 
			.01288461685861548 * t11 + t193 * 
			1.594797985588531e-4 * t21 - t197 * 1.5321088e-6 * 
			t105) - t31 * .03109 * t43 * (t7 * 2.982663080606168 *
			 t46 - t193 * 2.952166369166474 * t50 + t197 * 
			.6366208 * t147) + t73 * 2. * (t7 * .0142313679111053 
			* t76 - t193 * .001803273593935696 * t80 + t197 * 
			1.58841216e-5 * t186);
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
		t7 = 1 / t5 / t4;
		t8 = sigma * t7;
		t10 = t8 * .006349604207872798 + 1.;
		t11 = 1 / t10;
/* Computing 2nd power */
		d__1 = sigma;
		t14 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t4;
		t15 = d__1 * d__1;
		t18 = 1 / t2 / t15 / rho;
		t19 = t14 * t18;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		t24 = t8 * .00322115421465387 * t11 + .8094 + t19 * 
			3.016150199764335e-5 * t21;
		t27 = 1 / rho;
		t28 = pow_dd(&t27, &c_b2);
		t30 = t28 * .1606016560364007 + 1.;
		t31 = rho * t30;
		t32 = pow_dd(&t27, &c_b4);
		t35 = sqrt(t27);
/* Computing 2nd power */
		d__1 = t28;
		t37 = d__1 * d__1;
		t39 = t32 * 12.48219874679732 + t28 * 4.844076716063854 + t35 
			* 2.326004811900819 + t37 * .3819082618690966;
		t42 = 32.1646831778707 / t39 + 1.;
		t43 = log(t42);
		t45 = t8 * .3174802103936399 + 1.;
		t46 = 1 / t45;
/* Computing 2nd power */
		d__1 = t45;
		t49 = d__1 * d__1;
		t50 = 1 / t49;
		t53 = t8 * .745665770151542 * t46 + .1737 - t19 * 
			.2506537333502856 * t50;
		t54 = t43 * t53;
		t58 = t28 * .1325688999052018 + 1.;
		t64 = t32 * 5.98255043577108 + t28 * 2.225569421150687 + t35 *
			 .8004286349993634 + t37 * .1897004325747559;
		t67 = 16.0818243221511 / t64 + 1.;
		t68 = log(t67);
		t73 = rho * -.062182 * t58 * t68 + t31 * .03109 * t43;
		t75 = t8 * .009524406311809197 + 1.;
		t76 = 1 / t75;
/* Computing 2nd power */
		d__1 = t75;
		t79 = d__1 * d__1;
		t80 = 1 / t79;
		t83 = t8 * .007115683955552651 * t76 + .9454 - t19 * 
			4.169320658943715e-4 * t80;
		zk[i__] = t3 * -.7385587663820224 * t24 - t31 * .03109 * t54 
			+ t73 * t83;
		t87 = t4 * rho;
		t89 = 1 / t5 / t87;
		t90 = sigma * t89;
		t95 = 1 / t2 / t15 / t4;
		t96 = t14 * t95;
		t99 = t14 * sigma;
/* Computing 2nd power */
		d__1 = t15;
		t100 = d__1 * d__1;
		t102 = 1 / t100 / rho;
		t103 = t99 * t102;
		t105 = 1 / t20 / t10;
		t108 = t90 * -.01717948914482064 * t11 - t96 * 
			2.126397314118042e-4 * t21 + t103 * 
			2.042811733333333e-6 * t105;
		t111 = t30 * t43;
		t114 = 1 / t37;
		t115 = t27 * t114;
/* Computing 2nd power */
		d__1 = t39;
		t118 = d__1 * d__1;
		t119 = 1 / t118;
		t120 = t31 * t119;
/* Computing 2nd power */
		d__1 = t32;
		t121 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t121;
		t122 = d__1 * d__1;
		t123 = t122 * t32;
		t124 = 1 / t123;
		t125 = 1 / t4;
		t126 = t124 * t125;
		t128 = t114 * t125;
		t130 = 1 / t35;
		t131 = t130 * t125;
		t133 = 1 / t28;
		t134 = t133 * t125;
		t136 = t126 * -4.160732915599108 - t128 * 3.229384477375903 - 
			t131 * 2.326004811900819 - t134 * .5092110158254621;
		t137 = 1 / t42;
		t138 = t136 * t137;
		t139 = t138 * t53;
		t147 = 1 / t49 / t45;
		t150 = t90 * -3.976884107474891 * t46 + t96 * 
			3.936221825555298 * t50 - t103 * .8488277333333333 * 
			t147;
		t151 = t43 * t150;
		t156 = t128 * t68;
/* Computing 2nd power */
		d__1 = t64;
		t158 = d__1 * d__1;
		t159 = 1 / t158;
		t160 = t58 * t159;
		t165 = t126 * -.99709173929518 - t128 * .7418564737168958 - 
			t131 * .4002143174996817 - t134 * .1264669550498372;
		t166 = 1 / t67;
		t168 = t160 * t165 * t166;
		t176 = t119 * t136 * t137;
		t179 = t58 * -.062182 * t68 + rho * (t156 * 
			.002747799777968419 + t168 * 1.) + t111 * .03109 - 
			t115 * .001664368495390566 * t43 - t31 * .5 * t176;
		t186 = 1 / t79 / t75;
		t189 = t90 * -.01897515721480707 * t76 + t96 * 
			.002404364791914262 * t80 - t103 * 2.11788288e-5 * 
			t186;
		vrhoa[i__] = t2 * -.9847450218426965 * t24 - t3 * 
			.3692793831910112 * t108 - t111 * .03109 * t53 + t115 
			* .001664368495390566 * t54 + t120 * .5 * t139 - t31 *
			 .015545 * t151 + t179 * t83 + t73 * t189;
		t193 = sigma * t18;
		t196 = 1 / t100;
		t197 = t14 * t196;
		t200 = t7 * .01288461685861548 * t11 + t193 * 
			1.594797985588531e-4 * t21 - t197 * 1.5321088e-6 * 
			t105;
		t209 = t7 * 2.982663080606168 * t46 - t193 * 
			2.952166369166474 * t50 + t197 * .6366208 * t147;
		t210 = t43 * t209;
		t219 = t7 * .0142313679111053 * t76 - t193 * 
			.001803273593935696 * t80 + t197 * 1.58841216e-5 * 
			t186;
		vsigmaaa[i__] = t3 * -.7385587663820224 * t200 - t31 * .03109 
			* t210 + t73 * 2. * t219;
		t224 = sigma / t5 / t15;
		t230 = t14 / t2 / t15 / t87;
		t233 = t100 * t4;
		t235 = t99 / t233;
/* Computing 2nd power */
		d__1 = t14;
		t238 = d__1 * d__1;
		t242 = t238 / t5 / t100 / t15;
/* Computing 2nd power */
		d__1 = t20;
		t243 = d__1 * d__1;
		t244 = 1 / t243;
		t252 = 1 / t87;
		t254 = 1 / t37 / t27;
		t255 = t252 * t254;
		t262 = 1 / t118 / t39;
/* Computing 2nd power */
		d__1 = t136;
		t264 = d__1 * d__1;
		t271 = t30 * t119;
/* Computing 2nd power */
		d__1 = t49;
		t280 = d__1 * d__1;
		t281 = 1 / t280;
		t290 = 1 / t15;
		t291 = 1 / t123 / t27 * t290;
		t293 = t124 * t252;
		t295 = t254 * t290;
		t297 = t114 * t252;
		t301 = 1 / t35 / t27 * t290;
		t303 = t130 * t252;
		t307 = 1 / t28 / t27 * t290;
		t309 = t133 * t252;
		t311 = t291 * -6.934554859331846 + t293 * 16.64293166239643 - 
			t295 * 4.305845969834537 + t297 * 12.91753790950361 - 
			t301 * 2.326004811900819 + t303 * 9.304019247603276 - 
			t307 * .3394740105503081 + t309 * 2.036844063301848;
		t321 = t156 * .005495599555936838;
		t322 = t168 * 2.;
		t333 = log(29.60857464321668 / (t32 * 8.157414703487641 + t28 
			* 2.247591863577616 + t35 * .4300972471276643 + t37 * 
			.1911512595127338) + 1.);
		t335 = (t28 * .06901399211255825 + 1.) * t333 * t125;
		t340 = t128 * .08837926660346786 * t159 * t165 * t166;
		t342 = t297 * .005495599555936838 * t68;
		t344 = t295 * .001831866518645613 * t68;
/* Computing 2nd power */
		d__1 = t158;
		t345 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t165;
		t348 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t67;
		t349 = d__1 * d__1;
		t353 = t58 * 16.0818243221511 / t345 * t348 / t349;
		t365 = t160 * 1. * (t291 * -.8309097827459833 + t293 * 
			1.99418347859036 - t295 * .4945709824779306 + t297 * 
			1.483712947433792 - t301 * .2001071587498409 + t303 * 
			.8004286349993634 - t307 * .04215565168327908 + t309 *
			 .2529339100996745) * t166;
		t371 = t58 * 2. / t158 / t64 * t348 * t166;
/* Computing 2nd power */
		d__1 = t118;
		t380 = d__1 * d__1;
		t381 = 1 / t380;
/* Computing 2nd power */
		d__1 = t42;
		t383 = d__1 * d__1;
		t384 = 1 / t383;
		t410 = t230 * t80;
		t412 = t235 * t186;
/* Computing 2nd power */
		d__1 = t79;
		t414 = d__1 * d__1;
		t415 = 1 / t414;
		t417 = t242 * 1.613726165595571e-6 * t415;
		s1 = t3 * -.3692793831910112 * (t224 * .1259829203953514 * 
			t11 + t230 * .002111660829546542 * t21 - t235 * 
			5.117251128888889e-5 * t105 + t242 * 
			2.075367356458441e-7 * t244) - t111 * .06218 * t150 + 
			t255 * .002219157993854088 * t54 - t115 * 
			.1070677706909338 * t119 * t139 - t31 * 1. * t262 * 
			t264 * t137 * t53 + t115 * .003328736990781132 * t151 
			+ t271 * 2. * t139 - t31 * .015545 * t43 * (t224 * 
			29.1638167881492 * t46 - t230 * 56.59258047384578 * 
			t50 + t235 * 28.60873955555556 * t147 - t242 * 
			4.31177611786597 * t281) + t120 * .5 * t311 * t137 * 
			t53;
		v2rhoa2[i__] = s1 - .6564966812284644 / t5 * t24 - t2 * 
			1.969490043685393 * t108 + (t321 + t322 + rho * (t335 
			* .03377399869956914 - t340 - t342 + t344 + t353 + 
			t365 - t371) - t255 * .002219157993854088 * t43 + 
			t115 * .1070677706909338 * t176 - t271 * 2. * t138 - 
			t31 * 16.08234158893535 * t381 * t264 * t384 - t31 * 
			.5 * t119 * t311 * t137 + t31 * 1. * t262 * t264 * 
			t137) * t83 + t31 * 16.08234158893535 * t381 * t264 * 
			t384 * t53 + t120 * 1. * t138 * t150 + t179 * 4. * 
			t189 + t73 * (t224 * .1391511529085852 * t76 - t410 * 
			.02452558687152736 + t412 * 3.903992832e-4 - t417) + (
			t321 + t322 + rho * (t335 * -.03377399869956914 - 
			t340 - t342 + t344 + t353 + t365 - t371)) * t83 + t73 
			* (t410 * -.006893578397489445 + t412 * 2.35087872e-4 
			- t417);
		t433 = sigma * t95;
		t436 = t14 * t102;
		t442 = t99 / t5 / t100 / t87;
		t471 = t433 * t80;
		t473 = t436 * t186;
		t476 = t442 * 1.210294624196678e-6 * t415;
		v2rhoasigmaaa[i__] = t2 * -.9847450218426965 * t200 - t3 * 
			.3692793831910112 * (t89 * -.06871795657928257 * t11 
			- t433 * .001264786025042201 * t21 + t436 * 
			3.531516586666667e-5 * t105 - t442 * 
			1.556525517343831e-7 * t244) - t111 * .03109 * t209 + 
			t115 * .001664368495390566 * t210 + t120 * .5 * t138 *
			 t209 - t31 * .015545 * t43 * (t89 * 
			-15.90753642989956 * t46 + t433 * 36.54010261705139 * 
			t50 - t436 * 20.18331306666667 * t147 + t442 * 
			3.233832088399478 * t281) + t179 * 2. * t219 + t73 * (
			t89 * -.07590062885922828 * t76 + t471 * 
			.01478764296577413 - t473 * 2.610312192e-4 + t476) + 
			t73 * (t471 * .005170183798117084 - t473 * 
			1.76315904e-4 + t476);
		t485 = sigma * t196;
		t490 = t14 / t5 / t233;
		v2sigmaaa2[i__] = t3 * -.7385587663820224 * (t18 * 
			3.106703245462379e-4 * t21 - t485 * 2.03579392e-5 * 
			t105 + t490 * 1.167394138007873e-7 * t244) - t31 * 
			.03109 * t43 * (t18 * -15.59641148612265 * t50 + t485 
			* 12.5910016 * t147 - t490 * 2.425374066299608 * t281)
			 + t73 * 4. * (t18 * -.003877637848587813 * t80 + 
			t485 * 1.32236928e-4 * t186 - t490 * 
			9.077209681475088e-7 * t415);
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
} /* rks_xc_b97__ */


/* xc_xc_hcth120.f -- translated by f2c (version 20050501).
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

/* :XC_HCTH120subrstart */
/*    Generated: Sat May  8 18:39:42 GMT 2004 */
/* Subroutine */ int uks_xc_hcth120__(integer *ideriv, integer *npt, 
	doublereal *rhoa1, doublereal *rhob1, doublereal *sigmaaa1, 
	doublereal *sigmabb1, doublereal *sigmaab1, doublereal *zk, 
	doublereal *vrhoa, doublereal *vrhob, doublereal *vsigmaaa, 
	doublereal *vsigmabb, doublereal *vsigmaab, doublereal *v2rhoa2, 
	doublereal *v2rhob2, doublereal *v2rhoab, doublereal *v2rhoasigmaaa, 
	doublereal *v2rhoasigmaab, doublereal *v2rhoasigmabb, doublereal *
	v2rhobsigmabb, doublereal *v2rhobsigmaab, doublereal *v2rhobsigmaaa, 
	doublereal *v2sigmaaa2, doublereal *v2sigmaaaab, doublereal *
	v2sigmaaabb, doublereal *v2sigmaab2, doublereal *v2sigmaabbb, 
	doublereal *v2sigmabb2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t4, t5, t6, t7, t8, t9, t10, t20, t12, 
	    t21, t14, t15, t25, t32, t27, t19, t36, t37, t44, t45, t49, t52, 
	    t54, t60, t62, t66, t74, t16, t17, t22, t29, t34, t38, t39, t46, 
	    t47, t50, t51, t56, t64, t68, t76, t84, t86, t87, t90, t92, t96, 
	    t97, t11, t18, t24, t26, t35, t101, t102, t130, t131, t114, t142, 
	    t107, t126, t109, t118, t119, t127, t134, t136, t144, t148, t156, 
	    t164, t165, t166, t169, t172, t174, t180, t182, t193, t196, t197, 
	    t198, t199, t202, t203, t205, t206, t207, t208, t209, t211, t227, 
	    t243, t246, t250, t251, t260, t261, t41, t48, t59, t63, t67, t71, 
	    t75, t78, t79, t93, t98, t104, t111, t113, t122, t129, t133, t159,
	     t168, t171, t13, t23, rho, t28, t31, t40, t43, t58, t61, t65, 
	    t69, t73, t77, t80, t81, t85, t89, t100, t103, t106, t108, t117, 
	    t120, t123, t138, t141, t145, t149, t153, t157, t160, t161, t176, 
	    t179, t184, t189, t192, t194, t210, t213, t216, t218, t223, t226, 
	    t230, t231, t233, t240, t247, t252, t255, t257, t262, t265, t269, 
	    t272, t275, t278, t283, t289, t296, t298, t304, t307, t308, t311, 
	    t312, t314, t315, t318, t329, t330, t344, t352, t353, t354, t355, 
	    t356, t359, t360, t363, t367, t370, t375, t376, t380, t381, t394, 
	    t395, t396, t400, t404, t409, t410, t413, t414, rhoa, rhob, t421, 
	    t438, t442, t445, t459, t462, t465, t470, t477, t480, t483, t486, 
	    t491, t497, t504, t506, t512, t515, t516, t519, t520, t522, t523, 
	    t526, t537, t538, t552, t565, t602, t605, t608, t614, t648, t651, 
	    t654, t660, t33, t116, t128, t132, t146, t147, t162, t163, t217, 
	    t229, t236, t238, t249, t263, t270, t299, t319, t328, t331, t334, 
	    t338, t115, t214, t271, t277, t282, t288, t292, t295, t301, t313, 
	    t316, t317, t323, t326, t332, t347, t348, t351, t357, t358, t361, 
	    t362, t366, t369, t372, t373, t377, t378, t379, t382, t383, t388, 
	    t389, t390, t391, t392, t393, t397, t401, t406, t407, sigma, t408,
	     t412, t415, t417, t418, t422, t423, t428, t429, t435, t436, t437,
	     t439, t440, t441, t444, t452, t455, t469, t473, t479, t485, t490,
	     t496, t500, t503, t509, t521, t524, t525, t531, t534, t539, t540,
	     t555, t556, t559, t562, t567, t568, t569, t571, t572, t573, t574,
	     t575, t583, t586, t598, t613, t617, t630, t631, t644, t659, t663,
	     t676, t677, t690, t696, t699, t702, t705, t710, t715, t722, t739,
	     t745, t751, t755, t758, t760, t766, t774, t776, t784, t786, t787,
	     t792, t813, t821, t822, t824, t825, t841, t848, t850, t852, t853,
	     t854, t856, t876, t877, t879, t883, t884, t885, t889, t892, t894,
	     t895, t896, t897, t898, t901, t902, t903, t904, t905, t906, t907,
	     t910, t912, t913, t916, t917, t918, t924, t928, t929, t932, t933,
	     t940, t943, t945, t949, t954, t955, t957, t959, t963, t965, t971,
	     t973, t977, t979, t984, t991, t993, t994, t1005, t1006, t1009, 
	    t1010, t1011, t1014, t1017, t1018, t1030, t1046, t1049, t1060, 
	    t1063, t1065, t1068, t1071, t1074, t1082, t1084, t1085, t1088, 
	    t1090, t1104, t1107, t1110, t1114, t1120, t1126, t1130, t1133, 
	    t1135, t1141, t1143, t1144, t1151, t1153, t1169, t1178, t1199, 
	    t1207, t1208, t1210, t1211, t1222, t1225, t1228, sigmaaa, sigmaab,
	     sigmabb, t1234, t1245, t1247, t1251, t1255, t1258, t1261, t1266, 
	    t1268, t1270, t1291, t1330, t1346, t1362, t1365, t1375, t1393, 
	    t1396, t1399, t1402, t1407, t1463, t1475, t1478, t1493, t1496, 
	    t1499, t1502, t1507, t1561, t1564, t1567, t1571, t1605, t1617, 
	    t1620, t1623, t1627;


/*     A.D. Boese, N.L. Doltsinis, N.C. Handy, M. Sprik */
/*     New generalized gradient approximation functionals */
/*     J. Chem. Phys. 112 (2000) 1670-1678. */


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
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t27 = t14 * sigmabb / t25;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t36 = t32 / t5 / t25 / t4;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t44 = 1 / rhob;
		    t45 = pow_dd(&t44, &c_b2);
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t60 = log(32.1646831778707 / (t49 * 11.12037486309468 + 
			    t45 * 3.844746237447211 + t52 * 1.644733775567609 
			    + t54 * .2405871291288192) + 1.);
		    t62 = t8 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhob * (1.09163 - t8 * 
			    .00298886 / t10 + t19 * 8.125328e-5 / t20 - t27 * 
			    2.6287744e-7 / t20 / t10 + t36 * 2.9996288e-10 / 
			    t37) - rhob * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (.489508 - t8 * .0521398 / t62 + t19 *
			     .01731668 / t66 - t27 * .01593976 / t66 / t62 + 
			    t36 * .003976496 / t74);
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
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t27 = t14 * sigmaaa / t25;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t36 = t32 / t5 / t25 / t4;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t44 = 1 / rhoa;
		    t45 = pow_dd(&t44, &c_b2);
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t60 = log(32.1646831778707 / (t49 * 11.12037486309468 + 
			    t45 * 3.844746237447211 + t52 * 1.644733775567609 
			    + t54 * .2405871291288192) + 1.);
		    t62 = t8 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    zk[i__] = t2 * -.9305257363491 * rhoa * (1.09163 - t8 * 
			    .00298886 / t10 + t19 * 8.125328e-5 / t20 - t27 * 
			    2.6287744e-7 / t20 / t10 + t36 * 2.9996288e-10 / 
			    t37) - rhoa * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (.489508 - t8 * .0521398 / t62 + t19 *
			     .01731668 / t66 - t27 * .01593976 / t66 / t62 + 
			    t36 * .003976496 / t74);
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
/* Computing 2nd power */
		    d__1 = t17;
		    t27 = d__1 * d__1;
		    t29 = t16 * sigmaaa / t27;
/* Computing 2nd power */
		    d__1 = t16;
		    t34 = d__1 * d__1;
		    t38 = t34 / t7 / t27 / t6;
/* Computing 2nd power */
		    d__1 = t22;
		    t39 = d__1 * d__1;
		    t46 = 1 / rhoa;
		    t47 = pow_dd(&t46, &c_b2);
		    t50 = rhoa * (t47 * .1274696188700087 + 1.);
		    t51 = pow_dd(&t46, &c_b4);
		    t54 = sqrt(t46);
/* Computing 2nd power */
		    d__1 = t47;
		    t56 = d__1 * d__1;
		    t62 = log(32.1646831778707 / (t51 * 11.12037486309468 + 
			    t47 * 3.844746237447211 + t54 * 1.644733775567609 
			    + t56 * .2405871291288192) + 1.);
		    t64 = t10 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t64;
		    t68 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t68;
		    t76 = d__1 * d__1;
		    t84 = pow_dd(&rhob, &c_b2);
/* Computing 2nd power */
		    d__1 = rhob;
		    t86 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t84;
		    t87 = d__1 * d__1;
		    t90 = sigmabb / t87 / t86;
		    t92 = t90 * .004 + 1.;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t96 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t86;
		    t97 = d__1 * d__1;
		    t101 = t96 / t84 / t97 / rhob;
/* Computing 2nd power */
		    d__1 = t92;
		    t102 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t97;
		    t107 = d__1 * d__1;
		    t109 = t96 * sigmabb / t107;
/* Computing 2nd power */
		    d__1 = t96;
		    t114 = d__1 * d__1;
		    t118 = t114 / t87 / t107 / t86;
/* Computing 2nd power */
		    d__1 = t102;
		    t119 = d__1 * d__1;
		    t126 = 1 / rhob;
		    t127 = pow_dd(&t126, &c_b2);
		    t130 = rhob * (t127 * .1274696188700087 + 1.);
		    t131 = pow_dd(&t126, &c_b4);
		    t134 = sqrt(t126);
/* Computing 2nd power */
		    d__1 = t127;
		    t136 = d__1 * d__1;
		    t142 = log(32.1646831778707 / (t131 * 11.12037486309468 + 
			    t127 * 3.844746237447211 + t134 * 
			    1.644733775567609 + t136 * .2405871291288192) + 
			    1.);
		    t144 = t90 * .2 + 1.;
/* Computing 2nd power */
		    d__1 = t144;
		    t148 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t148;
		    t156 = d__1 * d__1;
		    t164 = rhoa + rhob;
		    t165 = 1 / t164;
		    t166 = pow_dd(&t165, &c_b2);
		    t169 = pow_dd(&t165, &c_b4);
		    t172 = sqrt(t165);
/* Computing 2nd power */
		    d__1 = t166;
		    t174 = d__1 * d__1;
		    t180 = log(16.0818243221511 / (t169 * 5.98255043577108 + 
			    t166 * 2.225569421150687 + t172 * 
			    .8004286349993634 + t174 * .1897004325747559) + 
			    1.);
		    t182 = (t166 * .1325688999052018 + 1.) * .062182 * t180;
		    t193 = log(29.60857464321668 / (t169 * 8.157414703487641 
			    + t166 * 2.247591863577616 + t172 * 
			    .4300972471276643 + t174 * .1911512595127338) + 
			    1.);
		    t196 = rhoa - rhob * 1.;
		    t197 = t196 * t165;
		    t198 = t197 + 1.;
		    t199 = pow_dd(&t198, &c_b2);
		    t202 = 1. - t197 * 1.;
		    t203 = pow_dd(&t202, &c_b2);
		    t205 = t199 * t198 + t203 * t202 - 2.;
/* Computing 2nd power */
		    d__1 = t196;
		    t206 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t206;
		    t207 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t164;
		    t208 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t208;
		    t209 = d__1 * d__1;
		    t211 = t207 / t209;
		    t227 = log(32.1646831778707 / (t169 * 11.12037486309468 + 
			    t166 * 3.844746237447211 + t172 * 
			    1.644733775567609 + t174 * .2405871291288192) + 
			    1.);
		    t243 = t10 * .5 + t90 * .5;
		    t246 = t10 * .003 + 1. + t90 * .003;
/* Computing 2nd power */
		    d__1 = t243;
		    t250 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t246;
		    t251 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t250;
		    t260 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t251;
		    t261 = d__1 * d__1;
		    s1 = t4 * -.9305257363491 * rhoa * (1.09163 - t10 * 
			    .00298886 / t12 + t21 * 8.125328e-5 / t22 - t29 * 
			    2.6287744e-7 / t22 / t12 + t38 * 2.9996288e-10 / 
			    t39) - t50 * .03109 * t62 * (.489508 - t10 * 
			    .0521398 / t64 + t21 * .01731668 / t68 - t29 * 
			    .01593976 / t68 / t64 + t38 * .003976496 / t76);
		    s2 = s1 - t84 * .9305257363491 * rhob * (1.09163 - t90 * 
			    .00298886 / t92 + t101 * 8.125328e-5 / t102 - 
			    t109 * 2.6287744e-7 / t102 / t92 + t118 * 
			    2.9996288e-10 / t119);
		    zk[i__] = s2 - t130 * .03109 * t142 * (.489508 - t90 * 
			    .0521398 / t144 + t101 * .01731668 / t148 - t109 *
			     .01593976 / t148 / t144 + t118 * .003976496 / 
			    t156) + (t164 * (-t182 + (t166 * 
			    .06901399211255825 + 1.) * .03799574853701528 * 
			    t193 * t205 * (1. - t211 * 1.) + ((t166 * 
			    .1274696188700087 + 1.) * -.03109 * t227 + t182) *
			     1.923661050931536 * t205 * t211) + t50 * .03109 *
			     t62 + t130 * .03109 * t142) * (t243 * .04157892 /
			     t246 + .51473 - t250 * 8.894628e-4 / t251 + t250 
			    * 4.9917168e-6 * t243 / t251 / t246 - t260 * 
			    1.46751264e-8 / t261);
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
		    t16 = t15 * rhob;
		    t18 = 1 / t2 / t16;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t14 * sigmabb;
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t26 = 1 / t25;
		    t27 = t24 * t26;
		    t29 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t35 = 1 / t5 / t25 / t4;
		    t36 = t32 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t38 = 1 / t37;
		    t41 = 1.09163 - t8 * .00298886 * t11 + t19 * 8.125328e-5 *
			     t21 - t27 * 2.6287744e-7 * t29 + t36 * 
			    2.9996288e-10 * t38;
		    t44 = 1 / rhob;
		    t45 = pow_dd(&t44, &c_b2);
		    t47 = t45 * .1274696188700087 + 1.;
		    t48 = rhob * t47;
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t56 = t49 * 11.12037486309468 + t45 * 3.844746237447211 + 
			    t52 * 1.644733775567609 + t54 * .2405871291288192;
		    t59 = 32.1646831778707 / t56 + 1.;
		    t60 = log(t59);
		    t62 = t8 * .2 + 1.;
		    t63 = 1 / t62;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t71 = 1 / t66 / t62;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t78 = .489508 - t8 * .0521398 * t63 + t19 * .01731668 * 
			    t67 - t27 * .01593976 * t71 + t36 * .003976496 * 
			    t75;
		    t79 = t60 * t78;
		    zk[i__] = t3 * -.9305257363491 * t41 - t48 * .03109 * t79;
		    vrhoa[i__] = 0.;
		    t84 = t4 * rhob;
		    t87 = sigmabb / t5 / t84;
		    t90 = t15 * t4;
		    t93 = t14 / t2 / t90;
		    t98 = t24 / t25 / rhob;
		    t104 = t32 / t5 / t25 / t84;
		    t111 = t32 * sigmabb / t2 / t25 / t90;
		    t113 = 1 / t37 / t10;
		    t122 = 1 / t54;
/* Computing 2nd power */
		    d__1 = t56;
		    t126 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t129 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t129;
		    t130 = d__1 * d__1;
		    t133 = 1 / t4;
		    t159 = 1 / t74 / t62;
		    vrhob[i__] = t2 * -1.2407009817988 * t41 - t3 * 
			    .9305257363491 * (t87 * .007970293333333333 * t11 
			    - t93 * 4.65232e-4 * t21 + t98 * 
			    3.836422826666667e-6 * t29 - t104 * 
			    1.161168213333333e-8 * t38 + t111 * 
			    1.279841621333333e-11 * t113) - t47 * .03109 * 
			    t60 * t78 + t44 * .001321010150222857 * t122 * 
			    t79 + t48 * 1. / t126 * (-1.853395810515781 / 
			    t130 / t49 * t133 - t122 * 1.28158207914907 * 
			    t133 - .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * .1390394666666667 * t63 - 
			    t93 * .12016352 * t67 + t98 * .1459892053333333 * 
			    t71 - t104 * .06791957333333333 * t75 + t111 * 
			    .008483191466666667 * t159);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t168 = sigmabb * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.00298886 * 
			    t11 + t168 * 1.74462e-4 * t21 - t171 * 
			    1.43865856e-6 * t29 + t174 * 4.3543808e-9 * t38 - 
			    t180 * 4.79940608e-12 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.0521398 * t63 + t168 * .04506132 * 
			    t67 - t171 * .054745952 * t71 + t174 * .02546984 *
			     t75 - t180 * .0031811968 * t159);
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
		    t16 = t15 * rhoa;
		    t18 = 1 / t2 / t16;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t14 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t26 = 1 / t25;
		    t27 = t24 * t26;
		    t29 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t35 = 1 / t5 / t25 / t4;
		    t36 = t32 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t38 = 1 / t37;
		    t41 = 1.09163 - t8 * .00298886 * t11 + t19 * 8.125328e-5 *
			     t21 - t27 * 2.6287744e-7 * t29 + t36 * 
			    2.9996288e-10 * t38;
		    t44 = 1 / rhoa;
		    t45 = pow_dd(&t44, &c_b2);
		    t47 = t45 * .1274696188700087 + 1.;
		    t48 = rhoa * t47;
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t56 = t49 * 11.12037486309468 + t45 * 3.844746237447211 + 
			    t52 * 1.644733775567609 + t54 * .2405871291288192;
		    t59 = 32.1646831778707 / t56 + 1.;
		    t60 = log(t59);
		    t62 = t8 * .2 + 1.;
		    t63 = 1 / t62;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t71 = 1 / t66 / t62;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t78 = .489508 - t8 * .0521398 * t63 + t19 * .01731668 * 
			    t67 - t27 * .01593976 * t71 + t36 * .003976496 * 
			    t75;
		    t79 = t60 * t78;
		    zk[i__] = t3 * -.9305257363491 * t41 - t48 * .03109 * t79;
		    t84 = t4 * rhoa;
		    t87 = sigmaaa / t5 / t84;
		    t90 = t15 * t4;
		    t93 = t14 / t2 / t90;
		    t98 = t24 / t25 / rhoa;
		    t104 = t32 / t5 / t25 / t84;
		    t111 = t32 * sigmaaa / t2 / t25 / t90;
		    t113 = 1 / t37 / t10;
		    t122 = 1 / t54;
/* Computing 2nd power */
		    d__1 = t56;
		    t126 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t49;
		    t129 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t129;
		    t130 = d__1 * d__1;
		    t133 = 1 / t4;
		    t159 = 1 / t74 / t62;
		    vrhoa[i__] = t2 * -1.2407009817988 * t41 - t3 * 
			    .9305257363491 * (t87 * .007970293333333333 * t11 
			    - t93 * 4.65232e-4 * t21 + t98 * 
			    3.836422826666667e-6 * t29 - t104 * 
			    1.161168213333333e-8 * t38 + t111 * 
			    1.279841621333333e-11 * t113) - t47 * .03109 * 
			    t60 * t78 + t44 * .001321010150222857 * t122 * 
			    t79 + t48 * 1. / t126 * (-1.853395810515781 / 
			    t130 / t49 * t133 - t122 * 1.28158207914907 * 
			    t133 - .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * .1390394666666667 * t63 - 
			    t93 * .12016352 * t67 + t98 * .1459892053333333 * 
			    t71 - t104 * .06791957333333333 * t75 + t111 * 
			    .008483191466666667 * t159);
		    vrhob[i__] = 0.;
		    t168 = sigmaaa * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.00298886 * 
			    t11 + t168 * 1.74462e-4 * t21 - t171 * 
			    1.43865856e-6 * t29 + t174 * 4.3543808e-9 * t38 - 
			    t180 * 4.79940608e-12 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.0521398 * t63 + t168 * .04506132 * 
			    t67 - t171 * .054745952 * t71 + t174 * .02546984 *
			     t75 - t180 * .0031811968 * t159);
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
		    t18 = t17 * rhoa;
		    t20 = 1 / t4 / t18;
		    t21 = t16 * t20;
/* Computing 2nd power */
		    d__1 = t12;
		    t22 = d__1 * d__1;
		    t23 = 1 / t22;
		    t26 = t16 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t17;
		    t27 = d__1 * d__1;
		    t28 = 1 / t27;
		    t29 = t26 * t28;
		    t31 = 1 / t22 / t12;
/* Computing 2nd power */
		    d__1 = t16;
		    t34 = d__1 * d__1;
		    t37 = 1 / t7 / t27 / t6;
		    t38 = t34 * t37;
/* Computing 2nd power */
		    d__1 = t22;
		    t39 = d__1 * d__1;
		    t40 = 1 / t39;
		    t43 = 1.09163 - t10 * .00298886 * t13 + t21 * 8.125328e-5 
			    * t23 - t29 * 2.6287744e-7 * t31 + t38 * 
			    2.9996288e-10 * t40;
		    t46 = 1 / rhoa;
		    t47 = pow_dd(&t46, &c_b2);
		    t49 = t47 * .1274696188700087 + 1.;
		    t50 = rhoa * t49;
		    t51 = pow_dd(&t46, &c_b4);
		    t54 = sqrt(t46);
/* Computing 2nd power */
		    d__1 = t47;
		    t56 = d__1 * d__1;
		    t58 = t51 * 11.12037486309468 + t47 * 3.844746237447211 + 
			    t54 * 1.644733775567609 + t56 * .2405871291288192;
		    t61 = 32.1646831778707 / t58 + 1.;
		    t62 = log(t61);
		    t64 = t10 * .2 + 1.;
		    t65 = 1 / t64;
/* Computing 2nd power */
		    d__1 = t64;
		    t68 = d__1 * d__1;
		    t69 = 1 / t68;
		    t73 = 1 / t68 / t64;
/* Computing 2nd power */
		    d__1 = t68;
		    t76 = d__1 * d__1;
		    t77 = 1 / t76;
		    t80 = .489508 - t10 * .0521398 * t65 + t21 * .01731668 * 
			    t69 - t29 * .01593976 * t73 + t38 * .003976496 * 
			    t77;
		    t81 = t62 * t80;
		    t84 = pow_dd(&rhob, &c_b2);
		    t85 = t84 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t86 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t84;
		    t87 = d__1 * d__1;
		    t89 = 1 / t87 / t86;
		    t90 = sigmabb * t89;
		    t92 = t90 * .004 + 1.;
		    t93 = 1 / t92;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t96 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t86;
		    t97 = d__1 * d__1;
		    t98 = t97 * rhob;
		    t100 = 1 / t84 / t98;
		    t101 = t96 * t100;
/* Computing 2nd power */
		    d__1 = t92;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
		    t106 = t96 * sigmabb;
/* Computing 2nd power */
		    d__1 = t97;
		    t107 = d__1 * d__1;
		    t108 = 1 / t107;
		    t109 = t106 * t108;
		    t111 = 1 / t102 / t92;
/* Computing 2nd power */
		    d__1 = t96;
		    t114 = d__1 * d__1;
		    t117 = 1 / t87 / t107 / t86;
		    t118 = t114 * t117;
/* Computing 2nd power */
		    d__1 = t102;
		    t119 = d__1 * d__1;
		    t120 = 1 / t119;
		    t123 = 1.09163 - t90 * .00298886 * t93 + t101 * 
			    8.125328e-5 * t103 - t109 * 2.6287744e-7 * t111 + 
			    t118 * 2.9996288e-10 * t120;
		    t126 = 1 / rhob;
		    t127 = pow_dd(&t126, &c_b2);
		    t129 = t127 * .1274696188700087 + 1.;
		    t130 = rhob * t129;
		    t131 = pow_dd(&t126, &c_b4);
		    t134 = sqrt(t126);
/* Computing 2nd power */
		    d__1 = t127;
		    t136 = d__1 * d__1;
		    t138 = t131 * 11.12037486309468 + t127 * 
			    3.844746237447211 + t134 * 1.644733775567609 + 
			    t136 * .2405871291288192;
		    t141 = 32.1646831778707 / t138 + 1.;
		    t142 = log(t141);
		    t144 = t90 * .2 + 1.;
		    t145 = 1 / t144;
/* Computing 2nd power */
		    d__1 = t144;
		    t148 = d__1 * d__1;
		    t149 = 1 / t148;
		    t153 = 1 / t148 / t144;
/* Computing 2nd power */
		    d__1 = t148;
		    t156 = d__1 * d__1;
		    t157 = 1 / t156;
		    t160 = .489508 - t90 * .0521398 * t145 + t101 * .01731668 
			    * t149 - t109 * .01593976 * t153 + t118 * 
			    .003976496 * t157;
		    t161 = t142 * t160;
		    t164 = rhoa + rhob;
		    t165 = 1 / t164;
		    t166 = pow_dd(&t165, &c_b2);
		    t168 = t166 * .1325688999052018 + 1.;
		    t169 = pow_dd(&t165, &c_b4);
		    t172 = sqrt(t165);
/* Computing 2nd power */
		    d__1 = t166;
		    t174 = d__1 * d__1;
		    t176 = t169 * 5.98255043577108 + t166 * 2.225569421150687 
			    + t172 * .8004286349993634 + t174 * 
			    .1897004325747559;
		    t179 = 16.0818243221511 / t176 + 1.;
		    t180 = log(t179);
		    t182 = t168 * .062182 * t180;
		    t184 = t166 * .06901399211255825 + 1.;
		    t189 = t169 * 8.157414703487641 + t166 * 
			    2.247591863577616 + t172 * .4300972471276643 + 
			    t174 * .1911512595127338;
		    t192 = 29.60857464321668 / t189 + 1.;
		    t193 = log(t192);
		    t194 = t184 * t193;
		    t196 = rhoa - rhob * 1.;
		    t197 = t196 * t165;
		    t198 = t197 + 1.;
		    t199 = pow_dd(&t198, &c_b2);
		    t202 = 1. - t197 * 1.;
		    t203 = pow_dd(&t202, &c_b2);
		    t205 = t199 * t198 + t203 * t202 - 2.;
/* Computing 2nd power */
		    d__1 = t196;
		    t206 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t206;
		    t207 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t164;
		    t208 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t208;
		    t209 = d__1 * d__1;
		    t210 = 1 / t209;
		    t211 = t207 * t210;
		    t213 = 1. - t211 * 1.;
		    t216 = t194 * .03799574853701528 * t205 * t213;
		    t218 = t166 * .1274696188700087 + 1.;
		    t223 = t169 * 11.12037486309468 + t166 * 
			    3.844746237447211 + t172 * 1.644733775567609 + 
			    t174 * .2405871291288192;
		    t226 = 32.1646831778707 / t223 + 1.;
		    t227 = log(t226);
		    t230 = t218 * -.03109 * t227 + t182;
		    t231 = t230 * t205;
		    t233 = t231 * 1.923661050931536 * t211;
		    t240 = t164 * (-t182 + t216 + t233) + t50 * .03109 * t62 
			    + t130 * .03109 * t142;
		    t243 = t10 * .5 + t90 * .5;
		    t246 = t10 * .003 + 1. + t90 * .003;
		    t247 = 1 / t246;
/* Computing 2nd power */
		    d__1 = t243;
		    t250 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t246;
		    t251 = d__1 * d__1;
		    t252 = 1 / t251;
		    t255 = t250 * t243;
		    t257 = 1 / t251 / t246;
/* Computing 2nd power */
		    d__1 = t250;
		    t260 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t251;
		    t261 = d__1 * d__1;
		    t262 = 1 / t261;
		    t265 = t243 * .04157892 * t247 + .51473 - t250 * 
			    8.894628e-4 * t252 + t255 * 4.9917168e-6 * t257 - 
			    t260 * 1.46751264e-8 * t262;
		    zk[i__] = t5 * -.9305257363491 * t43 - t50 * .03109 * t81 
			    - t85 * .9305257363491 * t123 - t130 * .03109 * 
			    t161 + t240 * t265;
		    t269 = t6 * rhoa;
		    t272 = sigmaaa / t7 / t269;
		    t275 = t17 * t6;
		    t278 = t16 / t4 / t275;
		    t283 = t26 / t27 / rhoa;
		    t289 = t34 / t7 / t27 / t269;
		    t296 = t34 * sigmaaa / t4 / t27 / t275;
		    t298 = 1 / t39 / t12;
		    t304 = t49 * t62;
		    t307 = 1 / t56;
		    t308 = t46 * t307;
/* Computing 2nd power */
		    d__1 = t58;
		    t311 = d__1 * d__1;
		    t312 = 1 / t311;
/* Computing 2nd power */
		    d__1 = t51;
		    t314 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t314;
		    t315 = d__1 * d__1;
		    t318 = 1 / t6;
		    t329 = -1.853395810515781 / t315 / t51 * t318 - t307 * 
			    1.28158207914907 * t318 - .8223668877838045 / t54 
			    * t318 - .1603914194192128 / t47 * t318;
		    t330 = 1 / t61;
		    t344 = 1 / t76 / t64;
		    t352 = 1 / t208;
		    t353 = 1 / t174 * t352;
		    t354 = t353 * t180;
		    t355 = t354 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t176;
		    t356 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t169;
		    t359 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t359;
		    t360 = d__1 * d__1;
		    t363 = 1 / t360 / t169 * t352;
		    t367 = 1 / t172 * t352;
		    t370 = 1 / t166 * t352;
		    t375 = t168 / t356 * (t363 * -.99709173929518 - t353 * 
			    .7418564737168958 - t367 * .4002143174996817 - 
			    t370 * .1264669550498372) / t179;
		    t376 = t375 * 1.;
		    t380 = t353 * 8.740794299481065e-4 * t193 * t205 * t213;
/* Computing 2nd power */
		    d__1 = t189;
		    t381 = d__1 * d__1;
		    t394 = t184 * 1.124999956683108 / t381 * (t363 * 
			    -1.35956911724794 - t353 * .7491972878592054 - 
			    t367 * .2150486235638321 - t370 * 
			    .1274341730084892) / t192 * t205 * t213;
		    t395 = t196 * t352;
		    t396 = t395 * 1.;
		    t400 = t165 * 1.;
		    t404 = t199 * 1.333333333333333 * (t165 - t396) + t203 * 
			    1.333333333333333 * (-t400 + t395);
		    t409 = t206 * t196 * t210;
		    t410 = t409 * 4.;
		    t413 = t207 / t209 / t164;
		    t414 = t413 * 4.;
/* Computing 2nd power */
		    d__1 = t223;
		    t421 = d__1 * d__1;
		    t438 = (t353 * .001321010150222857 * t227 + t218 * 1. / 
			    t421 * (t363 * -1.853395810515781 - t353 * 
			    1.28158207914907 - t367 * .8223668877838045 - 
			    t370 * .1603914194192128) / t226 - t354 * 
			    .002747799777968419 - t375 * 1.) * 
			    1.923661050931536 * t205 * t211;
		    t442 = t231 * t409;
		    t445 = t231 * 7.694644203726145 * t413;
		    t459 = t243 * t252;
		    t462 = t250 * t257;
		    t465 = t255 * t262;
		    t470 = t260 / t261 / t246;
		    s1 = t4 * -1.2407009817988 * t43 - t5 * .9305257363491 * (
			    t272 * .007970293333333333 * t13 - t278 * 
			    4.65232e-4 * t23 + t283 * 3.836422826666667e-6 * 
			    t31 - t289 * 1.161168213333333e-8 * t40 + t296 * 
			    1.279841621333333e-11 * t298) - t304 * .03109 * 
			    t80 + t308 * .001321010150222857 * t81;
		    vrhoa[i__] = s1 + t50 * 1. * t312 * t329 * t330 * t80 - 
			    t50 * .03109 * t62 * (t272 * .1390394666666667 * 
			    t65 - t278 * .12016352 * t69 + t283 * 
			    .1459892053333333 * t73 - t289 * 
			    .06791957333333333 * t77 + t296 * 
			    .008483191466666667 * t344) + (-t182 + t216 + 
			    t233 + t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t404 * t213 + t194 * 
			    .03799574853701528 * t205 * (-t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t404 * t211 + t442 * 
			    7.694644203726145 - t445) + t304 * .03109 - t308 *
			     .001321010150222857 * t62 - t50 * 1. * t312 * 
			    t329 * t330) * t265 + t240 * (t272 * -.05543856 * 
			    t247 + t459 * .00270453216 * t272 - t462 * 
			    3.4198272e-5 * t272 + t465 * 1.98068544e-7 * t272 
			    - t470 * 4.696040448e-10 * t272);
		    t477 = t86 * rhob;
		    t480 = sigmabb / t87 / t477;
		    t483 = t97 * t86;
		    t486 = t96 / t84 / t483;
		    t491 = t106 / t107 / rhob;
		    t497 = t114 / t87 / t107 / t477;
		    t504 = t114 * sigmabb / t84 / t107 / t483;
		    t506 = 1 / t119 / t92;
		    t512 = t129 * t142;
		    t515 = 1 / t136;
		    t516 = t126 * t515;
/* Computing 2nd power */
		    d__1 = t138;
		    t519 = d__1 * d__1;
		    t520 = 1 / t519;
/* Computing 2nd power */
		    d__1 = t131;
		    t522 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t522;
		    t523 = d__1 * d__1;
		    t526 = 1 / t86;
		    t537 = -1.853395810515781 / t523 / t131 * t526 - t515 * 
			    1.28158207914907 * t526 - .8223668877838045 / 
			    t134 * t526 - .1603914194192128 / t127 * t526;
		    t538 = 1 / t141;
		    t552 = 1 / t156 / t144;
		    t565 = t199 * 1.333333333333333 * (-t400 - t396) + t203 * 
			    1.333333333333333 * (t165 + t395);
		    s1 = t84 * -1.2407009817988 * t123 - t85 * .9305257363491 
			    * (t480 * .007970293333333333 * t93 - t486 * 
			    4.65232e-4 * t103 + t491 * 3.836422826666667e-6 * 
			    t111 - t497 * 1.161168213333333e-8 * t120 + t504 *
			     1.279841621333333e-11 * t506) - t512 * .03109 * 
			    t160 + t516 * .001321010150222857 * t161;
		    vrhob[i__] = s1 + t130 * 1. * t520 * t537 * t538 * t160 - 
			    t130 * .03109 * t142 * (t480 * .1390394666666667 *
			     t145 - t486 * .12016352 * t149 + t491 * 
			    .1459892053333333 * t153 - t497 * 
			    .06791957333333333 * t157 + t504 * 
			    .008483191466666667 * t552) + (-t182 + t216 + 
			    t233 + t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t565 * t213 + t194 * 
			    .03799574853701528 * t205 * (t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t565 * t211 - t442 * 
			    7.694644203726145 - t445) + t512 * .03109 - t516 *
			     .001321010150222857 * t142 - t130 * 1. * t520 * 
			    t537 * t538) * t265 + t240 * (t480 * -.05543856 * 
			    t247 + t459 * .00270453216 * t480 - t462 * 
			    3.4198272e-5 * t480 + t465 * 1.98068544e-7 * t480 
			    - t470 * 4.696040448e-10 * t480);
		    t602 = sigmaaa * t20;
		    t605 = t16 * t28;
		    t608 = t26 * t37;
		    t614 = t34 / t4 / t27 / t18;
		    vsigmaaa[i__] = t5 * -.9305257363491 * (t9 * -.00298886 * 
			    t13 + t602 * 1.74462e-4 * t23 - t605 * 
			    1.43865856e-6 * t31 + t608 * 4.3543808e-9 * t40 - 
			    t614 * 4.79940608e-12 * t298) - t50 * .03109 * 
			    t62 * (t9 * -.0521398 * t65 + t602 * .04506132 * 
			    t69 - t605 * .054745952 * t73 + t608 * .02546984 *
			     t77 - t614 * .0031811968 * t344) + t240 * (t9 * 
			    .02078946 * t247 - t459 * .00101419956 * t9 + 
			    t462 * 1.2824352e-5 * t9 - t465 * 7.4275704e-8 * 
			    t9 + t470 * 1.761015168e-10 * t9);
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t660 = t114 / t84 / t107 / t98;
		    vsigmabb[i__] = t85 * -.9305257363491 * (t89 * -.00298886 
			    * t93 + t648 * 1.74462e-4 * t103 - t651 * 
			    1.43865856e-6 * t111 + t654 * 4.3543808e-9 * t120 
			    - t660 * 4.79940608e-12 * t506) - t130 * .03109 * 
			    t142 * (t89 * -.0521398 * t145 + t648 * .04506132 
			    * t149 - t651 * .054745952 * t153 + t654 * 
			    .02546984 * t157 - t660 * .0031811968 * t552) + 
			    t240 * (t89 * .02078946 * t247 - t459 * 
			    .00101419956 * t89 + t462 * 1.2824352e-5 * t89 - 
			    t465 * 7.4275704e-8 * t89 + t470 * 
			    1.761015168e-10 * t89);
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
		    t16 = t15 * rhob;
		    t18 = 1 / t2 / t16;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t14 * sigmabb;
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t26 = 1 / t25;
		    t27 = t24 * t26;
		    t29 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t33 = t25 * t4;
		    t35 = 1 / t5 / t33;
		    t36 = t32 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t38 = 1 / t37;
		    t41 = 1.09163 - t8 * .00298886 * t11 + t19 * 8.125328e-5 *
			     t21 - t27 * 2.6287744e-7 * t29 + t36 * 
			    2.9996288e-10 * t38;
		    t44 = 1 / rhob;
		    t45 = pow_dd(&t44, &c_b2);
		    t47 = t45 * .1274696188700087 + 1.;
		    t48 = rhob * t47;
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t56 = t49 * 11.12037486309468 + t45 * 3.844746237447211 + 
			    t52 * 1.644733775567609 + t54 * .2405871291288192;
		    t59 = 32.1646831778707 / t56 + 1.;
		    t60 = log(t59);
		    t62 = t8 * .2 + 1.;
		    t63 = 1 / t62;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t71 = 1 / t66 / t62;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t78 = .489508 - t8 * .0521398 * t63 + t19 * .01731668 * 
			    t67 - t27 * .01593976 * t71 + t36 * .003976496 * 
			    t75;
		    t79 = t60 * t78;
		    zk[i__] = t3 * -.9305257363491 * t41 - t48 * .03109 * t79;
		    vrhoa[i__] = 0.;
		    t84 = t4 * rhob;
		    t87 = sigmabb / t5 / t84;
		    t90 = t15 * t4;
		    t93 = t14 / t2 / t90;
		    t98 = t24 / t25 / rhob;
		    t104 = t32 / t5 / t25 / t84;
		    t107 = t32 * sigmabb;
		    t111 = t107 / t2 / t25 / t90;
		    t113 = 1 / t37 / t10;
		    t116 = t87 * .007970293333333333 * t11 - t93 * 4.65232e-4 
			    * t21 + t98 * 3.836422826666667e-6 * t29 - t104 * 
			    1.161168213333333e-8 * t38 + t111 * 
			    1.279841621333333e-11 * t113;
		    t119 = t47 * t60;
		    t122 = 1 / t54;
		    t123 = t44 * t122;
/* Computing 2nd power */
		    d__1 = t56;
		    t126 = d__1 * d__1;
		    t127 = 1 / t126;
		    t128 = t48 * t127;
/* Computing 2nd power */
		    d__1 = t49;
		    t129 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t129;
		    t130 = d__1 * d__1;
		    t131 = t130 * t49;
		    t132 = 1 / t131;
		    t133 = 1 / t4;
		    t138 = 1 / t52;
		    t141 = 1 / t45;
		    t144 = t132 * -1.853395810515781 * t133 - t122 * 
			    1.28158207914907 * t133 - t138 * 
			    .8223668877838045 * t133 - t141 * 
			    .1603914194192128 * t133;
		    t145 = 1 / t59;
		    t146 = t144 * t145;
		    t147 = t146 * t78;
		    t159 = 1 / t74 / t62;
		    t162 = t87 * .1390394666666667 * t63 - t93 * .12016352 * 
			    t67 + t98 * .1459892053333333 * t71 - t104 * 
			    .06791957333333333 * t75 + t111 * 
			    .008483191466666667 * t159;
		    t163 = t60 * t162;
		    vrhob[i__] = t2 * -1.2407009817988 * t41 - t3 * 
			    .9305257363491 * t116 - t119 * .03109 * t78 + 
			    t123 * .001321010150222857 * t79 + t128 * 1. * 
			    t147 - t48 * .03109 * t163;
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t168 = sigmabb * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t179 = 1 / t2 / t25 / t16;
		    t180 = t32 * t179;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.00298886 * 
			    t11 + t168 * 1.74462e-4 * t21 - t171 * 
			    1.43865856e-6 * t29 + t174 * 4.3543808e-9 * t38 - 
			    t180 * 4.79940608e-12 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.0521398 * t63 + t168 * .04506132 * 
			    t67 - t171 * .054745952 * t71 + t174 * .02546984 *
			     t75 - t180 * .0031811968 * t159);
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t207 = sigmabb / t5 / t15;
		    t210 = t15 * t84;
		    t213 = t14 / t2 / t210;
		    t217 = t24 / t33;
		    t223 = t32 / t5 / t25 / t15;
		    t229 = t107 / t2 / t25 / t210;
/* Computing 2nd power */
		    d__1 = t25;
		    t233 = d__1 * d__1;
		    t236 = t32 * t14 / t233 / t4;
		    t238 = 1 / t37 / t20;
		    t249 = 1 / t84;
		    t251 = 1 / t54 / t44;
/* Computing 2nd power */
		    d__1 = t144;
		    t263 = d__1 * d__1;
		    t270 = 1 / t15;
/* Computing 2nd power */
		    d__1 = t126;
		    t296 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t59;
		    t299 = d__1 * d__1;
		    t319 = 1 / t74 / t66;
		    s1 = -.4135669939329333 / t5 * t41 - t2 * 2.4814019635976 
			    * t116 - t3 * .9305257363491 * (t207 * 
			    -.02922440888888889 * t11 + t213 * 
			    .003031485795555556 * t21 - t217 * 
			    4.445275477333333e-5 * t29 + t223 * 
			    2.582351553422222e-7 * t38 - t229 * 
			    6.788757367466667e-10 * t113 + t236 * 
			    6.825821980444444e-13 * t238) + t47 * 2. * t127 * 
			    t147 - t119 * .06218 * t162 + t249 * 
			    8.806734334819047e-4 * t251 * t79;
		    v2rhob2[i__] = s1 - t123 * .08497974591333914 * t127 * 
			    t147 + t123 * .002642020300445714 * t163 - t48 * 
			    2. / t126 / t56 * t263 * t145 * t78 + t128 * 1. * 
			    (-1.544496508763151 / t131 / t44 * t270 + t132 * 
			    3.706791621031562 * t249 - t251 * 
			    .854388052766047 * t270 + t122 * 
			    2.563164158298141 * t249 - .4111834438919023 / 
			    t52 / t44 * t270 + t138 * 1.644733775567609 * 
			    t249 - .05346380647307093 / t45 / t44 * t270 + 
			    t141 * .3207828388384256 * t249) * t145 * t78 + 
			    t48 * 32.1646831778707 / t296 * t263 / t299 * t78 
			    + t128 * 2. * t146 * t162 - t48 * .03109 * t60 * (
			    t207 * -.5098113777777778 * t63 + t213 * 
			    .8351900088888889 * t67 - t217 * 
			    1.442077269333333 * t71 + t223 * 
			    1.025977750755556 * t75 - t229 * .2664875008 * 
			    t159 + t236 * .02262184391111111 * t319);
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t328 = sigmabb * t26;
		    t331 = t14 * t35;
		    t334 = t24 * t179;
		    t338 = t32 / t233;
		    v2sigmabb2[i__] = t3 * -.9305257363491 * (t18 * 
			    1.8641744e-4 * t21 - t328 * 4.27301312e-6 * t29 + 
			    t331 * 3.032704512e-8 * t38 - t334 * 
			    8.886771712e-11 * t113 + t338 * 9.59881216e-14 * 
			    t238) - t48 * .03109 * t60 * (t18 * .05548928 * 
			    t67 - t328 * .127516432 * t71 + t331 * 
			    .1092570912 * t75 - t334 * .0331006592 * t159 + 
			    t338 * .0031811968 * t319);
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
		    t16 = t15 * rhoa;
		    t18 = 1 / t2 / t16;
		    t19 = t14 * t18;
/* Computing 2nd power */
		    d__1 = t10;
		    t20 = d__1 * d__1;
		    t21 = 1 / t20;
		    t24 = t14 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t15;
		    t25 = d__1 * d__1;
		    t26 = 1 / t25;
		    t27 = t24 * t26;
		    t29 = 1 / t20 / t10;
/* Computing 2nd power */
		    d__1 = t14;
		    t32 = d__1 * d__1;
		    t33 = t25 * t4;
		    t35 = 1 / t5 / t33;
		    t36 = t32 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
		    t38 = 1 / t37;
		    t41 = 1.09163 - t8 * .00298886 * t11 + t19 * 8.125328e-5 *
			     t21 - t27 * 2.6287744e-7 * t29 + t36 * 
			    2.9996288e-10 * t38;
		    t44 = 1 / rhoa;
		    t45 = pow_dd(&t44, &c_b2);
		    t47 = t45 * .1274696188700087 + 1.;
		    t48 = rhoa * t47;
		    t49 = pow_dd(&t44, &c_b4);
		    t52 = sqrt(t44);
/* Computing 2nd power */
		    d__1 = t45;
		    t54 = d__1 * d__1;
		    t56 = t49 * 11.12037486309468 + t45 * 3.844746237447211 + 
			    t52 * 1.644733775567609 + t54 * .2405871291288192;
		    t59 = 32.1646831778707 / t56 + 1.;
		    t60 = log(t59);
		    t62 = t8 * .2 + 1.;
		    t63 = 1 / t62;
/* Computing 2nd power */
		    d__1 = t62;
		    t66 = d__1 * d__1;
		    t67 = 1 / t66;
		    t71 = 1 / t66 / t62;
/* Computing 2nd power */
		    d__1 = t66;
		    t74 = d__1 * d__1;
		    t75 = 1 / t74;
		    t78 = .489508 - t8 * .0521398 * t63 + t19 * .01731668 * 
			    t67 - t27 * .01593976 * t71 + t36 * .003976496 * 
			    t75;
		    t79 = t60 * t78;
		    zk[i__] = t3 * -.9305257363491 * t41 - t48 * .03109 * t79;
		    t84 = t4 * rhoa;
		    t87 = sigmaaa / t5 / t84;
		    t90 = t15 * t4;
		    t93 = t14 / t2 / t90;
		    t98 = t24 / t25 / rhoa;
		    t104 = t32 / t5 / t25 / t84;
		    t107 = t32 * sigmaaa;
		    t111 = t107 / t2 / t25 / t90;
		    t113 = 1 / t37 / t10;
		    t116 = t87 * .007970293333333333 * t11 - t93 * 4.65232e-4 
			    * t21 + t98 * 3.836422826666667e-6 * t29 - t104 * 
			    1.161168213333333e-8 * t38 + t111 * 
			    1.279841621333333e-11 * t113;
		    t119 = t47 * t60;
		    t122 = 1 / t54;
		    t123 = t44 * t122;
/* Computing 2nd power */
		    d__1 = t56;
		    t126 = d__1 * d__1;
		    t127 = 1 / t126;
		    t128 = t48 * t127;
/* Computing 2nd power */
		    d__1 = t49;
		    t129 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t129;
		    t130 = d__1 * d__1;
		    t131 = t130 * t49;
		    t132 = 1 / t131;
		    t133 = 1 / t4;
		    t138 = 1 / t52;
		    t141 = 1 / t45;
		    t144 = t132 * -1.853395810515781 * t133 - t122 * 
			    1.28158207914907 * t133 - t138 * 
			    .8223668877838045 * t133 - t141 * 
			    .1603914194192128 * t133;
		    t145 = 1 / t59;
		    t146 = t144 * t145;
		    t147 = t146 * t78;
		    t159 = 1 / t74 / t62;
		    t162 = t87 * .1390394666666667 * t63 - t93 * .12016352 * 
			    t67 + t98 * .1459892053333333 * t71 - t104 * 
			    .06791957333333333 * t75 + t111 * 
			    .008483191466666667 * t159;
		    t163 = t60 * t162;
		    vrhoa[i__] = t2 * -1.2407009817988 * t41 - t3 * 
			    .9305257363491 * t116 - t119 * .03109 * t78 + 
			    t123 * .001321010150222857 * t79 + t128 * 1. * 
			    t147 - t48 * .03109 * t163;
		    vrhob[i__] = 0.;
		    t168 = sigmaaa * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t179 = 1 / t2 / t25 / t16;
		    t180 = t32 * t179;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.00298886 * 
			    t11 + t168 * 1.74462e-4 * t21 - t171 * 
			    1.43865856e-6 * t29 + t174 * 4.3543808e-9 * t38 - 
			    t180 * 4.79940608e-12 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.0521398 * t63 + t168 * .04506132 * 
			    t67 - t171 * .054745952 * t71 + t174 * .02546984 *
			     t75 - t180 * .0031811968 * t159);
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t207 = sigmaaa / t5 / t15;
		    t210 = t15 * t84;
		    t213 = t14 / t2 / t210;
		    t217 = t24 / t33;
		    t223 = t32 / t5 / t25 / t15;
		    t229 = t107 / t2 / t25 / t210;
/* Computing 2nd power */
		    d__1 = t25;
		    t233 = d__1 * d__1;
		    t236 = t32 * t14 / t233 / t4;
		    t238 = 1 / t37 / t20;
		    t249 = 1 / t84;
		    t251 = 1 / t54 / t44;
/* Computing 2nd power */
		    d__1 = t144;
		    t263 = d__1 * d__1;
		    t270 = 1 / t15;
/* Computing 2nd power */
		    d__1 = t126;
		    t296 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t59;
		    t299 = d__1 * d__1;
		    t319 = 1 / t74 / t66;
		    s1 = -.4135669939329333 / t5 * t41 - t2 * 2.4814019635976 
			    * t116 - t3 * .9305257363491 * (t207 * 
			    -.02922440888888889 * t11 + t213 * 
			    .003031485795555556 * t21 - t217 * 
			    4.445275477333333e-5 * t29 + t223 * 
			    2.582351553422222e-7 * t38 - t229 * 
			    6.788757367466667e-10 * t113 + t236 * 
			    6.825821980444444e-13 * t238) + t47 * 2. * t127 * 
			    t147 - t119 * .06218 * t162 + t249 * 
			    8.806734334819047e-4 * t251 * t79;
		    v2rhoa2[i__] = s1 - t123 * .08497974591333914 * t127 * 
			    t147 + t123 * .002642020300445714 * t163 - t48 * 
			    2. / t126 / t56 * t263 * t145 * t78 + t128 * 1. * 
			    (-1.544496508763151 / t131 / t44 * t270 + t132 * 
			    3.706791621031562 * t249 - t251 * 
			    .854388052766047 * t270 + t122 * 
			    2.563164158298141 * t249 - .4111834438919023 / 
			    t52 / t44 * t270 + t138 * 1.644733775567609 * 
			    t249 - .05346380647307093 / t45 / t44 * t270 + 
			    t141 * .3207828388384256 * t249) * t145 * t78 + 
			    t48 * 32.1646831778707 / t296 * t263 / t299 * t78 
			    + t128 * 2. * t146 * t162 - t48 * .03109 * t60 * (
			    t207 * -.5098113777777778 * t63 + t213 * 
			    .8351900088888889 * t67 - t217 * 
			    1.442077269333333 * t71 + t223 * 
			    1.025977750755556 * t75 - t229 * .2664875008 * 
			    t159 + t236 * .02262184391111111 * t319);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t328 = sigmaaa * t26;
		    t331 = t14 * t35;
		    t334 = t24 * t179;
		    t338 = t32 / t233;
		    v2sigmaaa2[i__] = t3 * -.9305257363491 * (t18 * 
			    1.8641744e-4 * t21 - t328 * 4.27301312e-6 * t29 + 
			    t331 * 3.032704512e-8 * t38 - t334 * 
			    8.886771712e-11 * t113 + t338 * 9.59881216e-14 * 
			    t238) - t48 * .03109 * t60 * (t18 * .05548928 * 
			    t67 - t328 * .127516432 * t71 + t331 * 
			    .1092570912 * t75 - t334 * .0331006592 * t159 + 
			    t338 * .0031811968 * t319);
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
		    t18 = t17 * rhoa;
		    t20 = 1 / t4 / t18;
		    t21 = t16 * t20;
/* Computing 2nd power */
		    d__1 = t12;
		    t22 = d__1 * d__1;
		    t23 = 1 / t22;
		    t26 = t16 * sigmaaa;
/* Computing 2nd power */
		    d__1 = t17;
		    t27 = d__1 * d__1;
		    t28 = 1 / t27;
		    t29 = t26 * t28;
		    t31 = 1 / t22 / t12;
/* Computing 2nd power */
		    d__1 = t16;
		    t34 = d__1 * d__1;
		    t35 = t27 * t6;
		    t37 = 1 / t7 / t35;
		    t38 = t34 * t37;
/* Computing 2nd power */
		    d__1 = t22;
		    t39 = d__1 * d__1;
		    t40 = 1 / t39;
		    t43 = 1.09163 - t10 * .00298886 * t13 + t21 * 8.125328e-5 
			    * t23 - t29 * 2.6287744e-7 * t31 + t38 * 
			    2.9996288e-10 * t40;
		    t46 = 1 / rhoa;
		    t47 = pow_dd(&t46, &c_b2);
		    t49 = t47 * .1274696188700087 + 1.;
		    t50 = rhoa * t49;
		    t51 = pow_dd(&t46, &c_b4);
		    t54 = sqrt(t46);
/* Computing 2nd power */
		    d__1 = t47;
		    t56 = d__1 * d__1;
		    t58 = t51 * 11.12037486309468 + t47 * 3.844746237447211 + 
			    t54 * 1.644733775567609 + t56 * .2405871291288192;
		    t61 = 32.1646831778707 / t58 + 1.;
		    t62 = log(t61);
		    t64 = t10 * .2 + 1.;
		    t65 = 1 / t64;
/* Computing 2nd power */
		    d__1 = t64;
		    t68 = d__1 * d__1;
		    t69 = 1 / t68;
		    t73 = 1 / t68 / t64;
/* Computing 2nd power */
		    d__1 = t68;
		    t76 = d__1 * d__1;
		    t77 = 1 / t76;
		    t80 = .489508 - t10 * .0521398 * t65 + t21 * .01731668 * 
			    t69 - t29 * .01593976 * t73 + t38 * .003976496 * 
			    t77;
		    t81 = t62 * t80;
		    t84 = pow_dd(&rhob, &c_b2);
		    t85 = t84 * rhob;
/* Computing 2nd power */
		    d__1 = rhob;
		    t86 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t84;
		    t87 = d__1 * d__1;
		    t89 = 1 / t87 / t86;
		    t90 = sigmabb * t89;
		    t92 = t90 * .004 + 1.;
		    t93 = 1 / t92;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t96 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t86;
		    t97 = d__1 * d__1;
		    t98 = t97 * rhob;
		    t100 = 1 / t84 / t98;
		    t101 = t96 * t100;
/* Computing 2nd power */
		    d__1 = t92;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
		    t106 = t96 * sigmabb;
/* Computing 2nd power */
		    d__1 = t97;
		    t107 = d__1 * d__1;
		    t108 = 1 / t107;
		    t109 = t106 * t108;
		    t111 = 1 / t102 / t92;
/* Computing 2nd power */
		    d__1 = t96;
		    t114 = d__1 * d__1;
		    t115 = t107 * t86;
		    t117 = 1 / t87 / t115;
		    t118 = t114 * t117;
/* Computing 2nd power */
		    d__1 = t102;
		    t119 = d__1 * d__1;
		    t120 = 1 / t119;
		    t123 = 1.09163 - t90 * .00298886 * t93 + t101 * 
			    8.125328e-5 * t103 - t109 * 2.6287744e-7 * t111 + 
			    t118 * 2.9996288e-10 * t120;
		    t126 = 1 / rhob;
		    t127 = pow_dd(&t126, &c_b2);
		    t129 = t127 * .1274696188700087 + 1.;
		    t130 = rhob * t129;
		    t131 = pow_dd(&t126, &c_b4);
		    t134 = sqrt(t126);
/* Computing 2nd power */
		    d__1 = t127;
		    t136 = d__1 * d__1;
		    t138 = t131 * 11.12037486309468 + t127 * 
			    3.844746237447211 + t134 * 1.644733775567609 + 
			    t136 * .2405871291288192;
		    t141 = 32.1646831778707 / t138 + 1.;
		    t142 = log(t141);
		    t144 = t90 * .2 + 1.;
		    t145 = 1 / t144;
/* Computing 2nd power */
		    d__1 = t144;
		    t148 = d__1 * d__1;
		    t149 = 1 / t148;
		    t153 = 1 / t148 / t144;
/* Computing 2nd power */
		    d__1 = t148;
		    t156 = d__1 * d__1;
		    t157 = 1 / t156;
		    t160 = .489508 - t90 * .0521398 * t145 + t101 * .01731668 
			    * t149 - t109 * .01593976 * t153 + t118 * 
			    .003976496 * t157;
		    t161 = t142 * t160;
		    t164 = rhoa + rhob;
		    t165 = 1 / t164;
		    t166 = pow_dd(&t165, &c_b2);
		    t168 = t166 * .1325688999052018 + 1.;
		    t169 = pow_dd(&t165, &c_b4);
		    t172 = sqrt(t165);
/* Computing 2nd power */
		    d__1 = t166;
		    t174 = d__1 * d__1;
		    t176 = t169 * 5.98255043577108 + t166 * 2.225569421150687 
			    + t172 * .8004286349993634 + t174 * 
			    .1897004325747559;
		    t179 = 16.0818243221511 / t176 + 1.;
		    t180 = log(t179);
		    t182 = t168 * .062182 * t180;
		    t184 = t166 * .06901399211255825 + 1.;
		    t189 = t169 * 8.157414703487641 + t166 * 
			    2.247591863577616 + t172 * .4300972471276643 + 
			    t174 * .1911512595127338;
		    t192 = 29.60857464321668 / t189 + 1.;
		    t193 = log(t192);
		    t194 = t184 * t193;
		    t196 = rhoa - rhob * 1.;
		    t197 = t196 * t165;
		    t198 = t197 + 1.;
		    t199 = pow_dd(&t198, &c_b2);
		    t202 = 1. - t197 * 1.;
		    t203 = pow_dd(&t202, &c_b2);
		    t205 = t199 * t198 + t203 * t202 - 2.;
/* Computing 2nd power */
		    d__1 = t196;
		    t206 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t206;
		    t207 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t164;
		    t208 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t208;
		    t209 = d__1 * d__1;
		    t210 = 1 / t209;
		    t211 = t207 * t210;
		    t213 = 1. - t211 * 1.;
		    t214 = t205 * t213;
		    t216 = t194 * .03799574853701528 * t214;
		    t218 = t166 * .1274696188700087 + 1.;
		    t223 = t169 * 11.12037486309468 + t166 * 
			    3.844746237447211 + t172 * 1.644733775567609 + 
			    t174 * .2405871291288192;
		    t226 = 32.1646831778707 / t223 + 1.;
		    t227 = log(t226);
		    t230 = t218 * -.03109 * t227 + t182;
		    t231 = t230 * t205;
		    t233 = t231 * 1.923661050931536 * t211;
		    t240 = t164 * (-t182 + t216 + t233) + t50 * .03109 * t62 
			    + t130 * .03109 * t142;
		    t243 = t10 * .5 + t90 * .5;
		    t246 = t10 * .003 + 1. + t90 * .003;
		    t247 = 1 / t246;
/* Computing 2nd power */
		    d__1 = t243;
		    t250 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t246;
		    t251 = d__1 * d__1;
		    t252 = 1 / t251;
		    t255 = t250 * t243;
		    t257 = 1 / t251 / t246;
/* Computing 2nd power */
		    d__1 = t250;
		    t260 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t251;
		    t261 = d__1 * d__1;
		    t262 = 1 / t261;
		    t265 = t243 * .04157892 * t247 + .51473 - t250 * 
			    8.894628e-4 * t252 + t255 * 4.9917168e-6 * t257 - 
			    t260 * 1.46751264e-8 * t262;
		    zk[i__] = t5 * -.9305257363491 * t43 - t50 * .03109 * t81 
			    - t85 * .9305257363491 * t123 - t130 * .03109 * 
			    t161 + t240 * t265;
		    t269 = t6 * rhoa;
		    t271 = 1 / t7 / t269;
		    t272 = sigmaaa * t271;
		    t275 = t17 * t6;
		    t277 = 1 / t4 / t275;
		    t278 = t16 * t277;
		    t282 = 1 / t27 / rhoa;
		    t283 = t26 * t282;
		    t288 = 1 / t7 / t27 / t269;
		    t289 = t34 * t288;
		    t292 = t34 * sigmaaa;
		    t295 = 1 / t4 / t27 / t275;
		    t296 = t292 * t295;
		    t298 = 1 / t39 / t12;
		    t301 = t272 * .007970293333333333 * t13 - t278 * 
			    4.65232e-4 * t23 + t283 * 3.836422826666667e-6 * 
			    t31 - t289 * 1.161168213333333e-8 * t40 + t296 * 
			    1.279841621333333e-11 * t298;
		    t304 = t49 * t62;
		    t307 = 1 / t56;
		    t308 = t46 * t307;
/* Computing 2nd power */
		    d__1 = t58;
		    t311 = d__1 * d__1;
		    t312 = 1 / t311;
		    t313 = t50 * t312;
/* Computing 2nd power */
		    d__1 = t51;
		    t314 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t314;
		    t315 = d__1 * d__1;
		    t316 = t315 * t51;
		    t317 = 1 / t316;
		    t318 = 1 / t6;
		    t323 = 1 / t54;
		    t326 = 1 / t47;
		    t329 = t317 * -1.853395810515781 * t318 - t307 * 
			    1.28158207914907 * t318 - t323 * 
			    .8223668877838045 * t318 - t326 * 
			    .1603914194192128 * t318;
		    t330 = 1 / t61;
		    t331 = t329 * t330;
		    t332 = t331 * t80;
		    t344 = 1 / t76 / t64;
		    t347 = t272 * .1390394666666667 * t65 - t278 * .12016352 *
			     t69 + t283 * .1459892053333333 * t73 - t289 * 
			    .06791957333333333 * t77 + t296 * 
			    .008483191466666667 * t344;
		    t348 = t62 * t347;
		    t351 = 1 / t174;
		    t352 = 1 / t208;
		    t353 = t351 * t352;
		    t354 = t353 * t180;
		    t355 = t354 * .002747799777968419;
/* Computing 2nd power */
		    d__1 = t176;
		    t356 = d__1 * d__1;
		    t357 = 1 / t356;
		    t358 = t168 * t357;
/* Computing 2nd power */
		    d__1 = t169;
		    t359 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t359;
		    t360 = d__1 * d__1;
		    t361 = t360 * t169;
		    t362 = 1 / t361;
		    t363 = t362 * t352;
		    t366 = 1 / t172;
		    t367 = t366 * t352;
		    t369 = 1 / t166;
		    t370 = t369 * t352;
		    t372 = t363 * -.99709173929518 - t353 * .7418564737168958 
			    - t367 * .4002143174996817 - t370 * 
			    .1264669550498372;
		    t373 = 1 / t179;
		    t375 = t358 * t372 * t373;
		    t376 = t375 * 1.;
		    t377 = t193 * t205;
		    t378 = t377 * t213;
		    t379 = t353 * t378;
		    t380 = t379 * 8.740794299481065e-4;
/* Computing 2nd power */
		    d__1 = t189;
		    t381 = d__1 * d__1;
		    t382 = 1 / t381;
		    t383 = t184 * t382;
		    t388 = t363 * -1.35956911724794 - t353 * 
			    .7491972878592054 - t367 * .2150486235638321 - 
			    t370 * .1274341730084892;
		    t389 = t383 * t388;
		    t390 = 1 / t192;
		    t391 = t390 * t205;
		    t392 = t391 * t213;
		    t393 = t389 * t392;
		    t394 = t393 * 1.124999956683108;
		    t395 = t196 * t352;
		    t396 = t395 * 1.;
		    t397 = t165 - t396;
		    t400 = t165 * 1.;
		    t401 = -t400 + t395;
		    t404 = t199 * 1.333333333333333 * t397 + t203 * 
			    1.333333333333333 * t401;
		    t406 = t194 * t404 * t213;
		    t407 = t406 * .03799574853701528;
		    t408 = t206 * t196;
		    t409 = t408 * t210;
		    t410 = t409 * 4.;
		    t412 = 1 / t209 / t164;
		    t413 = t207 * t412;
		    t414 = t413 * 4.;
		    t415 = -t410 + t414;
		    t417 = t194 * t205 * t415;
		    t418 = t417 * .03799574853701528;
/* Computing 2nd power */
		    d__1 = t223;
		    t421 = d__1 * d__1;
		    t422 = 1 / t421;
		    t423 = t218 * t422;
		    t428 = t363 * -1.853395810515781 - t353 * 
			    1.28158207914907 - t367 * .8223668877838045 - 
			    t370 * .1603914194192128;
		    t429 = 1 / t226;
		    t435 = t353 * .001321010150222857 * t227 + t423 * 1. * 
			    t428 * t429 - t354 * .002747799777968419 - t375 * 
			    1.;
		    t436 = t435 * t205;
		    t437 = t436 * t211;
		    t438 = t437 * 1.923661050931536;
		    t439 = t230 * t404;
		    t440 = t439 * t211;
		    t441 = t440 * 1.923661050931536;
		    t442 = t231 * t409;
		    t444 = t231 * t413;
		    t445 = t444 * 7.694644203726145;
		    t452 = t312 * t329 * t330;
		    t455 = -t182 + t216 + t233 + t164 * (t355 + t376 - t380 - 
			    t394 + t407 + t418 + t438 + t441 + t442 * 
			    7.694644203726145 - t445) + t304 * .03109 - t308 *
			     .001321010150222857 * t62 - t50 * 1. * t452;
		    t459 = t243 * t252;
		    t462 = t250 * t257;
		    t465 = t255 * t262;
		    t469 = 1 / t261 / t246;
		    t470 = t260 * t469;
		    t473 = t272 * -.05543856 * t247 + t459 * .00270453216 * 
			    t272 - t462 * 3.4198272e-5 * t272 + t465 * 
			    1.98068544e-7 * t272 - t470 * 4.696040448e-10 * 
			    t272;
		    vrhoa[i__] = t4 * -1.2407009817988 * t43 - t5 * 
			    .9305257363491 * t301 - t304 * .03109 * t80 + 
			    t308 * .001321010150222857 * t81 + t313 * 1. * 
			    t332 - t50 * .03109 * t348 + t455 * t265 + t240 * 
			    t473;
		    t477 = t86 * rhob;
		    t479 = 1 / t87 / t477;
		    t480 = sigmabb * t479;
		    t483 = t97 * t86;
		    t485 = 1 / t84 / t483;
		    t486 = t96 * t485;
		    t490 = 1 / t107 / rhob;
		    t491 = t106 * t490;
		    t496 = 1 / t87 / t107 / t477;
		    t497 = t114 * t496;
		    t500 = t114 * sigmabb;
		    t503 = 1 / t84 / t107 / t483;
		    t504 = t500 * t503;
		    t506 = 1 / t119 / t92;
		    t509 = t480 * .007970293333333333 * t93 - t486 * 
			    4.65232e-4 * t103 + t491 * 3.836422826666667e-6 * 
			    t111 - t497 * 1.161168213333333e-8 * t120 + t504 *
			     1.279841621333333e-11 * t506;
		    t512 = t129 * t142;
		    t515 = 1 / t136;
		    t516 = t126 * t515;
/* Computing 2nd power */
		    d__1 = t138;
		    t519 = d__1 * d__1;
		    t520 = 1 / t519;
		    t521 = t130 * t520;
/* Computing 2nd power */
		    d__1 = t131;
		    t522 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t522;
		    t523 = d__1 * d__1;
		    t524 = t523 * t131;
		    t525 = 1 / t524;
		    t526 = 1 / t86;
		    t531 = 1 / t134;
		    t534 = 1 / t127;
		    t537 = t525 * -1.853395810515781 * t526 - t515 * 
			    1.28158207914907 * t526 - t531 * 
			    .8223668877838045 * t526 - t534 * 
			    .1603914194192128 * t526;
		    t538 = 1 / t141;
		    t539 = t537 * t538;
		    t540 = t539 * t160;
		    t552 = 1 / t156 / t144;
		    t555 = t480 * .1390394666666667 * t145 - t486 * .12016352 
			    * t149 + t491 * .1459892053333333 * t153 - t497 * 
			    .06791957333333333 * t157 + t504 * 
			    .008483191466666667 * t552;
		    t556 = t142 * t555;
		    t559 = -t400 - t396;
		    t562 = t165 + t395;
		    t565 = t199 * 1.333333333333333 * t559 + t203 * 
			    1.333333333333333 * t562;
		    t567 = t194 * t565 * t213;
		    t568 = t567 * .03799574853701528;
		    t569 = t410 + t414;
		    t571 = t194 * t205 * t569;
		    t572 = t571 * .03799574853701528;
		    t573 = t230 * t565;
		    t574 = t573 * t211;
		    t575 = t574 * 1.923661050931536;
		    t583 = t520 * t537 * t538;
		    t586 = -t182 + t216 + t233 + t164 * (t355 + t376 - t380 - 
			    t394 + t568 + t572 + t438 + t575 - t442 * 
			    7.694644203726145 - t445) + t512 * .03109 - t516 *
			     .001321010150222857 * t142 - t130 * 1. * t583;
		    t598 = t480 * -.05543856 * t247 + t459 * .00270453216 * 
			    t480 - t462 * 3.4198272e-5 * t480 + t465 * 
			    1.98068544e-7 * t480 - t470 * 4.696040448e-10 * 
			    t480;
		    vrhob[i__] = t84 * -1.2407009817988 * t123 - t85 * 
			    .9305257363491 * t509 - t512 * .03109 * t160 + 
			    t516 * .001321010150222857 * t161 + t521 * 1. * 
			    t540 - t130 * .03109 * t556 + t586 * t265 + t240 *
			     t598;
		    t602 = sigmaaa * t20;
		    t605 = t16 * t28;
		    t608 = t26 * t37;
		    t613 = 1 / t4 / t27 / t18;
		    t614 = t34 * t613;
		    t617 = t9 * -.00298886 * t13 + t602 * 1.74462e-4 * t23 - 
			    t605 * 1.43865856e-6 * t31 + t608 * 4.3543808e-9 *
			     t40 - t614 * 4.79940608e-12 * t298;
		    t630 = t9 * -.0521398 * t65 + t602 * .04506132 * t69 - 
			    t605 * .054745952 * t73 + t608 * .02546984 * t77 
			    - t614 * .0031811968 * t344;
		    t631 = t62 * t630;
		    t644 = t9 * .02078946 * t247 - t459 * .00101419956 * t9 + 
			    t462 * 1.2824352e-5 * t9 - t465 * 7.4275704e-8 * 
			    t9 + t470 * 1.761015168e-10 * t9;
		    vsigmaaa[i__] = t5 * -.9305257363491 * t617 - t50 * 
			    .03109 * t631 + t240 * t644;
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t659 = 1 / t84 / t107 / t98;
		    t660 = t114 * t659;
		    t663 = t89 * -.00298886 * t93 + t648 * 1.74462e-4 * t103 
			    - t651 * 1.43865856e-6 * t111 + t654 * 
			    4.3543808e-9 * t120 - t660 * 4.79940608e-12 * 
			    t506;
		    t676 = t89 * -.0521398 * t145 + t648 * .04506132 * t149 - 
			    t651 * .054745952 * t153 + t654 * .02546984 * 
			    t157 - t660 * .0031811968 * t552;
		    t677 = t142 * t676;
		    t690 = t89 * .02078946 * t247 - t459 * .00101419956 * t89 
			    + t462 * 1.2824352e-5 * t89 - t465 * 7.4275704e-8 
			    * t89 + t470 * 1.761015168e-10 * t89;
		    vsigmabb[i__] = t85 * -.9305257363491 * t663 - t130 * 
			    .03109 * t677 + t240 * t690;
		    t696 = sigmaaa / t7 / t17;
		    t699 = t17 * t269;
		    t702 = t16 / t4 / t699;
		    t705 = t243 * t257;
		    t710 = t250 * t262;
		    t715 = t255 * t469;
		    t722 = t260 / t261 / t251;
		    t739 = t26 / t35;
		    t745 = t34 / t7 / t27 / t17;
		    t751 = t292 / t4 / t27 / t699;
/* Computing 2nd power */
		    d__1 = t27;
		    t755 = d__1 * d__1;
		    t758 = t34 * t16 / t755 / t6;
		    t760 = 1 / t39 / t22;
		    t766 = t49 * t312;
		    t774 = 1 / t311 / t58;
/* Computing 2nd power */
		    d__1 = t329;
		    t776 = d__1 * d__1;
		    t784 = 1 / t269;
		    t786 = 1 / t56 / t46;
		    t787 = t784 * t786;
		    t792 = 1 / t17;
		    t813 = -1.544496508763151 / t316 / t46 * t792 + t317 * 
			    3.706791621031562 * t784 - t786 * 
			    .854388052766047 * t792 + t307 * 
			    2.563164158298141 * t784 - .4111834438919023 / 
			    t54 / t46 * t792 + t323 * 1.644733775567609 * 
			    t784 - .05346380647307093 / t47 / t46 * t792 + 
			    t326 * .3207828388384256 * t784;
/* Computing 2nd power */
		    d__1 = t311;
		    t821 = d__1 * d__1;
		    t822 = 1 / t821;
/* Computing 2nd power */
		    d__1 = t61;
		    t824 = d__1 * d__1;
		    t825 = 1 / t824;
		    t841 = 1 / t76 / t68;
		    t848 = t375 * 2.;
		    t850 = t354 * .005495599555936838;
		    t852 = t379 * .001748158859896213;
		    t853 = t393 * 2.249999913366216;
		    t854 = t437 * 3.847322101863073;
		    t856 = t444 * 15.38928840745229;
		    t876 = t408 * t412;
		    t877 = t231 * t876;
		    t879 = t436 * t409;
		    t883 = 1 / t174 / t165 * t210;
		    t884 = t883 * t180;
		    t885 = t884 * .001831866518645613;
/* Computing 2nd power */
		    d__1 = t388;
		    t889 = d__1 * d__1;
		    t892 = t184 * 2.249999913366216 / t381 / t189 * t889 * 
			    t392;
		    t894 = 1 / t208 / t164;
		    t895 = t351 * t894;
		    t896 = t895 * t180;
		    t897 = t896 * .005495599555936838;
		    t898 = t439 * t413;
		    t901 = t436 * 15.38928840745229 * t413;
		    t902 = t206 * t210;
		    t903 = t231 * t902;
		    t904 = t903 * 23.08393261117844;
/* Computing 2nd power */
		    d__1 = t199;
		    t905 = d__1 * d__1;
		    t906 = 1 / t905;
/* Computing 2nd power */
		    d__1 = t397;
		    t907 = d__1 * d__1;
		    t910 = t352 * 2.;
		    t912 = t196 * 2. * t894;
		    t913 = -t910 + t912;
/* Computing 2nd power */
		    d__1 = t203;
		    t916 = d__1 * d__1;
		    t917 = 1 / t916;
/* Computing 2nd power */
		    d__1 = t401;
		    t918 = d__1 * d__1;
		    t924 = t906 * .4444444444444444 * t907 + t199 * 
			    1.333333333333333 * t913 + t917 * 
			    .4444444444444444 * t918 - t203 * 
			    1.333333333333333 * t913;
		    t928 = t902 * 12.;
		    t929 = t876 * 32.;
		    t932 = t207 / t209 / t208;
		    t933 = t932 * 20.;
		    t940 = t389 * t390 * t404 * t213;
		    t943 = t389 * t391 * t415;
/* Computing 2nd power */
		    d__1 = t381;
		    t945 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t192;
		    t949 = d__1 * d__1;
		    t954 = t184 * 33.30964519106732 / t945 * t889 / t949 * 
			    t205 * t213;
		    t955 = t439 * t409;
		    t957 = t877 * -61.55715362980916 + t879 * 
			    15.38928840745229 + t885 + t892 - t897 - t898 * 
			    15.38928840745229 - t901 + t904 + t194 * 
			    .03799574853701528 * t924 * t213 + t194 * 
			    .03799574853701528 * t205 * (-t928 + t929 - t933) 
			    - t940 * 2.249999913366216 - t943 * 
			    2.249999913366216 - t954 + t955 * 
			    15.38928840745229;
		    t959 = t435 * t404 * t211;
		    t963 = 1 / t361 / t165 * t210;
		    t965 = t362 * t894;
		    t971 = 1 / t172 / t165 * t210;
		    t973 = t366 * t894;
		    t977 = 1 / t166 / t165 * t210;
		    t979 = t369 * t894;
		    t984 = t383 * 1.124999956683108 * (t963 * 
			    -1.132974264373283 + t965 * 2.71913823449588 - 
			    t883 * .4994648585728036 + t895 * 
			    1.498394575718411 - t971 * .1075243117819161 + 
			    t973 * .4300972471276643 - t977 * 
			    .04247805766949639 + t979 * .2548683460169784) * 
			    t392;
/* Computing 2nd power */
		    d__1 = t372;
		    t991 = d__1 * d__1;
		    t993 = t168 / t356 / t176 * t991 * t373;
		    t994 = t993 * 2.;
		    t1005 = t358 * (t963 * -.8309097827459833 + t965 * 
			    1.99418347859036 - t883 * .4945709824779306 + 
			    t895 * 1.483712947433792 - t971 * 
			    .2001071587498409 + t973 * .8004286349993634 - 
			    t977 * .04215565168327908 + t979 * 
			    .2529339100996745) * t373;
		    t1006 = t1005 * 1.;
		    t1009 = t353 * t357 * t372 * t373;
		    t1010 = t1009 * .08837926660346786;
/* Computing 2nd power */
		    d__1 = t356;
		    t1011 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t179;
		    t1014 = d__1 * d__1;
		    t1017 = t168 / t1011 * t991 / t1014;
		    t1018 = t1017 * 16.0818243221511;
/* Computing 2nd power */
		    d__1 = t428;
		    t1030 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t421;
		    t1046 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t226;
		    t1049 = d__1 * d__1;
		    t1060 = t883 * 8.806734334819047e-4 * t227 - t895 * 
			    .002642020300445714 * t227 - t353 * 
			    .08497974591333914 * t422 * t428 * t429 - t218 * 
			    2. / t421 / t223 * t1030 * t429 + t423 * 1. * (
			    t963 * -1.544496508763151 + t965 * 
			    3.706791621031562 - t883 * .854388052766047 + 
			    t895 * 2.563164158298141 - t971 * 
			    .4111834438919023 + t973 * 1.644733775567609 - 
			    t977 * .05346380647307093 + t979 * 
			    .3207828388384256) * t429 + t218 * 
			    32.1646831778707 / t1046 * t1030 / t1049 - t884 * 
			    .001831866518645613 + t896 * .005495599555936838 
			    + t1009 * .08837926660346786 + t993 * 2. - t1005 *
			     1. - t1017 * 16.0818243221511;
		    t1063 = t1060 * 1.923661050931536 * t205 * t211;
		    t1065 = t231 * 38.47322101863073 * t932;
		    t1068 = t353 * t193 * t404 * t213;
		    t1071 = t353 * t377 * t415;
		    t1074 = t883 * 5.827196199654043e-4 * t378;
		    t1082 = t353 * .05176049209143758 * t382 * t388 * t390 * 
			    t214;
		    t1084 = t895 * .001748158859896213 * t378;
		    t1085 = t959 * 3.847322101863073 - t984 + t194 * 
			    .07599149707403056 * t404 * t415 - t994 + t1006 - 
			    t1010 + t1018 + t1063 + t1065 - t1068 * 
			    .001748158859896213 - t1071 * .001748158859896213 
			    - t1074 + t230 * 1.923661050931536 * t924 * t211 
			    + t1082 + t1084;
		    t1088 = t848 + t406 * .07599149707403056 + t850 + t417 * 
			    .07599149707403056 - t852 - t853 + t854 + t442 * 
			    15.38928840745229 - t856 + t440 * 
			    3.847322101863073 - t50 * 32.1646831778707 * t822 
			    * t776 * t825 - t50 * 1. * t312 * t813 * t330 - 
			    t787 * 8.806734334819047e-4 * t62 + t308 * 
			    .08497974591333914 * t452 + t50 * 2. * t774 * 
			    t776 * t330 - t766 * 2. * t331 + t164 * (t957 + 
			    t1085);
		    s1 = t455 * 2. * t473 + t240 * (t696 * .20327472 * t247 - 
			    t702 * .00404955136 * t252 + t705 * 
			    1.3446790656e-4 * t702 - t459 * .00991661792 * 
			    t696 - t710 * 1.613032704e-6 * t702 + t462 * 
			    1.25393664e-4 * t696 + t715 * 8.8427483136e-9 * 
			    t702 - t465 * 7.26251328e-7 * t696 - t722 * 
			    1.8784161792e-11 * t702 + t470 * 1.7218814976e-9 *
			     t696) - .4135669939329333 / t7 * t43 - t4 * 
			    2.4814019635976 * t301 - t5 * .9305257363491 * (
			    t696 * -.02922440888888889 * t13 + t702 * 
			    .003031485795555556 * t23 - t739 * 
			    4.445275477333333e-5 * t31 + t745 * 
			    2.582351553422222e-7 * t40 - t751 * 
			    6.788757367466667e-10 * t298 + t758 * 
			    6.825821980444444e-13 * t760) + t766 * 2. * t332 
			    - t304 * .06218 * t347 + t308 * 
			    .002642020300445714 * t348;
		    v2rhoa2[i__] = s1 - t50 * 2. * t774 * t776 * t330 * t80 - 
			    t308 * .08497974591333914 * t312 * t332 + t787 * 
			    8.806734334819047e-4 * t81 + t313 * 1. * t813 * 
			    t330 * t80 + t313 * 2. * t331 * t347 + t50 * 
			    32.1646831778707 * t822 * t776 * t825 * t80 - t50 
			    * .03109 * t62 * (t696 * -.5098113777777778 * t65 
			    + t702 * .8351900088888889 * t69 - t739 * 
			    1.442077269333333 * t73 + t745 * 
			    1.025977750755556 * t77 - t751 * .2664875008 * 
			    t344 + t758 * .02262184391111111 * t841) + t1088 *
			     t265;
		    t1090 = t129 * t520;
		    t1104 = sigmabb / t87 / t97;
		    t1107 = t97 * t477;
		    t1110 = t96 / t84 / t1107;
		    t1114 = t106 / t115;
		    t1120 = t114 / t87 / t107 / t97;
		    t1126 = t500 / t84 / t107 / t1107;
/* Computing 2nd power */
		    d__1 = t107;
		    t1130 = d__1 * d__1;
		    t1133 = t114 * t96 / t1130 / t86;
		    t1135 = 1 / t119 / t102;
		    t1141 = 1 / t477;
		    t1143 = 1 / t136 / t126;
		    t1144 = t1141 * t1143;
		    t1151 = 1 / t519 / t138;
/* Computing 2nd power */
		    d__1 = t537;
		    t1153 = d__1 * d__1;
		    t1169 = 1 / t156 / t148;
		    t1178 = 1 / t97;
		    t1199 = -1.544496508763151 / t524 / t126 * t1178 + t525 * 
			    3.706791621031562 * t1141 - t1143 * 
			    .854388052766047 * t1178 + t515 * 
			    2.563164158298141 * t1141 - .4111834438919023 / 
			    t134 / t126 * t1178 + t531 * 1.644733775567609 * 
			    t1141 - .05346380647307093 / t127 / t126 * t1178 
			    + t534 * .3207828388384256 * t1141;
/* Computing 2nd power */
		    d__1 = t519;
		    t1207 = d__1 * d__1;
		    t1208 = 1 / t1207;
/* Computing 2nd power */
		    d__1 = t141;
		    t1210 = d__1 * d__1;
		    t1211 = 1 / t1210;
/* Computing 2nd power */
		    d__1 = t559;
		    t1222 = d__1 * d__1;
		    t1225 = t910 + t912;
/* Computing 2nd power */
		    d__1 = t562;
		    t1228 = d__1 * d__1;
		    t1234 = t906 * .4444444444444444 * t1222 + t199 * 
			    1.333333333333333 * t1225 + t917 * 
			    .4444444444444444 * t1228 - t203 * 
			    1.333333333333333 * t1225;
		    t1245 = t877 * 61.55715362980916 - t879 * 
			    15.38928840745229 + t885 + t892 - t897 - t901 + 
			    t904 - t954 - t984 - t994 + t1006 + t194 * 
			    .03799574853701528 * t1234 * t213 + t194 * 
			    .07599149707403056 * t565 * t569 + t194 * 
			    .03799574853701528 * t205 * (-t928 - t929 - t933);
		    t1247 = t353 * t377 * t569;
		    t1251 = t353 * t193 * t565 * t213;
		    t1255 = t389 * t390 * t565 * t213;
		    t1258 = t389 * t391 * t569;
		    t1261 = t435 * t565 * t211;
		    t1266 = t573 * t413;
		    t1268 = t573 * t409;
		    t1270 = -t1010 - t1247 * .001748158859896213 + t1018 - 
			    t1251 * .001748158859896213 - t1255 * 
			    2.249999913366216 - t1258 * 2.249999913366216 + 
			    t1063 + t1065 + t1261 * 3.847322101863073 + t230 *
			     1.923661050931536 * t1234 * t211 - t1266 * 
			    15.38928840745229 - t1268 * 15.38928840745229 - 
			    t1074 + t1082 + t1084;
		    t1291 = t848 + t850 - t852 - t853 + t854 - t442 * 
			    15.38928840745229 - t856 + t574 * 
			    3.847322101863073 + t567 * .07599149707403056 + 
			    t571 * .07599149707403056 + t164 * (t1245 + t1270)
			     - t1090 * 2. * t539 - t1144 * 
			    8.806734334819047e-4 * t142 + t516 * 
			    .08497974591333914 * t583 + t130 * 2. * t1151 * 
			    t1153 * t538 - t130 * 1. * t520 * t1199 * t538 - 
			    t130 * 32.1646831778707 * t1208 * t1153 * t1211;
		    s1 = t1090 * 2. * t540 + t516 * .002642020300445714 * 
			    t556 - .4135669939329333 / t87 * t123 - t84 * 
			    2.4814019635976 * t509 - t512 * .06218 * t555 - 
			    t85 * .9305257363491 * (t1104 * 
			    -.02922440888888889 * t93 + t1110 * 
			    .003031485795555556 * t103 - t1114 * 
			    4.445275477333333e-5 * t111 + t1120 * 
			    2.582351553422222e-7 * t120 - t1126 * 
			    6.788757367466667e-10 * t506 + t1133 * 
			    6.825821980444444e-13 * t1135) + t1144 * 
			    8.806734334819047e-4 * t161 - t516 * 
			    .08497974591333914 * t520 * t540;
		    s2 = s1 - t130 * 2. * t1151 * t1153 * t538 * t160 - t130 *
			     .03109 * t142 * (t1104 * -.5098113777777778 * 
			    t145 + t1110 * .8351900088888889 * t149 - t1114 * 
			    1.442077269333333 * t153 + t1120 * 
			    1.025977750755556 * t157 - t1126 * .2664875008 * 
			    t552 + t1133 * .02262184391111111 * t1169) + t521 
			    * 1. * t1199 * t538 * t160;
		    v2rhob2[i__] = s2 + t521 * 2. * t539 * t555 + t130 * 
			    32.1646831778707 * t1208 * t1153 * t1211 * t160 + 
			    t1291 * t265 + t586 * 2. * t598 + t240 * (t1104 * 
			    .20327472 * t247 - t1110 * .00404955136 * t252 + 
			    t705 * 1.3446790656e-4 * t1110 - t459 * 
			    .00991661792 * t1104 - t710 * 1.613032704e-6 * 
			    t1110 + t462 * 1.25393664e-4 * t1104 + t715 * 
			    8.8427483136e-9 * t1110 - t465 * 7.26251328e-7 * 
			    t1104 - t722 * 1.8784161792e-11 * t1110 + t470 * 
			    1.7218814976e-9 * t1104);
		    t1330 = t906 * .4444444444444444 * t397 * t559 + t199 * 
			    2.666666666666667 * t196 * t894 + t917 * 
			    .4444444444444444 * t401 * t562 - t203 * 
			    2.666666666666667 * t196 * t894;
		    t1346 = t885 + t892 - t897 - t898 * 7.694644203726145 - 
			    t901 + t194 * .03799574853701528 * t1330 * t213 + 
			    t194 * .03799574853701528 * t205 * (t928 - t933) 
			    - t903 * 23.08393261117844 - t940 * 
			    1.124999956683108 - t943 * 1.124999956683108 + 
			    t194 * .03799574853701528 * t404 * t569 - t954 - 
			    t955 * 7.694644203726145 + t959 * 
			    1.923661050931536 - t984 - t994 + t1006;
		    t1362 = -t1010 - t1247 * 8.740794299481065e-4 + t1018 - 
			    t1251 * 8.740794299481065e-4 - t1255 * 
			    1.124999956683108 - t1258 * 1.124999956683108 + 
			    t1063 + t1065 + t1261 * 1.923661050931536 - t1266 
			    * 7.694644203726145 - t1068 * 
			    8.740794299481065e-4 - t1071 * 
			    8.740794299481065e-4 + t1268 * 7.694644203726145 
			    - t1074 + t1082 + t1084 + t194 * 
			    .03799574853701528 * t565 * t415 + t230 * 
			    1.923661050931536 * t1330 * t211;
		    t1365 = t850 + t848 - t852 - t853 + t568 + t572 + t854 + 
			    t575 - t856 + t407 + t418 + t441 + t164 * (t1346 
			    + t1362);
		    t1375 = t271 * sigmabb * t479;
		    v2rhoab[i__] = t1365 * t265 + t455 * t598 + t586 * t473 + 
			    t240 * (t272 * -.00404955136 * t252 * sigmabb * 
			    t479 + t705 * 1.3446790656e-4 * sigmaaa * t1375 - 
			    t710 * 1.613032704e-6 * sigmaaa * t1375 + t715 * 
			    8.8427483136e-9 * sigmaaa * t1375 - t722 * 
			    1.8784161792e-11 * sigmaaa * t1375);
		    t1393 = sigmaaa * t277;
		    t1396 = t16 * t282;
		    t1399 = t26 * t288;
		    t1402 = t34 * t295;
		    t1407 = t292 / t755 / rhoa;
		    s1 = t4 * -1.2407009817988 * t617 - t5 * .9305257363491 * 
			    (t271 * .007970293333333333 * t13 - t1393 * 
			    9.623451733333333e-4 * t23 + t1396 * 
			    1.523112448e-5 * t31 - t1399 * 
			    9.248380245333333e-8 * t40 + t1402 * 
			    2.497789952e-10 * t298 - t1407 * 
			    2.559683242666667e-13 * t760) - t304 * .03109 * 
			    t630 + t308 * .001321010150222857 * t631;
		    v2rhoasigmaaa[i__] = s1 + t313 * 1. * t331 * t630 - t50 * 
			    .03109 * t62 * (t271 * .1390394666666667 * t65 - 
			    t1393 * .2681349333333333 * t69 + t1396 * 
			    .486033024 * t73 - t1399 * .3592718165333333 * 
			    t77 + t1402 * .096751616 * t344 - t1407 * 
			    .008483191466666667 * t841) + t455 * t644 + t240 *
			     (t271 * -.05543856 * t247 + t1393 * .00151858176 
			    * t252 - t705 * 5.042546496e-5 * t1393 + t459 * 
			    .00270453216 * t271 + t710 * 6.04887264e-7 * 
			    t1393 - t462 * 3.4198272e-5 * t271 - t715 * 
			    3.3160306176e-9 * t1393 + t465 * 1.98068544e-7 * 
			    t271 + t722 * 7.044060672e-12 * t1393 - t470 * 
			    4.696040448e-10 * t271);
		    v2rhoasigmaab[i__] = 0.;
		    t1463 = t272 * t89;
		    v2rhoasigmabb[i__] = t455 * t690 + t240 * (t272 * 
			    .00151858176 * t252 * t89 - t705 * 5.042546496e-5 
			    * t1463 + t710 * 6.04887264e-7 * t1463 - t715 * 
			    3.3160306176e-9 * t1463 + t722 * 7.044060672e-12 *
			     t1463);
		    t1475 = t252 * t9;
		    t1478 = t480 * t9;
		    v2rhobsigmaaa[i__] = t586 * t644 + t240 * (t480 * 
			    .00151858176 * t1475 - t705 * 5.042546496e-5 * 
			    t1478 + t710 * 6.04887264e-7 * t1478 - t715 * 
			    3.3160306176e-9 * t1478 + t722 * 7.044060672e-12 *
			     t1478);
		    v2rhobsigmaab[i__] = 0.;
		    t1493 = sigmabb * t485;
		    t1496 = t96 * t490;
		    t1499 = t106 * t496;
		    t1502 = t114 * t503;
		    t1507 = t500 / t1130 / rhob;
		    s1 = t84 * -1.2407009817988 * t663 - t85 * .9305257363491 
			    * (t479 * .007970293333333333 * t93 - t1493 * 
			    9.623451733333333e-4 * t103 + t1496 * 
			    1.523112448e-5 * t111 - t1499 * 
			    9.248380245333333e-8 * t120 + t1502 * 
			    2.497789952e-10 * t506 - t1507 * 
			    2.559683242666667e-13 * t1135) - t512 * .03109 * 
			    t676 + t516 * .001321010150222857 * t677;
		    v2rhobsigmabb[i__] = s1 + t521 * 1. * t539 * t676 - t130 *
			     .03109 * t142 * (t479 * .1390394666666667 * t145 
			    - t1493 * .2681349333333333 * t149 + t1496 * 
			    .486033024 * t153 - t1499 * .3592718165333333 * 
			    t157 + t1502 * .096751616 * t552 - t1507 * 
			    .008483191466666667 * t1169) + t586 * t690 + t240 
			    * (t479 * -.05543856 * t247 + t1493 * 
			    .00151858176 * t252 - t705 * 5.042546496e-5 * 
			    t1493 + t459 * .00270453216 * t479 + t710 * 
			    6.04887264e-7 * t1493 - t462 * 3.4198272e-5 * 
			    t479 - t715 * 3.3160306176e-9 * t1493 + t465 * 
			    1.98068544e-7 * t479 + t722 * 7.044060672e-12 * 
			    t1493 - t470 * 4.696040448e-10 * t479);
		    t1561 = sigmaaa * t28;
		    t1564 = t16 * t37;
		    t1567 = t26 * t613;
		    t1571 = t34 / t755;
		    v2sigmaaa2[i__] = t5 * -.9305257363491 * (t20 * 
			    1.8641744e-4 * t23 - t1561 * 4.27301312e-6 * t31 
			    + t1564 * 3.032704512e-8 * t40 - t1567 * 
			    8.886771712e-11 * t298 + t1571 * 9.59881216e-14 * 
			    t760) - t50 * .03109 * t62 * (t20 * .05548928 * 
			    t69 - t1561 * .127516432 * t73 + t1564 * 
			    .1092570912 * t77 - t1567 * .0331006592 * t344 + 
			    t1571 * .0031811968 * t841) + t240 * (t20 * 
			    -5.6946816e-4 * t252 + t705 * 1.890954936e-5 * 
			    t20 - t710 * 2.26832724e-7 * t20 + t715 * 
			    1.2435114816e-9 * t20 - t722 * 2.641522752e-12 * 
			    t20);
		    v2sigmaaaab[i__] = 0.;
		    t1605 = t9 * t89;
		    v2sigmaaabb[i__] = t240 * (t1475 * -5.6946816e-4 * t89 + 
			    t705 * 1.890954936e-5 * t1605 - t710 * 
			    2.26832724e-7 * t1605 + t715 * 1.2435114816e-9 * 
			    t1605 - t722 * 2.641522752e-12 * t1605);
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t1617 = sigmabb * t108;
		    t1620 = t96 * t117;
		    t1623 = t106 * t659;
		    t1627 = t114 / t1130;
		    v2sigmabb2[i__] = t85 * -.9305257363491 * (t100 * 
			    1.8641744e-4 * t103 - t1617 * 4.27301312e-6 * 
			    t111 + t1620 * 3.032704512e-8 * t120 - t1623 * 
			    8.886771712e-11 * t506 + t1627 * 9.59881216e-14 * 
			    t1135) - t130 * .03109 * t142 * (t100 * .05548928 
			    * t149 - t1617 * .127516432 * t153 + t1620 * 
			    .1092570912 * t157 - t1623 * .0331006592 * t552 + 
			    t1627 * .0031811968 * t1169) + t240 * (t100 * 
			    -5.6946816e-4 * t252 + t705 * 1.890954936e-5 * 
			    t100 - t710 * 2.26832724e-7 * t100 + t715 * 
			    1.2435114816e-9 * t100 - t722 * 2.641522752e-12 * 
			    t100);
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
} /* uks_xc_hcth120__ */

/* Subroutine */ int rks_xc_hcth120__(integer *ideriv, integer *npt, 
	doublereal *rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *
	vrhoa, doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *
	v2rhoasigmaaa, doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t4, t5, t7, t8, t10, t20, t11, t21, t14,
	     t15, t25, t32, t27, t19, t36, t37, t44, t45, t48, t49, t52, t54, 
	    t60, t62, t66, t74, t93, t16, t18, t24, t26, t29, t35, t38, t41, 
	    t47, t56, t59, t63, t67, t71, t75, t78, t79, t83, t89, t92, t98, 
	    t33, t100, t101, t120, t112, t104, t105, t113, t116, t123, t109, 
	    t126, t129, t134, t140, t147, t149, t155, t158, t159, t162, t163, 
	    t165, t166, t169, t170, t172, t175, t178, t180, t181, t195, t206, 
	    t238, t245, t248, t251, t257, t122, t128, t133, t139, t143, t146, 
	    t152, t164, t167, t168, t174, t177, t182, t183, t198, t199, t204, 
	    t207, t208, t213, t214, t216, t224, t227, t241, t256, t260, t273, 
	    t274, t287, t294, rho, t297, t300, t301, t304, t305, t310, t311, 
	    t316, t317, t320, t323, t325, t327, t334, t348, t361, t362, t364, 
	    t365, t368, t369, t371, t375, t377, t381, t383, t385, t390, t397, 
	    t399, t416, t417, t419, t420, t424, t425, t427, t429, t440, t442, 
	    t447, t450, t462, t463, t466, t470, t474, t492, t522, t525, t528, 
	    t531, t536, t569, t571, t573, t575, t578, t589, t592, t595, t599, 
	    sigma;


/*     A.D. Boese, N.L. Doltsinis, N.C. Handy, M. Sprik */
/*     New generalized gradient approximation functionals */
/*     J. Chem. Phys. 112 (2000) 1670-1678. */


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
/* Computing 2nd power */
		d__1 = t15;
		t25 = d__1 * d__1;
		t27 = t14 * sigma / t25;
/* Computing 2nd power */
		d__1 = t14;
		t32 = d__1 * d__1;
		t36 = t32 / t5 / t25 / t4;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
		t44 = 1 / rho;
		t45 = pow_dd(&t44, &c_b2);
		t48 = rho * (t45 * .1606016560364007 + 1.);
		t49 = pow_dd(&t44, &c_b4);
		t52 = sqrt(t44);
/* Computing 2nd power */
		d__1 = t45;
		t54 = d__1 * d__1;
		t60 = log(32.1646831778707 / (t49 * 12.48219874679732 + t45 * 
			4.844076716063854 + t52 * 2.326004811900819 + t54 * 
			.3819082618690966) + 1.);
		t62 = t8 * .3174802103936399 + 1.;
/* Computing 2nd power */
		d__1 = t62;
		t66 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t66;
		t74 = d__1 * d__1;
		t93 = log(16.0818243221511 / (t49 * 5.98255043577108 + t45 * 
			2.225569421150687 + t52 * .8004286349993634 + t54 * 
			.1897004325747559) + 1.);
		t100 = t8 * .009524406311809197 + 1.;
/* Computing 2nd power */
		d__1 = t100;
		t104 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t104;
		t112 = d__1 * d__1;
		zk[i__] = t2 * -.7385587663820224 * rho * (1.09163 - t8 * 
			.004744519508185673 / t10 + t19 * 
			2.047454356900042e-4 / t20 - t27 * 1.05150976e-6 / 
			t20 / t10 + t36 * 1.904645565053643e-9 / t37) - t48 * 
			.03109 * t60 * (.489508 - t8 * .08276677336941153 / 
			t62 + t19 * .0436352992925871 / t66 - t27 * .06375904 
			/ t66 / t62 + t36 * .02524917573418935 / t74) + (rho *
			 -.062182 * (t45 * .1325688999052018 + 1.) * t93 + 
			t48 * .03109 * t60) * (t8 * .06600242134770161 / t100 
			+ .51473 - t19 * .002241305809636867 / t104 + t27 * 
			1.99668672e-5 / t104 / t100 - t36 * 
			9.318124434050518e-8 / t112);
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
		t16 = t15 * rho;
		t18 = 1 / t2 / t16;
		t19 = t14 * t18;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		t24 = t14 * sigma;
/* Computing 2nd power */
		d__1 = t15;
		t25 = d__1 * d__1;
		t26 = 1 / t25;
		t27 = t24 * t26;
		t29 = 1 / t20 / t10;
/* Computing 2nd power */
		d__1 = t14;
		t32 = d__1 * d__1;
		t35 = 1 / t5 / t25 / t4;
		t36 = t32 * t35;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
		t38 = 1 / t37;
		t41 = 1.09163 - t8 * .004744519508185673 * t11 + t19 * 
			2.047454356900042e-4 * t21 - t27 * 1.05150976e-6 * 
			t29 + t36 * 1.904645565053643e-9 * t38;
		t44 = 1 / rho;
		t45 = pow_dd(&t44, &c_b2);
		t47 = t45 * .1606016560364007 + 1.;
		t48 = rho * t47;
		t49 = pow_dd(&t44, &c_b4);
		t52 = sqrt(t44);
/* Computing 2nd power */
		d__1 = t45;
		t54 = d__1 * d__1;
		t56 = t49 * 12.48219874679732 + t45 * 4.844076716063854 + t52 
			* 2.326004811900819 + t54 * .3819082618690966;
		t59 = 32.1646831778707 / t56 + 1.;
		t60 = log(t59);
		t62 = t8 * .3174802103936399 + 1.;
		t63 = 1 / t62;
/* Computing 2nd power */
		d__1 = t62;
		t66 = d__1 * d__1;
		t67 = 1 / t66;
		t71 = 1 / t66 / t62;
/* Computing 2nd power */
		d__1 = t66;
		t74 = d__1 * d__1;
		t75 = 1 / t74;
		t78 = .489508 - t8 * .08276677336941153 * t63 + t19 * 
			.0436352992925871 * t67 - t27 * .06375904 * t71 + t36 
			* .02524917573418935 * t75;
		t79 = t60 * t78;
		t83 = t45 * .1325688999052018 + 1.;
		t89 = t49 * 5.98255043577108 + t45 * 2.225569421150687 + t52 *
			 .8004286349993634 + t54 * .1897004325747559;
		t92 = 16.0818243221511 / t89 + 1.;
		t93 = log(t92);
		t98 = rho * -.062182 * t83 * t93 + t48 * .03109 * t60;
		t100 = t8 * .009524406311809197 + 1.;
		t101 = 1 / t100;
/* Computing 2nd power */
		d__1 = t100;
		t104 = d__1 * d__1;
		t105 = 1 / t104;
		t109 = 1 / t104 / t100;
/* Computing 2nd power */
		d__1 = t104;
		t112 = d__1 * d__1;
		t113 = 1 / t112;
		t116 = t8 * .06600242134770161 * t101 + .51473 - t19 * 
			.002241305809636867 * t105 + t27 * 1.99668672e-5 * 
			t109 - t36 * 9.318124434050518e-8 * t113;
		zk[i__] = t3 * -.7385587663820224 * t41 - t48 * .03109 * t79 
			+ t98 * t116;
		t120 = t4 * rho;
		t123 = sigma / t5 / t120;
		t126 = t15 * t4;
		t129 = t14 / t2 / t126;
		t134 = t24 / t25 / rho;
		t140 = t32 / t5 / t25 / t120;
		t147 = t32 * sigma / t2 / t25 / t126;
		t149 = 1 / t37 / t10;
		t155 = t47 * t60;
		t158 = 1 / t54;
		t159 = t44 * t158;
/* Computing 2nd power */
		d__1 = t56;
		t162 = d__1 * d__1;
		t163 = 1 / t162;
/* Computing 2nd power */
		d__1 = t49;
		t165 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t165;
		t166 = d__1 * d__1;
		t169 = 1 / t4;
		t170 = 1 / t166 / t49 * t169;
		t172 = t158 * t169;
		t175 = 1 / t52 * t169;
		t178 = 1 / t45 * t169;
		t180 = t170 * -4.160732915599108 - t172 * 3.229384477375903 - 
			t175 * 2.326004811900819 - t178 * .5092110158254621;
		t181 = 1 / t59;
		t195 = 1 / t74 / t62;
/* Computing 2nd power */
		d__1 = t89;
		t206 = d__1 * d__1;
		t238 = 1 / t112 / t100;
		s1 = t2 * -.9847450218426965 * t41 - t3 * .3692793831910112 * 
			(t123 * .02530410404365692 * t11 - t129 * 
			.002344622359538767 * t21 + t134 * 
			3.069138261333333e-5 * t29 - t140 * 
			1.474591714685894e-7 * t38 + t147 * 
			2.57999903879912e-10 * t149) - t155 * .03109 * t78 + 
			t159 * .001664368495390566 * t79;
		vrhoa[i__] = s1 + t48 * .5 * t163 * t180 * t181 * t78 - t48 * 
			.015545 * t60 * (t123 * .4414227913035281 * t63 - 
			t129 * .6055861931098544 * t67 + t134 * 
			1.167913642666667 * t71 - t140 * .8625248172685168 * 
			t75 + t147 * .1710104239862703 * t195) + (t83 * 
			-.062182 * t93 + rho * (t172 * .002747799777968419 * 
			t93 + t83 * 1. / t206 * (t170 * -.99709173929518 - 
			t172 * .7418564737168958 - t175 * .4002143174996817 - 
			t178 * .1264669550498372) / t92) + t155 * .03109 - 
			t159 * .001664368495390566 * t60 - t48 * .5 * t163 * 
			t180 * t181) * t116 + t98 * (t123 * 
			-.1760064569272043 * t101 + t129 * .0136299879940066 *
			 t105 - t134 * 2.73586176e-4 * t109 + t140 * 
			2.515313720859277e-6 * t113 - t147 * 
			9.466624338548721e-9 * t238);
		t245 = sigma * t18;
		t248 = t14 * t26;
		t251 = t24 * t35;
		t257 = t32 / t2 / t25 / t16;
		vsigmaaa[i__] = t3 * -.7385587663820224 * (t7 * 
			-.01897807803274269 * t11 + t245 * 
			.001758466769654075 * t21 - t248 * 2.301853696e-5 * 
			t29 + t251 * 1.105943786014421e-7 * t38 - t257 * 
			1.93499927909934e-10 * t149) - t48 * .03109 * t60 * (
			t7 * -.3310670934776461 * t63 + t245 * 
			.4541896448323908 * t67 - t248 * .875935232 * t71 + 
			t251 * .6468936129513876 * t75 - t257 * 
			.1282578179897027 * t195) + t98 * 2. * (t7 * 
			.1320048426954032 * t101 - t245 * .01022249099550495 *
			 t105 + t248 * 2.05189632e-4 * t109 - t251 * 
			1.886485290644458e-6 * t113 + t257 * 
			7.099968253911541e-9 * t238);
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
		t16 = t15 * rho;
		t18 = 1 / t2 / t16;
		t19 = t14 * t18;
/* Computing 2nd power */
		d__1 = t10;
		t20 = d__1 * d__1;
		t21 = 1 / t20;
		t24 = t14 * sigma;
/* Computing 2nd power */
		d__1 = t15;
		t25 = d__1 * d__1;
		t26 = 1 / t25;
		t27 = t24 * t26;
		t29 = 1 / t20 / t10;
/* Computing 2nd power */
		d__1 = t14;
		t32 = d__1 * d__1;
		t33 = t25 * t4;
		t35 = 1 / t5 / t33;
		t36 = t32 * t35;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
		t38 = 1 / t37;
		t41 = 1.09163 - t8 * .004744519508185673 * t11 + t19 * 
			2.047454356900042e-4 * t21 - t27 * 1.05150976e-6 * 
			t29 + t36 * 1.904645565053643e-9 * t38;
		t44 = 1 / rho;
		t45 = pow_dd(&t44, &c_b2);
		t47 = t45 * .1606016560364007 + 1.;
		t48 = rho * t47;
		t49 = pow_dd(&t44, &c_b4);
		t52 = sqrt(t44);
/* Computing 2nd power */
		d__1 = t45;
		t54 = d__1 * d__1;
		t56 = t49 * 12.48219874679732 + t45 * 4.844076716063854 + t52 
			* 2.326004811900819 + t54 * .3819082618690966;
		t59 = 32.1646831778707 / t56 + 1.;
		t60 = log(t59);
		t62 = t8 * .3174802103936399 + 1.;
		t63 = 1 / t62;
/* Computing 2nd power */
		d__1 = t62;
		t66 = d__1 * d__1;
		t67 = 1 / t66;
		t71 = 1 / t66 / t62;
/* Computing 2nd power */
		d__1 = t66;
		t74 = d__1 * d__1;
		t75 = 1 / t74;
		t78 = .489508 - t8 * .08276677336941153 * t63 + t19 * 
			.0436352992925871 * t67 - t27 * .06375904 * t71 + t36 
			* .02524917573418935 * t75;
		t79 = t60 * t78;
		t83 = t45 * .1325688999052018 + 1.;
		t89 = t49 * 5.98255043577108 + t45 * 2.225569421150687 + t52 *
			 .8004286349993634 + t54 * .1897004325747559;
		t92 = 16.0818243221511 / t89 + 1.;
		t93 = log(t92);
		t98 = rho * -.062182 * t83 * t93 + t48 * .03109 * t60;
		t100 = t8 * .009524406311809197 + 1.;
		t101 = 1 / t100;
/* Computing 2nd power */
		d__1 = t100;
		t104 = d__1 * d__1;
		t105 = 1 / t104;
		t109 = 1 / t104 / t100;
/* Computing 2nd power */
		d__1 = t104;
		t112 = d__1 * d__1;
		t113 = 1 / t112;
		t116 = t8 * .06600242134770161 * t101 + .51473 - t19 * 
			.002241305809636867 * t105 + t27 * 1.99668672e-5 * 
			t109 - t36 * 9.318124434050518e-8 * t113;
		zk[i__] = t3 * -.7385587663820224 * t41 - t48 * .03109 * t79 
			+ t98 * t116;
		t120 = t4 * rho;
		t122 = 1 / t5 / t120;
		t123 = sigma * t122;
		t126 = t15 * t4;
		t128 = 1 / t2 / t126;
		t129 = t14 * t128;
		t133 = 1 / t25 / rho;
		t134 = t24 * t133;
		t139 = 1 / t5 / t25 / t120;
		t140 = t32 * t139;
		t143 = t32 * sigma;
		t146 = 1 / t2 / t25 / t126;
		t147 = t143 * t146;
		t149 = 1 / t37 / t10;
		t152 = t123 * .02530410404365692 * t11 - t129 * 
			.002344622359538767 * t21 + t134 * 
			3.069138261333333e-5 * t29 - t140 * 
			1.474591714685894e-7 * t38 + t147 * 
			2.57999903879912e-10 * t149;
		t155 = t47 * t60;
		t158 = 1 / t54;
		t159 = t44 * t158;
/* Computing 2nd power */
		d__1 = t56;
		t162 = d__1 * d__1;
		t163 = 1 / t162;
		t164 = t48 * t163;
/* Computing 2nd power */
		d__1 = t49;
		t165 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t165;
		t166 = d__1 * d__1;
		t167 = t166 * t49;
		t168 = 1 / t167;
		t169 = 1 / t4;
		t170 = t168 * t169;
		t172 = t158 * t169;
		t174 = 1 / t52;
		t175 = t174 * t169;
		t177 = 1 / t45;
		t178 = t177 * t169;
		t180 = t170 * -4.160732915599108 - t172 * 3.229384477375903 - 
			t175 * 2.326004811900819 - t178 * .5092110158254621;
		t181 = 1 / t59;
		t182 = t180 * t181;
		t183 = t182 * t78;
		t195 = 1 / t74 / t62;
		t198 = t123 * .4414227913035281 * t63 - t129 * 
			.6055861931098544 * t67 + t134 * 1.167913642666667 * 
			t71 - t140 * .8625248172685168 * t75 + t147 * 
			.1710104239862703 * t195;
		t199 = t60 * t198;
		t204 = t172 * t93;
/* Computing 2nd power */
		d__1 = t89;
		t206 = d__1 * d__1;
		t207 = 1 / t206;
		t208 = t83 * t207;
		t213 = t170 * -.99709173929518 - t172 * .7418564737168958 - 
			t175 * .4002143174996817 - t178 * .1264669550498372;
		t214 = 1 / t92;
		t216 = t208 * t213 * t214;
		t224 = t163 * t180 * t181;
		t227 = t83 * -.062182 * t93 + rho * (t204 * 
			.002747799777968419 + t216 * 1.) + t155 * .03109 - 
			t159 * .001664368495390566 * t60 - t48 * .5 * t224;
		t238 = 1 / t112 / t100;
		t241 = t123 * -.1760064569272043 * t101 + t129 * 
			.0136299879940066 * t105 - t134 * 2.73586176e-4 * 
			t109 + t140 * 2.515313720859277e-6 * t113 - t147 * 
			9.466624338548721e-9 * t238;
		vrhoa[i__] = t2 * -.9847450218426965 * t41 - t3 * 
			.3692793831910112 * t152 - t155 * .03109 * t78 + t159 
			* .001664368495390566 * t79 + t164 * .5 * t183 - t48 *
			 .015545 * t199 + t227 * t116 + t98 * t241;
		t245 = sigma * t18;
		t248 = t14 * t26;
		t251 = t24 * t35;
		t256 = 1 / t2 / t25 / t16;
		t257 = t32 * t256;
		t260 = t7 * -.01897807803274269 * t11 + t245 * 
			.001758466769654075 * t21 - t248 * 2.301853696e-5 * 
			t29 + t251 * 1.105943786014421e-7 * t38 - t257 * 
			1.93499927909934e-10 * t149;
		t273 = t7 * -.3310670934776461 * t63 + t245 * 
			.4541896448323908 * t67 - t248 * .875935232 * t71 + 
			t251 * .6468936129513876 * t75 - t257 * 
			.1282578179897027 * t195;
		t274 = t60 * t273;
		t287 = t7 * .1320048426954032 * t101 - t245 * 
			.01022249099550495 * t105 + t248 * 2.05189632e-4 * 
			t109 - t251 * 1.886485290644458e-6 * t113 + t257 * 
			7.099968253911541e-9 * t238;
		vsigmaaa[i__] = t3 * -.7385587663820224 * t260 - t48 * .03109 
			* t274 + t98 * 2. * t287;
		t294 = sigma / t5 / t15;
		t297 = t15 * t120;
		t300 = t14 / t2 / t297;
		t301 = t300 * t105;
		t304 = t24 / t33;
		t305 = t304 * t109;
		t310 = t32 / t5 / t25 / t15;
		t311 = t310 * t113;
		t316 = t143 / t2 / t25 / t297;
		t317 = t316 * t238;
/* Computing 2nd power */
		d__1 = t25;
		t320 = d__1 * d__1;
		t323 = t32 * t14 / t320 / t4;
		t325 = 1 / t112 / t104;
		t327 = t323 * 1.202186354688e-9 * t325;
		t334 = t47 * t163;
		t348 = 1 / t37 / t20;
		t361 = 1 / t15;
		t362 = 1 / t167 / t44 * t361;
		t364 = 1 / t120;
		t365 = t168 * t364;
		t368 = 1 / t54 / t44;
		t369 = t368 * t361;
		t371 = t158 * t364;
		t375 = 1 / t52 / t44 * t361;
		t377 = t174 * t364;
		t381 = 1 / t45 / t44 * t361;
		t383 = t177 * t364;
		t385 = t362 * -6.934554859331846 + t365 * 16.64293166239643 - 
			t369 * 4.305845969834537 + t371 * 12.91753790950361 - 
			t375 * 2.326004811900819 + t377 * 9.304019247603276 - 
			t381 * .3394740105503081 + t383 * 2.036844063301848;
		t390 = t364 * t368;
		t397 = 1 / t162 / t56;
/* Computing 2nd power */
		d__1 = t180;
		t399 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t162;
		t416 = d__1 * d__1;
		t417 = 1 / t416;
/* Computing 2nd power */
		d__1 = t59;
		t419 = d__1 * d__1;
		t420 = 1 / t419;
		t424 = t216 * 2.;
		t425 = t204 * .005495599555936838;
		t427 = t369 * .001831866518645613 * t93;
		t429 = t371 * .005495599555936838 * t93;
		t440 = log(29.60857464321668 / (t49 * 8.157414703487641 + t45 
			* 2.247591863577616 + t52 * .4300972471276643 + t54 * 
			.1911512595127338) + 1.);
		t442 = (t45 * .06901399211255825 + 1.) * t440 * t169;
/* Computing 2nd power */
		d__1 = t213;
		t447 = d__1 * d__1;
		t450 = t83 * 2. / t206 / t89 * t447 * t214;
		t462 = t208 * 1. * (t362 * -.8309097827459833 + t365 * 
			1.99418347859036 - t369 * .4945709824779306 + t371 * 
			1.483712947433792 - t375 * .2001071587498409 + t377 * 
			.8004286349993634 - t381 * .04215565168327908 + t383 *
			 .2529339100996745) * t214;
/* Computing 2nd power */
		d__1 = t206;
		t463 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t92;
		t466 = d__1 * d__1;
		t470 = t83 * 16.0818243221511 / t463 * t447 / t466;
		t474 = t172 * .08837926660346786 * t207 * t213 * t214;
		t492 = 1 / t74 / t66;
		s1 = t227 * 4. * t241 + t98 * (t294 * 1.290714017466165 * 
			t101 - t301 * .140770165298137 + t305 * 
			.00415778512896 - t311 * 5.941411093198738e-5 + t317 *
			 4.259391834712889e-7 - t327) + t159 * 
			.003328736990781132 * t199 - t155 * .06218 * t198 + 
			t334 * 2. * t183 - t3 * .3692793831910112 * (t294 * 
			-.1855634296534841 * t11 + t300 * .030555462130222 * 
			t21 - t304 * 7.112440763733333e-4 * t29 + t310 * 
			6.558764115926639e-6 * t38 - t316 * 
			2.737055459168051e-8 * t149 + t323 * 
			4.368526067484444e-11 * t348) - .6564966812284644 / 
			t5 * t41 - t2 * 1.969490043685393 * t152 + t164 * .5 *
			 t385 * t181 * t78;
		s2 = s1 + t390 * .002219157993854088 * t79 - t159 * 
			.1070677706909338 * t163 * t183 - t48 * 1. * t397 * 
			t399 * t181 * t78 + (t390 * -.002219157993854088 * 
			t60 + t159 * .1070677706909338 * t224 - t48 * .5 * 
			t163 * t385 * t181 + t48 * 1. * t397 * t399 * t181 - 
			t48 * 16.08234158893535 * t417 * t399 * t420 + t424 + 
			t425 + rho * (t427 - t429 + t442 * .03377399869956914 
			- t450 + t462 + t470 - t474) - t334 * 2. * t182) * 
			t116;
		v2rhoa2[i__] = s2 - t48 * .015545 * t60 * (t294 * 
			-3.237100469559206 * t63 + t300 * 8.418187782887979 * 
			t67 - t304 * 23.07323630933333 * t71 + t310 * 
			26.05821057352538 * t75 - t316 * 10.7441027773375 * 
			t195 + t323 * 1.447798010311111 * t492) + t48 * 
			16.08234158893535 * t417 * t399 * t420 * t78 + t164 * 
			1. * t182 * t198 + (t425 + t424 + rho * (t427 - t429 
			- t442 * .03377399869956914 - t450 + t462 + t470 - 
			t474)) * t116 + t98 * (t301 * -.04081692000875529 + 
			t305 * .00215148650496 - t311 * 4.096847697901935e-5 
			+ t317 * 3.56517271655265e-7 - t327);
		t522 = sigma * t128;
		t525 = t14 * t133;
		t528 = t24 * t139;
		t531 = t32 * t146;
		t536 = t143 / t320 / rho;
		t569 = t522 * t105;
		t571 = t525 * t109;
		t573 = t528 * t113;
		t575 = t531 * t238;
		t578 = t536 * 9.01639766016e-10 * t325;
		s1 = t2 * -.9847450218426965 * t260 - t3 * .3692793831910112 *
			 (t122 * .1012164161746277 * t11 - t522 * 
			.01939966305835835 * t21 + t525 * 4.8739598336e-4 * 
			t29 - t528 * 4.697884329742095e-6 * t38 + t531 * 
			2.014091608794051e-8 * t149 - t536 * 
			3.276394550613333e-11 * t348) - t155 * .03109 * t273 
			+ t159 * .001664368495390566 * t274;
		v2rhoasigmaaa[i__] = s1 + t164 * .5 * t182 * t273 - t48 * 
			.015545 * t60 * (t122 * 1.765691165214113 * t63 - 
			t522 * 5.405261547501203 * t67 + t525 * 15.553056768 *
			 t71 - t528 * 18.24987070424126 * t75 + t531 * 
			7.801561447023719 * t195 - t536 * 1.085848507733333 * 
			t492) + t227 * 2. * t287 + t98 * (t122 * 
			-.7040258277088172 * t101 + t569 * .08513264198259285 
			- t571 * .00270795958272 + t573 * 
			4.078761261770162e-5 - t575 * 3.052544510956436e-7 + 
			t578) + t98 * (t569 * .03061269000656647 - t571 * 
			.00161361487872 + t573 * 3.072635773426451e-5 - t575 *
			 2.673879537414487e-7 + t578);
		t589 = sigma * t26;
		t592 = t14 * t35;
		t595 = t24 * t256;
		t599 = t32 / t320;
		v2sigmaaa2[i__] = t3 * -.7385587663820224 * (t18 * 
			.007515880215152465 * t21 - t589 * 2.7347283968e-4 * 
			t29 + t592 * 3.081035732900803e-6 * t38 - t595 * 
			1.433168735431565e-8 * t149 + t599 * 
			2.45729591296e-11 * t348) - t48 * .03109 * t60 * (t18 
			* 2.237187581296339 * t67 - t589 * 8.161051648 * t71 
			+ t592 * 11.09982857637539 * t75 - t595 * 
			5.338139813308978 * t195 + t599 * .8143863808 * t492) 
			+ t98 * 4. * (t18 * -.02295951750492485 * t105 + t589 
			* .00121021115904 * t109 - t592 * 
			2.304476830069838e-5 * t113 + t595 * 
			2.005409653060866e-7 * t238 - t599 * 
			6.76229824512e-10 * t325);
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
} /* rks_xc_hcth120__ */


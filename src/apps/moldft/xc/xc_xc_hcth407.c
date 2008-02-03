/* xc_xc_hcth407.f -- translated by f2c (version 20050501).
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

/* :XC_HCTH407subrstart */
/*    Generated: Sat May  8 18:40:30 GMT 2004 */
/* Subroutine */ int uks_xc_hcth407__(integer *ideriv, integer *npt, 
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
	     t745, t751, t755, t758, t760, t768, t774, t776, t781, t783, t784,
	     t792, t813, t821, t822, t824, t825, t841, t848, t849, t851, t852,
	     t853, t856, t877, t878, t879, t880, t883, t885, t886, t887, t888,
	     t889, t892, t895, t898, t901, t903, t904, t905, t908, t910, t911,
	     t914, t915, t916, t922, t928, t932, t933, t935, t937, t940, t942,
	     t947, t950, t952, t958, t960, t964, t966, t971, t972, t977, t983,
	     t987, t989, t990, t1001, t1002, t1003, t1006, t1009, t1010, 
	    t1022, t1038, t1041, t1052, t1055, t1056, t1061, t1063, t1064, 
	    t1065, t1068, t1069, t1075, t1076, t1078, t1079, t1084, t1085, 
	    t1088, t1097, t1100, t1103, t1107, t1113, t1119, t1123, t1126, 
	    t1128, t1139, t1140, t1142, t1143, t1144, t1154, t1157, t1160, 
	    t1166, t1171, t1173, t1184, t1186, t1189, t1194, t1197, sigmaaa, 
	    sigmaab, sigmabb, t1201, t1203, t1206, t1209, t1211, t1212, t1218,
	     t1225, t1246, t1255, t1272, t1323, t1352, t1362, t1365, t1375, 
	    t1393, t1396, t1399, t1402, t1407, t1463, t1475, t1478, t1493, 
	    t1496, t1499, t1502, t1507, t1561, t1564, t1567, t1571, t1605, 
	    t1617, t1620, t1623, t1627;


/*     A.D. Boese, N.C. Handy */
/*     A new parametriztion of exchange-correlation generalized gradient */
/*     approximation functionals */
/*     J. Chem. Phys. 114 (2001) 5497-5503. */


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
		    zk[i__] = t2 * -.9305257363491 * rhob * (1.08184 - t8 * 
			    .002073356 / t10 + t19 * 5.480992e-5 / t20 - t27 *
			     1.6825664e-7 / t20 / t10 + t36 * 5.858688e-10 / 
			    t37) - rhob * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (1.18777 - t8 * .480584 / t62 + t19 * 
			    .2246964 / t66 - t27 * .07343384 / t66 / t62 + 
			    t36 * .009996768 / t74);
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
		    zk[i__] = t2 * -.9305257363491 * rhoa * (1.08184 - t8 * 
			    .002073356 / t10 + t19 * 5.480992e-5 / t20 - t27 *
			     1.6825664e-7 / t20 / t10 + t36 * 5.858688e-10 / 
			    t37) - rhoa * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (1.18777 - t8 * .480584 / t62 + t19 * 
			    .2246964 / t66 - t27 * .07343384 / t66 / t62 + 
			    t36 * .009996768 / t74);
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
		    s1 = t4 * -.9305257363491 * rhoa * (1.08184 - t10 * 
			    .002073356 / t12 + t21 * 5.480992e-5 / t22 - t29 *
			     1.6825664e-7 / t22 / t12 + t38 * 5.858688e-10 / 
			    t39) - t50 * .03109 * t62 * (1.18777 - t10 * 
			    .480584 / t64 + t21 * .2246964 / t68 - t29 * 
			    .07343384 / t68 / t64 + t38 * .009996768 / t76);
		    s2 = s1 - t84 * .9305257363491 * rhob * (1.08184 - t90 * 
			    .002073356 / t92 + t101 * 5.480992e-5 / t102 - 
			    t109 * 1.6825664e-7 / t102 / t92 + t118 * 
			    5.858688e-10 / t119);
		    zk[i__] = s2 - t130 * .03109 * t142 * (1.18777 - t90 * 
			    .480584 / t144 + t101 * .2246964 / t148 - t109 * 
			    .07343384 / t148 / t144 + t118 * .009996768 / 
			    t156) + (t164 * (-t182 + (t166 * 
			    .06901399211255825 + 1.) * .03799574853701528 * 
			    t193 * t205 * (1. - t211 * 1.) + ((t166 * 
			    .1274696188700087 + 1.) * -.03109 * t227 + t182) *
			     1.923661050931536 * t205 * t211) + t50 * .03109 *
			     t62 + t130 * .03109 * t142) * (t243 * .02654244 /
			     t246 + .589076 - t250 * 6.919848e-4 / t251 + 
			    t250 * 9.1955736e-6 * t243 / t251 / t246 - t260 * 
			    5.44387392e-8 / t261);
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
		    t41 = 1.08184 - t8 * .002073356 * t11 + t19 * 5.480992e-5 
			    * t21 - t27 * 1.6825664e-7 * t29 + t36 * 
			    5.858688e-10 * t38;
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
		    t78 = 1.18777 - t8 * .480584 * t63 + t19 * .2246964 * t67 
			    - t27 * .07343384 * t71 + t36 * .009996768 * t75;
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
			    .9305257363491 * (t87 * .005528949333333333 * t11 
			    - t93 * 3.144353706666667e-4 * t21 + t98 * 
			    2.515331413333333e-6 * t29 - t104 * 
			    1.163347968e-8 * t38 + t111 * 2.49970688e-11 * 
			    t113) - t47 * .03109 * t60 * t78 + t44 * 
			    .001321010150222857 * t122 * t79 + t48 * 1. / 
			    t126 * (-1.853395810515781 / t130 / t49 * t133 - 
			    t122 * 1.28158207914907 * t133 - 
			    .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * 1.281557333333333 * t63 - 
			    t93 * 1.454692266666667 * t67 + t98 * .82714688 * 
			    t71 - t104 * .224126336 * t75 + t111 * 
			    .0213264384 * t159);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t168 = sigmabb * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.002073356 *
			     t11 + t168 * 1.17913264e-4 * t21 - t171 * 
			    9.4324928e-7 * t29 + t174 * 4.36255488e-9 * t38 - 
			    t180 * 9.3739008e-12 * t113) - t48 * .03109 * t60 
			    * (t7 * -.480584 * t63 + t168 * .5455096 * t67 - 
			    t171 * .31018008 * t71 + t174 * .084047376 * t75 
			    - t180 * .0079974144 * t159);
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
		    t41 = 1.08184 - t8 * .002073356 * t11 + t19 * 5.480992e-5 
			    * t21 - t27 * 1.6825664e-7 * t29 + t36 * 
			    5.858688e-10 * t38;
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
		    t78 = 1.18777 - t8 * .480584 * t63 + t19 * .2246964 * t67 
			    - t27 * .07343384 * t71 + t36 * .009996768 * t75;
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
			    .9305257363491 * (t87 * .005528949333333333 * t11 
			    - t93 * 3.144353706666667e-4 * t21 + t98 * 
			    2.515331413333333e-6 * t29 - t104 * 
			    1.163347968e-8 * t38 + t111 * 2.49970688e-11 * 
			    t113) - t47 * .03109 * t60 * t78 + t44 * 
			    .001321010150222857 * t122 * t79 + t48 * 1. / 
			    t126 * (-1.853395810515781 / t130 / t49 * t133 - 
			    t122 * 1.28158207914907 * t133 - 
			    .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * 1.281557333333333 * t63 - 
			    t93 * 1.454692266666667 * t67 + t98 * .82714688 * 
			    t71 - t104 * .224126336 * t75 + t111 * 
			    .0213264384 * t159);
		    vrhob[i__] = 0.;
		    t168 = sigmaaa * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.002073356 *
			     t11 + t168 * 1.17913264e-4 * t21 - t171 * 
			    9.4324928e-7 * t29 + t174 * 4.36255488e-9 * t38 - 
			    t180 * 9.3739008e-12 * t113) - t48 * .03109 * t60 
			    * (t7 * -.480584 * t63 + t168 * .5455096 * t67 - 
			    t171 * .31018008 * t71 + t174 * .084047376 * t75 
			    - t180 * .0079974144 * t159);
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
		    t43 = 1.08184 - t10 * .002073356 * t13 + t21 * 
			    5.480992e-5 * t23 - t29 * 1.6825664e-7 * t31 + 
			    t38 * 5.858688e-10 * t40;
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
		    t80 = 1.18777 - t10 * .480584 * t65 + t21 * .2246964 * 
			    t69 - t29 * .07343384 * t73 + t38 * .009996768 * 
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
		    t123 = 1.08184 - t90 * .002073356 * t93 + t101 * 
			    5.480992e-5 * t103 - t109 * 1.6825664e-7 * t111 + 
			    t118 * 5.858688e-10 * t120;
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
		    t160 = 1.18777 - t90 * .480584 * t145 + t101 * .2246964 * 
			    t149 - t109 * .07343384 * t153 + t118 * 
			    .009996768 * t157;
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
		    t265 = t243 * .02654244 * t247 + .589076 - t250 * 
			    6.919848e-4 * t252 + t255 * 9.1955736e-6 * t257 - 
			    t260 * 5.44387392e-8 * t262;
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
			    t272 * .005528949333333333 * t13 - t278 * 
			    3.144353706666667e-4 * t23 + t283 * 
			    2.515331413333333e-6 * t31 - t289 * 
			    1.163347968e-8 * t40 + t296 * 2.49970688e-11 * 
			    t298) - t304 * .03109 * t80 + t308 * 
			    .001321010150222857 * t81;
		    vrhoa[i__] = s1 + t50 * 1. * t312 * t329 * t330 * t80 - 
			    t50 * .03109 * t62 * (t272 * 1.281557333333333 * 
			    t65 - t278 * 1.454692266666667 * t69 + t283 * 
			    .82714688 * t73 - t289 * .224126336 * t77 + t296 *
			     .0213264384 * t344) + (-t182 + t216 + t233 + 
			    t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t404 * t213 + t194 * 
			    .03799574853701528 * t205 * (-t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t404 * t211 + t442 * 
			    7.694644203726145 - t445) + t304 * .03109 - t308 *
			     .001321010150222857 * t62 - t50 * 1. * t312 * 
			    t329 * t330) * t265 + t240 * (t272 * -.03538992 * 
			    t247 + t459 * .00205763232 * t272 - t462 * 
			    4.78540512e-5 * t272 + t465 * 5.110337088e-7 * 
			    t272 - t470 * 1.7420396544e-9 * t272);
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
			    * (t480 * .005528949333333333 * t93 - t486 * 
			    3.144353706666667e-4 * t103 + t491 * 
			    2.515331413333333e-6 * t111 - t497 * 
			    1.163347968e-8 * t120 + t504 * 2.49970688e-11 * 
			    t506) - t512 * .03109 * t160 + t516 * 
			    .001321010150222857 * t161;
		    vrhob[i__] = s1 + t130 * 1. * t520 * t537 * t538 * t160 - 
			    t130 * .03109 * t142 * (t480 * 1.281557333333333 *
			     t145 - t486 * 1.454692266666667 * t149 + t491 * 
			    .82714688 * t153 - t497 * .224126336 * t157 + 
			    t504 * .0213264384 * t552) + (-t182 + t216 + t233 
			    + t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t565 * t213 + t194 * 
			    .03799574853701528 * t205 * (t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t565 * t211 - t442 * 
			    7.694644203726145 - t445) + t512 * .03109 - t516 *
			     .001321010150222857 * t142 - t130 * 1. * t520 * 
			    t537 * t538) * t265 + t240 * (t480 * -.03538992 * 
			    t247 + t459 * .00205763232 * t480 - t462 * 
			    4.78540512e-5 * t480 + t465 * 5.110337088e-7 * 
			    t480 - t470 * 1.7420396544e-9 * t480);
		    t602 = sigmaaa * t20;
		    t605 = t16 * t28;
		    t608 = t26 * t37;
		    t614 = t34 / t4 / t27 / t18;
		    vsigmaaa[i__] = t5 * -.9305257363491 * (t9 * -.002073356 *
			     t13 + t602 * 1.17913264e-4 * t23 - t605 * 
			    9.4324928e-7 * t31 + t608 * 4.36255488e-9 * t40 - 
			    t614 * 9.3739008e-12 * t298) - t50 * .03109 * t62 
			    * (t9 * -.480584 * t65 + t602 * .5455096 * t69 - 
			    t605 * .31018008 * t73 + t608 * .084047376 * t77 
			    - t614 * .0079974144 * t344) + t240 * (t9 * 
			    .01327122 * t247 - t459 * 7.7161212e-4 * t9 + 
			    t462 * 1.79452692e-5 * t9 - t465 * 1.916376408e-7 
			    * t9 + t470 * 6.532648704e-10 * t9);
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t660 = t114 / t84 / t107 / t98;
		    vsigmabb[i__] = t85 * -.9305257363491 * (t89 * 
			    -.002073356 * t93 + t648 * 1.17913264e-4 * t103 - 
			    t651 * 9.4324928e-7 * t111 + t654 * 4.36255488e-9 
			    * t120 - t660 * 9.3739008e-12 * t506) - t130 * 
			    .03109 * t142 * (t89 * -.480584 * t145 + t648 * 
			    .5455096 * t149 - t651 * .31018008 * t153 + t654 *
			     .084047376 * t157 - t660 * .0079974144 * t552) + 
			    t240 * (t89 * .01327122 * t247 - t459 * 
			    7.7161212e-4 * t89 + t462 * 1.79452692e-5 * t89 - 
			    t465 * 1.916376408e-7 * t89 + t470 * 
			    6.532648704e-10 * t89);
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
		    t41 = 1.08184 - t8 * .002073356 * t11 + t19 * 5.480992e-5 
			    * t21 - t27 * 1.6825664e-7 * t29 + t36 * 
			    5.858688e-10 * t38;
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
		    t78 = 1.18777 - t8 * .480584 * t63 + t19 * .2246964 * t67 
			    - t27 * .07343384 * t71 + t36 * .009996768 * t75;
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
		    t116 = t87 * .005528949333333333 * t11 - t93 * 
			    3.144353706666667e-4 * t21 + t98 * 
			    2.515331413333333e-6 * t29 - t104 * 
			    1.163347968e-8 * t38 + t111 * 2.49970688e-11 * 
			    t113;
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
		    t162 = t87 * 1.281557333333333 * t63 - t93 * 
			    1.454692266666667 * t67 + t98 * .82714688 * t71 - 
			    t104 * .224126336 * t75 + t111 * .0213264384 * 
			    t159;
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
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.002073356 *
			     t11 + t168 * 1.17913264e-4 * t21 - t171 * 
			    9.4324928e-7 * t29 + t174 * 4.36255488e-9 * t38 - 
			    t180 * 9.3739008e-12 * t113) - t48 * .03109 * t60 
			    * (t7 * -.480584 * t63 + t168 * .5455096 * t67 - 
			    t171 * .31018008 * t71 + t174 * .084047376 * t75 
			    - t180 * .0079974144 * t159);
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
			    -.02027281422222222 * t11 + t213 * 
			    .002050399473777778 * t21 - t217 * 
			    2.934593729422222e-5 * t29 + t223 * 
			    2.162145348266667e-7 * t38 - t229 * 
			    8.546531191466667e-10 * t113 + t236 * 
			    1.333177002666667e-12 * t238) + t47 * 2. * t127 * 
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
			    t207 * -4.699043555555556 * t63 + t213 * 
			    9.896548266666667 * t67 - t217 * 
			    8.995993671111111 * t71 + t223 * 
			    3.938242261333333 * t75 - t229 * 
			    .7838151338666667 * t159 + t236 * .0568705024 * 
			    t319);
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
			    1.26206688e-4 * t21 - t328 * 2.829804672e-6 * t29 
			    + t331 * 2.4406656e-8 * t38 - t334 * 
			    1.0729648128e-10 * t113 + t338 * 1.87478016e-13 * 
			    t238) - t48 * .03109 * t60 * (t18 * .6416264 * 
			    t67 - t328 * .838564 * t71 + t331 * .438250176 * 
			    t75 - t334 * .0992275584 * t159 + t338 * 
			    .0079974144 * t319);
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
		    t41 = 1.08184 - t8 * .002073356 * t11 + t19 * 5.480992e-5 
			    * t21 - t27 * 1.6825664e-7 * t29 + t36 * 
			    5.858688e-10 * t38;
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
		    t78 = 1.18777 - t8 * .480584 * t63 + t19 * .2246964 * t67 
			    - t27 * .07343384 * t71 + t36 * .009996768 * t75;
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
		    t116 = t87 * .005528949333333333 * t11 - t93 * 
			    3.144353706666667e-4 * t21 + t98 * 
			    2.515331413333333e-6 * t29 - t104 * 
			    1.163347968e-8 * t38 + t111 * 2.49970688e-11 * 
			    t113;
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
		    t162 = t87 * 1.281557333333333 * t63 - t93 * 
			    1.454692266666667 * t67 + t98 * .82714688 * t71 - 
			    t104 * .224126336 * t75 + t111 * .0213264384 * 
			    t159;
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
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.002073356 *
			     t11 + t168 * 1.17913264e-4 * t21 - t171 * 
			    9.4324928e-7 * t29 + t174 * 4.36255488e-9 * t38 - 
			    t180 * 9.3739008e-12 * t113) - t48 * .03109 * t60 
			    * (t7 * -.480584 * t63 + t168 * .5455096 * t67 - 
			    t171 * .31018008 * t71 + t174 * .084047376 * t75 
			    - t180 * .0079974144 * t159);
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
			    -.02027281422222222 * t11 + t213 * 
			    .002050399473777778 * t21 - t217 * 
			    2.934593729422222e-5 * t29 + t223 * 
			    2.162145348266667e-7 * t38 - t229 * 
			    8.546531191466667e-10 * t113 + t236 * 
			    1.333177002666667e-12 * t238) + t47 * 2. * t127 * 
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
			    t207 * -4.699043555555556 * t63 + t213 * 
			    9.896548266666667 * t67 - t217 * 
			    8.995993671111111 * t71 + t223 * 
			    3.938242261333333 * t75 - t229 * 
			    .7838151338666667 * t159 + t236 * .0568705024 * 
			    t319);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t328 = sigmaaa * t26;
		    t331 = t14 * t35;
		    t334 = t24 * t179;
		    t338 = t32 / t233;
		    v2sigmaaa2[i__] = t3 * -.9305257363491 * (t18 * 
			    1.26206688e-4 * t21 - t328 * 2.829804672e-6 * t29 
			    + t331 * 2.4406656e-8 * t38 - t334 * 
			    1.0729648128e-10 * t113 + t338 * 1.87478016e-13 * 
			    t238) - t48 * .03109 * t60 * (t18 * .6416264 * 
			    t67 - t328 * .838564 * t71 + t331 * .438250176 * 
			    t75 - t334 * .0992275584 * t159 + t338 * 
			    .0079974144 * t319);
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
		    t43 = 1.08184 - t10 * .002073356 * t13 + t21 * 
			    5.480992e-5 * t23 - t29 * 1.6825664e-7 * t31 + 
			    t38 * 5.858688e-10 * t40;
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
		    t80 = 1.18777 - t10 * .480584 * t65 + t21 * .2246964 * 
			    t69 - t29 * .07343384 * t73 + t38 * .009996768 * 
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
		    t123 = 1.08184 - t90 * .002073356 * t93 + t101 * 
			    5.480992e-5 * t103 - t109 * 1.6825664e-7 * t111 + 
			    t118 * 5.858688e-10 * t120;
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
		    t160 = 1.18777 - t90 * .480584 * t145 + t101 * .2246964 * 
			    t149 - t109 * .07343384 * t153 + t118 * 
			    .009996768 * t157;
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
		    t265 = t243 * .02654244 * t247 + .589076 - t250 * 
			    6.919848e-4 * t252 + t255 * 9.1955736e-6 * t257 - 
			    t260 * 5.44387392e-8 * t262;
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
		    t301 = t272 * .005528949333333333 * t13 - t278 * 
			    3.144353706666667e-4 * t23 + t283 * 
			    2.515331413333333e-6 * t31 - t289 * 
			    1.163347968e-8 * t40 + t296 * 2.49970688e-11 * 
			    t298;
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
		    t347 = t272 * 1.281557333333333 * t65 - t278 * 
			    1.454692266666667 * t69 + t283 * .82714688 * t73 
			    - t289 * .224126336 * t77 + t296 * .0213264384 * 
			    t344;
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
		    t473 = t272 * -.03538992 * t247 + t459 * .00205763232 * 
			    t272 - t462 * 4.78540512e-5 * t272 + t465 * 
			    5.110337088e-7 * t272 - t470 * 1.7420396544e-9 * 
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
		    t509 = t480 * .005528949333333333 * t93 - t486 * 
			    3.144353706666667e-4 * t103 + t491 * 
			    2.515331413333333e-6 * t111 - t497 * 
			    1.163347968e-8 * t120 + t504 * 2.49970688e-11 * 
			    t506;
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
		    t555 = t480 * 1.281557333333333 * t145 - t486 * 
			    1.454692266666667 * t149 + t491 * .82714688 * 
			    t153 - t497 * .224126336 * t157 + t504 * 
			    .0213264384 * t552;
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
		    t598 = t480 * -.03538992 * t247 + t459 * .00205763232 * 
			    t480 - t462 * 4.78540512e-5 * t480 + t465 * 
			    5.110337088e-7 * t480 - t470 * 1.7420396544e-9 * 
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
		    t617 = t9 * -.002073356 * t13 + t602 * 1.17913264e-4 * 
			    t23 - t605 * 9.4324928e-7 * t31 + t608 * 
			    4.36255488e-9 * t40 - t614 * 9.3739008e-12 * t298;
		    t630 = t9 * -.480584 * t65 + t602 * .5455096 * t69 - t605 
			    * .31018008 * t73 + t608 * .084047376 * t77 - 
			    t614 * .0079974144 * t344;
		    t631 = t62 * t630;
		    t644 = t9 * .01327122 * t247 - t459 * 7.7161212e-4 * t9 + 
			    t462 * 1.79452692e-5 * t9 - t465 * 1.916376408e-7 
			    * t9 + t470 * 6.532648704e-10 * t9;
		    vsigmaaa[i__] = t5 * -.9305257363491 * t617 - t50 * 
			    .03109 * t631 + t240 * t644;
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t659 = 1 / t84 / t107 / t98;
		    t660 = t114 * t659;
		    t663 = t89 * -.002073356 * t93 + t648 * 1.17913264e-4 * 
			    t103 - t651 * 9.4324928e-7 * t111 + t654 * 
			    4.36255488e-9 * t120 - t660 * 9.3739008e-12 * 
			    t506;
		    t676 = t89 * -.480584 * t145 + t648 * .5455096 * t149 - 
			    t651 * .31018008 * t153 + t654 * .084047376 * 
			    t157 - t660 * .0079974144 * t552;
		    t677 = t142 * t676;
		    t690 = t89 * .01327122 * t247 - t459 * 7.7161212e-4 * t89 
			    + t462 * 1.79452692e-5 * t89 - t465 * 
			    1.916376408e-7 * t89 + t470 * 6.532648704e-10 * 
			    t89;
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
		    t768 = t49 * t312;
		    t774 = 1 / t311 / t58;
/* Computing 2nd power */
		    d__1 = t329;
		    t776 = d__1 * d__1;
		    t781 = 1 / t269;
		    t783 = 1 / t56 / t46;
		    t784 = t781 * t783;
		    t792 = 1 / t17;
		    t813 = -1.544496508763151 / t316 / t46 * t792 + t317 * 
			    3.706791621031562 * t781 - t783 * 
			    .854388052766047 * t792 + t307 * 
			    2.563164158298141 * t781 - .4111834438919023 / 
			    t54 / t46 * t792 + t323 * 1.644733775567609 * 
			    t781 - .05346380647307093 / t47 / t46 * t792 + 
			    t326 * .3207828388384256 * t781;
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
		    t849 = t354 * .005495599555936838;
		    t851 = t379 * .001748158859896213;
		    t852 = t393 * 2.249999913366216;
		    t853 = t437 * 3.847322101863073;
		    t856 = t444 * 15.38928840745229;
		    t877 = 1 / t208 / t164;
		    t878 = t351 * t877;
		    t879 = t878 * t180;
		    t880 = t879 * .005495599555936838;
		    t883 = 1 / t174 / t165 * t210;
		    t885 = t883 * 5.827196199654043e-4 * t378;
		    t886 = t883 * t180;
		    t887 = t886 * .001831866518645613;
		    t888 = t408 * t412;
		    t889 = t231 * t888;
		    t892 = t435 * t404 * t211;
		    t895 = t353 * t377 * t415;
		    t898 = t878 * .001748158859896213 * t378;
		    t901 = t353 * t193 * t404 * t213;
/* Computing 2nd power */
		    d__1 = t199;
		    t903 = d__1 * d__1;
		    t904 = 1 / t903;
/* Computing 2nd power */
		    d__1 = t397;
		    t905 = d__1 * d__1;
		    t908 = t352 * 2.;
		    t910 = t196 * 2. * t877;
		    t911 = -t908 + t910;
/* Computing 2nd power */
		    d__1 = t203;
		    t914 = d__1 * d__1;
		    t915 = 1 / t914;
/* Computing 2nd power */
		    d__1 = t401;
		    t916 = d__1 * d__1;
		    t922 = t904 * .4444444444444444 * t905 + t199 * 
			    1.333333333333333 * t911 + t915 * 
			    .4444444444444444 * t916 - t203 * 
			    1.333333333333333 * t911;
		    t928 = t389 * t390 * t404 * t213;
		    t932 = t353 * t357 * t372 * t373;
		    t933 = t932 * .08837926660346786;
		    t935 = t389 * t391 * t415;
/* Computing 2nd power */
		    d__1 = t381;
		    t937 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t388;
		    t940 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t192;
		    t942 = d__1 * d__1;
		    t947 = t184 * 33.30964519106732 / t937 * t940 / t942 * 
			    t205 * t213;
		    t950 = 1 / t361 / t165 * t210;
		    t952 = t362 * t877;
		    t958 = 1 / t172 / t165 * t210;
		    t960 = t366 * t877;
		    t964 = 1 / t166 / t165 * t210;
		    t966 = t369 * t877;
		    t971 = t383 * 1.124999956683108 * (t950 * 
			    -1.132974264373283 + t952 * 2.71913823449588 - 
			    t883 * .4994648585728036 + t878 * 
			    1.498394575718411 - t958 * .1075243117819161 + 
			    t960 * .4300972471276643 - t964 * 
			    .04247805766949639 + t966 * .2548683460169784) * 
			    t392;
		    t972 = -t880 - t885 + t887 - t889 * 61.55715362980916 + 
			    t892 * 3.847322101863073 - t895 * 
			    .001748158859896213 + t898 - t901 * 
			    .001748158859896213 + t194 * .03799574853701528 * 
			    t922 * t213 - t928 * 2.249999913366216 - t933 - 
			    t935 * 2.249999913366216 - t947 - t971;
		    t977 = t353 * .05176049209143758 * t382 * t388 * t390 * 
			    t214;
		    t983 = t184 * 2.249999913366216 / t381 / t189 * t940 * 
			    t392;
/* Computing 2nd power */
		    d__1 = t372;
		    t987 = d__1 * d__1;
		    t989 = t168 / t356 / t176 * t987 * t373;
		    t990 = t989 * 2.;
		    t1001 = t358 * (t950 * -.8309097827459833 + t952 * 
			    1.99418347859036 - t883 * .4945709824779306 + 
			    t878 * 1.483712947433792 - t958 * 
			    .2001071587498409 + t960 * .8004286349993634 - 
			    t964 * .04215565168327908 + t966 * 
			    .2529339100996745) * t373;
		    t1002 = t1001 * 1.;
/* Computing 2nd power */
		    d__1 = t356;
		    t1003 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t179;
		    t1006 = d__1 * d__1;
		    t1009 = t168 / t1003 * t987 / t1006;
		    t1010 = t1009 * 16.0818243221511;
/* Computing 2nd power */
		    d__1 = t428;
		    t1022 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t421;
		    t1038 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t226;
		    t1041 = d__1 * d__1;
		    t1052 = t883 * 8.806734334819047e-4 * t227 - t878 * 
			    .002642020300445714 * t227 - t353 * 
			    .08497974591333914 * t422 * t428 * t429 - t218 * 
			    2. / t421 / t223 * t1022 * t429 + t423 * 1. * (
			    t950 * -1.544496508763151 + t952 * 
			    3.706791621031562 - t883 * .854388052766047 + 
			    t878 * 2.563164158298141 - t958 * 
			    .4111834438919023 + t960 * 1.644733775567609 - 
			    t964 * .05346380647307093 + t966 * 
			    .3207828388384256) * t429 + t218 * 
			    32.1646831778707 / t1038 * t1022 / t1041 - t886 * 
			    .001831866518645613 + t879 * .005495599555936838 
			    + t932 * .08837926660346786 + t989 * 2. - t1001 * 
			    1. - t1009 * 16.0818243221511;
		    t1055 = t1052 * 1.923661050931536 * t205 * t211;
		    t1056 = t436 * t409;
		    t1061 = t439 * t409;
		    t1063 = t206 * t210;
		    t1064 = t1063 * 12.;
		    t1065 = t888 * 32.;
		    t1068 = t207 / t209 / t208;
		    t1069 = t1068 * 20.;
		    t1075 = t231 * 38.47322101863073 * t1068;
		    t1076 = t439 * t413;
		    t1078 = t231 * t1063;
		    t1079 = t1078 * 23.08393261117844;
		    t1084 = t436 * 15.38928840745229 * t413;
		    t1085 = t977 + t983 - t990 + t1002 + t1010 + t1055 + 
			    t1056 * 15.38928840745229 + t230 * 
			    1.923661050931536 * t922 * t211 + t1061 * 
			    15.38928840745229 + t194 * .03799574853701528 * 
			    t205 * (-t1064 + t1065 - t1069) + t1075 - t1076 * 
			    15.38928840745229 + t1079 + t194 * 
			    .07599149707403056 * t404 * t415 - t1084;
		    t1088 = t848 + t849 + t417 * .07599149707403056 - t851 - 
			    t852 + t853 + t406 * .07599149707403056 + t442 * 
			    15.38928840745229 - t856 + t440 * 
			    3.847322101863073 - t768 * 2. * t331 - t50 * 
			    32.1646831778707 * t822 * t776 * t825 - t50 * 1. *
			     t312 * t813 * t330 - t784 * 8.806734334819047e-4 
			    * t62 + t308 * .08497974591333914 * t452 + t50 * 
			    2. * t774 * t776 * t330 + t164 * (t972 + t1085);
		    s1 = t455 * 2. * t473 + t240 * (t696 * .12976304 * t247 - 
			    t702 * .00302662912 * t252 + t705 * 
			    1.6053292032e-4 * t702 - t459 * .00754465184 * 
			    t696 - t710 * 3.192632064e-6 * t702 + t462 * 
			    1.754648544e-4 * t696 + t715 * 2.56439568384e-8 * 
			    t702 - t465 * 1.8737902656e-6 * t696 - t722 * 
			    6.9681586176e-11 * t702 + t470 * 6.3874787328e-9 *
			     t696) - .4135669939329333 / t7 * t43 - t4 * 
			    2.4814019635976 * t301 - t5 * .9305257363491 * (
			    t696 * -.02027281422222222 * t13 + t702 * 
			    .002050399473777778 * t23 - t739 * 
			    2.934593729422222e-5 * t31 + t745 * 
			    2.162145348266667e-7 * t40 - t751 * 
			    8.546531191466667e-10 * t298 + t758 * 
			    1.333177002666667e-12 * t760) - t304 * .06218 * 
			    t347 + t768 * 2. * t332 + t308 * 
			    .002642020300445714 * t348;
		    v2rhoa2[i__] = s1 - t50 * 2. * t774 * t776 * t330 * t80 + 
			    t784 * 8.806734334819047e-4 * t81 - t308 * 
			    .08497974591333914 * t312 * t332 + t313 * 1. * 
			    t813 * t330 * t80 + t313 * 2. * t331 * t347 + t50 
			    * 32.1646831778707 * t822 * t776 * t825 * t80 - 
			    t50 * .03109 * t62 * (t696 * -4.699043555555556 * 
			    t65 + t702 * 9.896548266666667 * t69 - t739 * 
			    8.995993671111111 * t73 + t745 * 
			    3.938242261333333 * t77 - t751 * 
			    .7838151338666667 * t344 + t758 * .0568705024 * 
			    t841) + t1088 * t265;
		    t1097 = sigmabb / t87 / t97;
		    t1100 = t97 * t477;
		    t1103 = t96 / t84 / t1100;
		    t1107 = t106 / t115;
		    t1113 = t114 / t87 / t107 / t97;
		    t1119 = t500 / t84 / t107 / t1100;
/* Computing 2nd power */
		    d__1 = t107;
		    t1123 = d__1 * d__1;
		    t1126 = t114 * t96 / t1123 / t86;
		    t1128 = 1 / t119 / t102;
/* Computing 2nd power */
		    d__1 = t519;
		    t1139 = d__1 * d__1;
		    t1140 = 1 / t1139;
/* Computing 2nd power */
		    d__1 = t537;
		    t1142 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t141;
		    t1143 = d__1 * d__1;
		    t1144 = 1 / t1143;
/* Computing 2nd power */
		    d__1 = t559;
		    t1154 = d__1 * d__1;
		    t1157 = t908 + t910;
/* Computing 2nd power */
		    d__1 = t562;
		    t1160 = d__1 * d__1;
		    t1166 = t904 * .4444444444444444 * t1154 + t199 * 
			    1.333333333333333 * t1157 + t915 * 
			    .4444444444444444 * t1160 - t203 * 
			    1.333333333333333 * t1157;
		    t1171 = t435 * t565 * t211;
		    t1173 = -t880 - t885 + t887 + t889 * 61.55715362980916 + 
			    t898 - t933 - t947 - t971 + t977 + t983 - t990 + 
			    t1002 + t194 * .03799574853701528 * t1166 * t213 
			    + t1171 * 3.847322101863073;
		    t1184 = t573 * t413;
		    t1186 = t573 * t409;
		    t1189 = t389 * t391 * t569;
		    t1194 = t353 * t193 * t565 * t213;
		    t1197 = t353 * t377 * t569;
		    t1201 = t389 * t390 * t565 * t213;
		    t1203 = t194 * .07599149707403056 * t565 * t569 + t230 * 
			    1.923661050931536 * t1166 * t211 + t194 * 
			    .03799574853701528 * t205 * (-t1064 - t1065 - 
			    t1069) - t1184 * 15.38928840745229 + t1010 + 
			    t1055 - t1186 * 15.38928840745229 - t1189 * 
			    2.249999913366216 - t1056 * 15.38928840745229 - 
			    t1194 * .001748158859896213 - t1197 * 
			    .001748158859896213 - t1201 * 2.249999913366216 + 
			    t1075 + t1079 - t1084;
		    t1206 = t129 * t520;
		    t1209 = 1 / t477;
		    t1211 = 1 / t136 / t126;
		    t1212 = t1209 * t1211;
		    t1218 = 1 / t519 / t138;
		    t1225 = 1 / t97;
		    t1246 = -1.544496508763151 / t524 / t126 * t1225 + t525 * 
			    3.706791621031562 * t1209 - t1211 * 
			    .854388052766047 * t1225 + t515 * 
			    2.563164158298141 * t1209 - .4111834438919023 / 
			    t134 / t126 * t1225 + t531 * 1.644733775567609 * 
			    t1209 - .05346380647307093 / t127 / t126 * t1225 
			    + t534 * .3207828388384256 * t1209;
		    t1255 = t848 + t849 - t851 - t852 + t853 - t442 * 
			    15.38928840745229 - t856 + t567 * 
			    .07599149707403056 + t571 * .07599149707403056 + 
			    t574 * 3.847322101863073 + t164 * (t1173 + t1203) 
			    - t1206 * 2. * t539 - t1212 * 
			    8.806734334819047e-4 * t142 + t516 * 
			    .08497974591333914 * t583 + t130 * 2. * t1218 * 
			    t1142 * t538 - t130 * 1. * t520 * t1246 * t538 - 
			    t130 * 32.1646831778707 * t1140 * t1142 * t1144;
		    t1272 = 1 / t156 / t148;
		    s1 = -.4135669939329333 / t87 * t123 - t84 * 
			    2.4814019635976 * t509 - t85 * .9305257363491 * (
			    t1097 * -.02027281422222222 * t93 + t1103 * 
			    .002050399473777778 * t103 - t1107 * 
			    2.934593729422222e-5 * t111 + t1113 * 
			    2.162145348266667e-7 * t120 - t1119 * 
			    8.546531191466667e-10 * t506 + t1126 * 
			    1.333177002666667e-12 * t1128) - t512 * .06218 * 
			    t555 + t521 * 2. * t539 * t555 + t130 * 
			    32.1646831778707 * t1140 * t1142 * t1144 * t160 + 
			    t1255 * t265 + t1206 * 2. * t540;
		    s2 = s1 + t516 * .002642020300445714 * t556 - t130 * 
			    .03109 * t142 * (t1097 * -4.699043555555556 * 
			    t145 + t1103 * 9.896548266666667 * t149 - t1107 * 
			    8.995993671111111 * t153 + t1113 * 
			    3.938242261333333 * t157 - t1119 * 
			    .7838151338666667 * t552 + t1126 * .0568705024 * 
			    t1272) + t1212 * 8.806734334819047e-4 * t161;
		    v2rhob2[i__] = s2 - t516 * .08497974591333914 * t520 * 
			    t540 - t130 * 2. * t1218 * t1142 * t538 * t160 + 
			    t521 * 1. * t1246 * t538 * t160 + t586 * 2. * 
			    t598 + t240 * (t1097 * .12976304 * t247 - t1103 * 
			    .00302662912 * t252 + t705 * 1.6053292032e-4 * 
			    t1103 - t459 * .00754465184 * t1097 - t710 * 
			    3.192632064e-6 * t1103 + t462 * 1.754648544e-4 * 
			    t1097 + t715 * 2.56439568384e-8 * t1103 - t465 * 
			    1.8737902656e-6 * t1097 - t722 * 6.9681586176e-11 
			    * t1103 + t470 * 6.3874787328e-9 * t1097);
		    t1323 = -t880 - t885 + t887 + t892 * 1.923661050931536 - 
			    t895 * 8.740794299481065e-4 + t898 - t901 * 
			    8.740794299481065e-4 - t928 * 1.124999956683108 - 
			    t933 - t935 * 1.124999956683108 - t947 - t971 + 
			    t977 + t983 - t990 + t1002 + t1171 * 
			    1.923661050931536;
		    t1352 = t904 * .4444444444444444 * t397 * t559 + t199 * 
			    2.666666666666667 * t196 * t877 + t915 * 
			    .4444444444444444 * t401 * t562 - t203 * 
			    2.666666666666667 * t196 * t877;
		    t1362 = t1184 * -7.694644203726145 + t1010 + t1055 + 
			    t1186 * 7.694644203726145 - t1189 * 
			    1.124999956683108 + t194 * .03799574853701528 * 
			    t205 * (t1064 - t1069) - t1061 * 
			    7.694644203726145 - t1194 * 8.740794299481065e-4 
			    - t1197 * 8.740794299481065e-4 - t1201 * 
			    1.124999956683108 + t194 * .03799574853701528 * 
			    t404 * t569 + t1075 - t1076 * 7.694644203726145 - 
			    t1078 * 23.08393261117844 - t1084 + t194 * 
			    .03799574853701528 * t1352 * t213 + t194 * 
			    .03799574853701528 * t565 * t415 + t230 * 
			    1.923661050931536 * t1352 * t211;
		    t1365 = t849 + t848 - t851 - t852 + t568 + t572 + t853 + 
			    t575 - t856 + t407 + t418 + t441 + t164 * (t1323 
			    + t1362);
		    t1375 = t271 * sigmabb * t479;
		    v2rhoab[i__] = t1365 * t265 + t455 * t598 + t586 * t473 + 
			    t240 * (t272 * -.00302662912 * t252 * sigmabb * 
			    t479 + t705 * 1.6053292032e-4 * sigmaaa * t1375 - 
			    t710 * 3.192632064e-6 * sigmaaa * t1375 + t715 * 
			    2.56439568384e-8 * sigmaaa * t1375 - t722 * 
			    6.9681586176e-11 * sigmaaa * t1375);
		    t1393 = sigmaaa * t277;
		    t1396 = t16 * t282;
		    t1399 = t26 * t288;
		    t1402 = t34 * t295;
		    t1407 = t292 / t755 / rhoa;
		    s1 = t4 * -1.2407009817988 * t617 - t5 * .9305257363491 * 
			    (t271 * .005528949333333333 * t13 - t1393 * 
			    6.509865386666667e-4 * t23 + t1396 * 
			    1.006147720533333e-5 * t31 - t1399 * 
			    7.671789568e-8 * t40 + t1402 * 3.1112101888e-10 * 
			    t298 - t1407 * 4.99941376e-13 * t760) - t304 * 
			    .03109 * t630 + t308 * .001321010150222857 * t631;
		    v2rhoasigmaaa[i__] = s1 + t313 * 1. * t331 * t630 - t50 * 
			    .03109 * t62 * (t271 * 1.281557333333333 * t65 - 
			    t1393 * 3.165696 * t69 + t1396 * 
			    3.063317546666667 * t73 - t1399 * 1.392793472 * 
			    t77 + t1402 * .2859332608 * t344 - t1407 * 
			    .0213264384 * t841) + t455 * t644 + t240 * (t271 *
			     -.03538992 * t247 + t1393 * .00113498592 * t252 
			    - t705 * 6.019984512e-5 * t1393 + t459 * 
			    .00205763232 * t271 + t710 * 1.197237024e-6 * 
			    t1393 - t462 * 4.78540512e-5 * t271 - t715 * 
			    9.6164838144e-9 * t1393 + t465 * 5.110337088e-7 * 
			    t271 + t722 * 2.6130594816e-11 * t1393 - t470 * 
			    1.7420396544e-9 * t271);
		    v2rhoasigmaab[i__] = 0.;
		    t1463 = t272 * t89;
		    v2rhoasigmabb[i__] = t455 * t690 + t240 * (t272 * 
			    .00113498592 * t252 * t89 - t705 * 6.019984512e-5 
			    * t1463 + t710 * 1.197237024e-6 * t1463 - t715 * 
			    9.6164838144e-9 * t1463 + t722 * 2.6130594816e-11 
			    * t1463);
		    t1475 = t252 * t9;
		    t1478 = t480 * t9;
		    v2rhobsigmaaa[i__] = t586 * t644 + t240 * (t480 * 
			    .00113498592 * t1475 - t705 * 6.019984512e-5 * 
			    t1478 + t710 * 1.197237024e-6 * t1478 - t715 * 
			    9.6164838144e-9 * t1478 + t722 * 2.6130594816e-11 
			    * t1478);
		    v2rhobsigmaab[i__] = 0.;
		    t1493 = sigmabb * t485;
		    t1496 = t96 * t490;
		    t1499 = t106 * t496;
		    t1502 = t114 * t503;
		    t1507 = t500 / t1123 / rhob;
		    s1 = t84 * -1.2407009817988 * t663 - t85 * .9305257363491 
			    * (t479 * .005528949333333333 * t93 - t1493 * 
			    6.509865386666667e-4 * t103 + t1496 * 
			    1.006147720533333e-5 * t111 - t1499 * 
			    7.671789568e-8 * t120 + t1502 * 3.1112101888e-10 *
			     t506 - t1507 * 4.99941376e-13 * t1128) - t512 * 
			    .03109 * t676 + t516 * .001321010150222857 * t677;
		    v2rhobsigmabb[i__] = s1 + t521 * 1. * t539 * t676 - t130 *
			     .03109 * t142 * (t479 * 1.281557333333333 * t145 
			    - t1493 * 3.165696 * t149 + t1496 * 
			    3.063317546666667 * t153 - t1499 * 1.392793472 * 
			    t157 + t1502 * .2859332608 * t552 - t1507 * 
			    .0213264384 * t1272) + t586 * t690 + t240 * (t479 
			    * -.03538992 * t247 + t1493 * .00113498592 * t252 
			    - t705 * 6.019984512e-5 * t1493 + t459 * 
			    .00205763232 * t479 + t710 * 1.197237024e-6 * 
			    t1493 - t462 * 4.78540512e-5 * t479 - t715 * 
			    9.6164838144e-9 * t1493 + t465 * 5.110337088e-7 * 
			    t479 + t722 * 2.6130594816e-11 * t1493 - t470 * 
			    1.7420396544e-9 * t479);
		    t1561 = sigmaaa * t28;
		    t1564 = t16 * t37;
		    t1567 = t26 * t613;
		    t1571 = t34 / t755;
		    v2sigmaaa2[i__] = t5 * -.9305257363491 * (t20 * 
			    1.26206688e-4 * t23 - t1561 * 2.829804672e-6 * 
			    t31 + t1564 * 2.4406656e-8 * t40 - t1567 * 
			    1.0729648128e-10 * t298 + t1571 * 1.87478016e-13 *
			     t760) - t50 * .03109 * t62 * (t20 * .6416264 * 
			    t69 - t1561 * .838564 * t73 + t1564 * .438250176 *
			     t77 - t1567 * .0992275584 * t344 + t1571 * 
			    .0079974144 * t841) + t240 * (t20 * -4.2561972e-4 
			    * t252 + t705 * 2.257494192e-5 * t20 - t710 * 
			    4.48963884e-7 * t20 + t715 * 3.6061814304e-9 * 
			    t20 - t722 * 9.798973056e-12 * t20);
		    v2sigmaaaab[i__] = 0.;
		    t1605 = t9 * t89;
		    v2sigmaaabb[i__] = t240 * (t1475 * -4.2561972e-4 * t89 + 
			    t705 * 2.257494192e-5 * t1605 - t710 * 
			    4.48963884e-7 * t1605 + t715 * 3.6061814304e-9 * 
			    t1605 - t722 * 9.798973056e-12 * t1605);
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t1617 = sigmabb * t108;
		    t1620 = t96 * t117;
		    t1623 = t106 * t659;
		    t1627 = t114 / t1123;
		    v2sigmabb2[i__] = t85 * -.9305257363491 * (t100 * 
			    1.26206688e-4 * t103 - t1617 * 2.829804672e-6 * 
			    t111 + t1620 * 2.4406656e-8 * t120 - t1623 * 
			    1.0729648128e-10 * t506 + t1627 * 1.87478016e-13 *
			     t1128) - t130 * .03109 * t142 * (t100 * .6416264 
			    * t149 - t1617 * .838564 * t153 + t1620 * 
			    .438250176 * t157 - t1623 * .0992275584 * t552 + 
			    t1627 * .0079974144 * t1272) + t240 * (t100 * 
			    -4.2561972e-4 * t252 + t705 * 2.257494192e-5 * 
			    t100 - t710 * 4.48963884e-7 * t100 + t715 * 
			    3.6061814304e-9 * t100 - t722 * 9.798973056e-12 * 
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
} /* uks_xc_hcth407__ */

/* Subroutine */ int rks_xc_hcth407__(integer *ideriv, integer *npt, 
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
	    t207, t208, t213, t214, t216, t224, t227, t230, t232, t234, t235, 
	    t240, t241, t256, rho, t260, t273, t274, t287, t294, t297, t300, 
	    t301, t304, t305, t310, t311, t316, t317, t320, t323, t325, t327, 
	    t330, t346, t357, t358, t360, t361, t362, t381, t390, t391, t393, 
	    t394, t397, t398, t400, t404, t406, t410, t412, t414, t419, t426, 
	    t452, t453, t457, t460, t472, t473, t476, t480, t482, t484, t495, 
	    t497, t502, t525, t528, t531, t534, t539, t572, t574, t576, t578, 
	    t581, t592, t595, t598, t602, sigma;


/*     A.D. Boese, N.C. Handy */
/*     A new parametriztion of exchange-correlation generalized gradient */
/*     approximation functionals */
/*     J. Chem. Phys. 114 (2001) 5497-5503. */


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
		zk[i__] = t2 * -.7385587663820224 * rho * (1.08184 - t8 * 
			.003291247495504578 / t10 + t19 * 1.38112343902108e-4 
			/ t20 - t27 * 6.7302656e-7 / t20 / t10 + t36 * 
			3.720034997741387e-9 / t37) - t48 * .03109 * t60 * (
			1.18777 - t8 * .7628795471590852 / t62 + t19 * 
			.5661994483911968 / t66 - t27 * .29373536 / t66 / t62 
			+ t36 * .06347552015792813 / t74) + (rho * -.062182 * 
			(t45 * .1325688999052018 + 1.) * t93 + t48 * .03109 * 
			t60) * (t8 * .04213349717780282 / t100 + .589076 - 
			t19 * .001743692431454588 / t104 + t27 * 
			3.67822944e-5 / t104 / t100 - t36 * 
			3.456644474956098e-7 / t112);
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
		t41 = 1.08184 - t8 * .003291247495504578 * t11 + t19 * 
			1.38112343902108e-4 * t21 - t27 * 6.7302656e-7 * t29 
			+ t36 * 3.720034997741387e-9 * t38;
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
		t78 = 1.18777 - t8 * .7628795471590852 * t63 + t19 * 
			.5661994483911968 * t67 - t27 * .29373536 * t71 + t36 
			* .06347552015792813 * t75;
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
		t116 = t8 * .04213349717780282 * t101 + .589076 - t19 * 
			.001743692431454588 * t105 + t27 * 3.67822944e-5 * 
			t109 - t36 * 3.456644474956098e-7 * t113;
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
			(t123 * .01755331997602442 * t11 - t129 * 
			.001584654969337721 * t21 + t134 * 
			2.012265130666667e-5 * t29 - t140 * 
			1.477359830566614e-7 * t38 + t147 * 
			5.03909330668646e-10 * t149) - t155 * .03109 * t78 + 
			t159 * .001664368495390566 * t79;
		vrhoa[i__] = s1 + t48 * .5 * t163 * t180 * t181 * t78 - t48 * 
			.015545 * t60 * (t123 * 4.068690918181788 * t63 - 
			t129 * 7.331189631570478 * t67 + t134 * 6.61717504 * 
			t71 - t140 * 2.846227052321425 * t75 + t147 * 
			.4299140585511414 * t195) + (t83 * -.062182 * t93 + 
			rho * (t172 * .002747799777968419 * t93 + t83 * 1. / 
			t206 * (t170 * -.99709173929518 - t172 * 
			.7418564737168958 - t175 * .4002143174996817 - t178 * 
			.1264669550498372) / t92) + t155 * .03109 - t159 * 
			.001664368495390566 * t60 - t48 * .5 * t163 * t180 * 
			t181) * t116 + t98 * (t123 * -.1123559924741408 * 
			t101 + t129 * .01036981709164809 * t105 - t134 * 
			3.828324096e-4 * t109 + t140 * 6.489723575522644e-6 * 
			t113 - t147 * 3.51173188852824e-8 * t238);
		t245 = sigma * t18;
		t248 = t14 * t26;
		t251 = t24 * t35;
		t257 = t32 / t2 / t25 / t16;
		vsigmaaa[i__] = t3 * -.7385587663820224 * (t7 * 
			-.01316498998201831 * t11 + t245 * 
			.001188491227003291 * t21 - t248 * 1.509198848e-5 * 
			t29 + t251 * 1.10801987292496e-7 * t38 - t257 * 
			3.779319980014845e-10 * t149) - t48 * .03109 * t60 * (
			t7 * -3.051518188636341 * t63 + t245 * 
			5.498392223677858 * t67 - t248 * 4.96288128 * t71 + 
			t251 * 2.134670289241069 * t75 - t257 * 
			.3224355439133561 * t195) + t98 * 2. * (t7 * 
			.08426699435560563 * t101 - t245 * 
			.007777362818736071 * t105 + t248 * 2.871243072e-4 * 
			t109 - t251 * 4.867292681641983e-6 * t113 + t257 * 
			2.63379891639618e-8 * t238);
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
		t41 = 1.08184 - t8 * .003291247495504578 * t11 + t19 * 
			1.38112343902108e-4 * t21 - t27 * 6.7302656e-7 * t29 
			+ t36 * 3.720034997741387e-9 * t38;
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
		t78 = 1.18777 - t8 * .7628795471590852 * t63 + t19 * 
			.5661994483911968 * t67 - t27 * .29373536 * t71 + t36 
			* .06347552015792813 * t75;
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
		t116 = t8 * .04213349717780282 * t101 + .589076 - t19 * 
			.001743692431454588 * t105 + t27 * 3.67822944e-5 * 
			t109 - t36 * 3.456644474956098e-7 * t113;
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
		t152 = t123 * .01755331997602442 * t11 - t129 * 
			.001584654969337721 * t21 + t134 * 
			2.012265130666667e-5 * t29 - t140 * 
			1.477359830566614e-7 * t38 + t147 * 
			5.03909330668646e-10 * t149;
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
		t198 = t123 * 4.068690918181788 * t63 - t129 * 
			7.331189631570478 * t67 + t134 * 6.61717504 * t71 - 
			t140 * 2.846227052321425 * t75 + t147 * 
			.4299140585511414 * t195;
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
		t230 = t123 * .1123559924741408 * t101;
		t232 = t129 * .01036981709164809 * t105;
		t234 = t134 * 3.828324096e-4 * t109;
		t235 = t140 * t113;
		t238 = 1 / t112 / t100;
		t240 = t147 * 3.51173188852824e-8 * t238;
		t241 = -t230 + t232 - t234 + t235 * 6.489723575522644e-6 - 
			t240;
		vrhoa[i__] = t2 * -.9847450218426965 * t41 - t3 * 
			.3692793831910112 * t152 - t155 * .03109 * t78 + t159 
			* .001664368495390566 * t79 + t164 * .5 * t183 - t48 *
			 .015545 * t199 + t227 * t116 + t98 * t241;
		t245 = sigma * t18;
		t248 = t14 * t26;
		t251 = t24 * t35;
		t256 = 1 / t2 / t25 / t16;
		t257 = t32 * t256;
		t260 = t7 * -.01316498998201831 * t11 + t245 * 
			.001188491227003291 * t21 - t248 * 1.509198848e-5 * 
			t29 + t251 * 1.10801987292496e-7 * t38 - t257 * 
			3.779319980014845e-10 * t149;
		t273 = t7 * -3.051518188636341 * t63 + t245 * 
			5.498392223677858 * t67 - t248 * 4.96288128 * t71 + 
			t251 * 2.134670289241069 * t75 - t257 * 
			.3224355439133561 * t195;
		t274 = t60 * t273;
		t287 = t7 * .08426699435560563 * t101 - t245 * 
			.007777362818736071 * t105 + t248 * 2.871243072e-4 * 
			t109 - t251 * 4.867292681641983e-6 * t113 + t257 * 
			2.63379891639618e-8 * t238;
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
		t327 = t323 * 4.459621515264e-9 * t325;
		t330 = t47 * t163;
		t346 = 1 / t37 / t20;
/* Computing 2nd power */
		d__1 = t162;
		t357 = d__1 * d__1;
		t358 = 1 / t357;
/* Computing 2nd power */
		d__1 = t180;
		t360 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t59;
		t361 = d__1 * d__1;
		t362 = 1 / t361;
		t381 = 1 / t74 / t66;
		t390 = 1 / t15;
		t391 = 1 / t167 / t44 * t390;
		t393 = 1 / t120;
		t394 = t168 * t393;
		t397 = 1 / t54 / t44;
		t398 = t397 * t390;
		t400 = t158 * t393;
		t404 = 1 / t52 / t44 * t390;
		t406 = t174 * t393;
		t410 = 1 / t45 / t44 * t390;
		t412 = t177 * t393;
		t414 = t391 * -6.934554859331846 + t394 * 16.64293166239643 - 
			t398 * 4.305845969834537 + t400 * 12.91753790950361 - 
			t404 * 2.326004811900819 + t406 * 9.304019247603276 - 
			t410 * .3394740105503081 + t412 * 2.036844063301848;
		t419 = t393 * t397;
		t426 = 1 / t162 / t56;
		t452 = t204 * .005495599555936838;
		t453 = t216 * 2.;
/* Computing 2nd power */
		d__1 = t213;
		t457 = d__1 * d__1;
		t460 = t83 * 2. / t206 / t89 * t457 * t214;
		t472 = t208 * 1. * (t391 * -.8309097827459833 + t394 * 
			1.99418347859036 - t398 * .4945709824779306 + t400 * 
			1.483712947433792 - t404 * .2001071587498409 + t406 * 
			.8004286349993634 - t410 * .04215565168327908 + t412 *
			 .2529339100996745) * t214;
/* Computing 2nd power */
		d__1 = t206;
		t473 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t92;
		t476 = d__1 * d__1;
		t480 = t83 * 16.0818243221511 / t473 * t457 / t476;
		t482 = t398 * .001831866518645613 * t93;
		t484 = t400 * .005495599555936838 * t93;
		t495 = log(29.60857464321668 / (t49 * 8.157414703487641 + t45 
			* 2.247591863577616 + t52 * .4300972471276643 + t54 * 
			.1911512595127338) + 1.);
		t497 = (t45 * .06901399211255825 + 1.) * t495 * t169;
		t502 = t172 * .08837926660346786 * t207 * t213 * t214;
		s1 = t227 * 3. * t241 + t98 * (t294 * .8239439448103662 * 
			t101 - t301 * .1065518352468551 + t305 * 
			.00537596439552 - t311 * 1.286791061715555e-4 + t317 *
			 1.291426557904201e-6 - t327) + t330 * 2. * t183 - 
			t155 * .06218 * t198 - t3 * .3692793831910112 * (t294 
			* -.1287243464908457 * t11 + t300 * 
			.02066673166164795 * t21 - t304 * 
			4.695349967075556e-4 * t29 + t310 * 
			5.491506880554649e-6 * t38 - t316 * 
			3.445745456547828e-8 * t149 + t323 * 
			8.532332817066667e-11 * t346) - .6564966812284644 / 
			t5 * t41 - t2 * 1.969490043685393 * t152 + t48 * 
			16.08234158893535 * t358 * t360 * t362 * t78 + t164 * 
			1. * t182 * t198;
		s2 = s1 - t48 * .015545 * t60 * (t294 * -29.83706673333311 * 
			t63 + t300 * 99.75095585979163 * t67 - t304 * 
			143.9358987377778 * t71 + t310 * 100.0251185367385 * 
			t75 - t316 * 31.601445964313 * t195 + t323 * 
			3.6397121536 * t381) + t164 * .5 * t414 * t181 * t78 
			+ t419 * .002219157993854088 * t79 - t159 * 
			.1070677706909338 * t163 * t183;
		v2rhoa2[i__] = s2 - t48 * 1. * t426 * t360 * t181 * t78 + 
			t159 * .003328736990781132 * t199 + (t330 * -2. * 
			t182 + t48 * 1. * t426 * t360 * t181 - t419 * 
			.002219157993854088 * t60 + t159 * .1070677706909338 *
			 t224 - t48 * 16.08234158893535 * t358 * t360 * t362 
			- t48 * .5 * t163 * t414 * t181 + t452 + t453 + rho * 
			(-t460 + t472 + t480 + t482 - t484 + t497 * 
			.03377399869956914 - t502)) * t116 + (t452 + t453 + 
			rho * (-t460 + t472 + t480 + t482 - t484 - t497 * 
			.03377399869956914 - t502)) * t116 + t227 * (-t230 + 
			t232 - t234 + t235 * 6.489723575522644e-6 - t240) + 
			t98 * (t301 * -.03050650990810237 + t305 * 
			.00256852672512 - t311 * 8.108779995105606e-5 + t317 *
			 1.033899552745464e-6 - t327);
		t525 = sigma * t128;
		t528 = t14 * t133;
		t531 = t24 * t139;
		t534 = t32 * t146;
		t539 = t143 / t320 / rho;
		t572 = t525 * t105;
		t574 = t528 * t109;
		t576 = t531 * t113;
		t578 = t534 * t238;
		t581 = t539 * 3.344716136448e-9 * t325;
		s1 = t2 * -.9847450218426965 * t260 - t3 * .3692793831910112 *
			 (t122 * .07021327990409767 * t11 - t525 * 
			.01312306629222938 * t21 + t528 * 
			3.219672705706667e-4 * t29 - t531 * 
			3.897026185830995e-6 * t38 + t534 * 
			2.508722692810574e-8 * t149 - t539 * 6.3992496128e-11 
			* t346) - t155 * .03109 * t273 + t159 * 
			.001664368495390566 * t274;
		v2rhoasigmaaa[i__] = s1 + t164 * .5 * t182 * t273 - t48 * 
			.015545 * t60 * (t122 * 16.27476367272715 * t63 - 
			t525 * 63.81643244748801 * t67 + t528 * 
			98.02616149333333 * t71 - t531 * 70.74949832407171 * 
			t75 + t534 * 23.05621338540804 * t195 - t539 * 
			2.7297841152 * t381) + t227 * 2. * t287 + t98 * (t122 
			* -.4494239698965634 * t101 + t572 * 
			.06435915079766915 - t574 * .00345772468224 + t576 * 
			8.677474426538262e-5 - t578 * 9.158939401002274e-7 + 
			t581) + t98 * (t572 * .02287988243107678 - t574 * 
			.00192639504384 + t576 * 6.081584996329205e-5 - t578 *
			 7.754246645590978e-7 + t581);
		t592 = sigma * t26;
		t595 = t14 * t35;
		t598 = t24 * t256;
		t602 = t32 / t320;
		v2sigmaaa2[i__] = t3 * -.7385587663820224 * (t18 * 
			.00508833481115887 * t21 - t592 * 1.81107499008e-4 * 
			t29 + t595 * 2.479561690203262e-6 * t38 - t598 * 
			1.730369220407337e-8 * t149 + t602 * 4.7994372096e-11 
			* t346) - t48 * .03109 * t60 * (t18 * 
			25.86875544090457 * t67 - t592 * 53.668096 * t71 + 
			t595 * 44.52344258608951 * t75 - t598 * 
			16.0024178634026 * t195 + t602 * 2.0473380864 * t381) 
			+ t98 * 4. * (t18 * -.01715991182330758 * t105 + t592 
			* .00144479628288 * t109 - t595 * 
			4.561188747246904e-5 * t113 + t598 * 
			5.815684984193233e-7 * t238 - t602 * 
			2.508537102336e-9 * t325);
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
} /* rks_xc_hcth407__ */


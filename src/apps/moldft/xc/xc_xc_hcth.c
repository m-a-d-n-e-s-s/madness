/* xc_xc_hcth.f -- translated by f2c (version 20050501).
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

/* :XC_HCTHsubrstart */
/*    Generated: Fri Feb 14 10:46:20 GMT 2003 */
/* Subroutine */ int uks_xc_hcth__(integer *ideriv, integer *npt, doublereal *
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
	    t14, t15, t25, t32, t27, t19, t36, t37, t44, t45, t49, t52, t54, 
	    t60, t62, t66, t74, t16, t17, t22, t29, t34, t38, t39, t46, t47, 
	    t50, t51, t56, t64, t68, t76, t84, t86, t87, t90, t92, t96, t97, 
	    t11, t18, t24, t26, t35, t101, t102, t130, t131, t114, t142, t107,
	     t126, t109, t118, t119, t127, t134, t136, t144, t148, t156, t164,
	     t165, t166, t169, t172, t174, t180, t182, t193, t196, t197, t198,
	     t199, t202, t203, t205, t206, t207, t208, t209, t211, t227, t243,
	     t246, t250, t251, t260, t261, t41, t48, t59, t63, t67, t71, t75, 
	    t78, t79, t93, t98, t104, t111, t113, t122, t129, t133, t159, 
	    t168, t171, t13, t23, rho, t28, t31, t40, t43, t58, t61, t65, t69,
	     t73, t77, t80, t81, t85, t89, t100, t103, t106, t108, t117, t120,
	     t123, t138, t141, t145, t149, t153, t157, t160, t161, t176, t179,
	     t184, t189, t192, t194, t210, t213, t216, t218, t223, t226, t230,
	     t231, t233, t240, t247, t252, t255, t257, t262, t265, t269, t272,
	     t275, t278, t283, t289, t296, t298, t304, t307, t308, t311, t312,
	     t314, t315, t318, t329, t330, t344, t352, t353, t354, t355, t356,
	     t359, t360, t363, t367, t370, t375, t376, t380, t381, t394, t395,
	     t396, t400, t404, t409, t410, t413, t414, rhoa, rhob, t421, t438,
	     t442, t445, t459, t462, t465, t470, t477, t480, t483, t486, t491,
	     t497, t504, t506, t512, t515, t516, t519, t520, t522, t523, t526,
	     t537, t538, t552, t565, t602, t605, t608, t614, t648, t651, t654,
	     t660, t33, t116, t128, t132, t146, t147, t162, t163, t217, t229, 
	    t236, t238, t249, t263, t270, t299, t319, t328, t331, t334, t338, 
	    t115, t214, t271, t277, t282, t288, t292, t295, t301, t313, t316, 
	    t317, t323, t326, t332, t347, t348, t351, t357, t358, t361, t362, 
	    t366, t369, t372, t373, t377, t378, t379, t382, t383, t388, t389, 
	    t390, t391, t392, t393, t397, t401, t406, t407, sigma, t408, t412,
	     t415, t417, t418, t422, t423, t428, t429, t435, t436, t437, t439,
	     t440, t441, t444, t452, t455, t469, t473, t479, t485, t490, t496,
	     t500, t503, t509, t521, t524, t525, t531, t534, t539, t540, t555,
	     t556, t559, t562, t567, t568, t569, t571, t572, t573, t574, t575,
	     t583, t586, t598, t613, t617, t630, t631, t644, t659, t663, t676,
	     t677, t690, t694, t697, t700, t704, t710, t716, t720, t723, t725,
	     t736, t739, t740, t741, t742, t744, t751, t754, t755, t759, t761,
	     t765, t767, t771, t773, t777, t778, t782, t784, t785, t787, t788,
	     t789, t790, t791, t794, t795, t796, t799, t802, t803, t815, t831,
	     t834, t845, t848, t849, t850, t851, t854, t856, t857, t860, t861,
	     t862, t868, t872, t875, t877, t882, t884, t888, t895, t896, t908,
	     t910, t913, t916, t922, t923, t926, t928, t930, t931, t932, t933,
	     t934, t936, t937, t940, t941, t946, t948, t954, t955, t960, t963,
	     t967, t984, t990, t991, t996, t1001, t1004, t1005, t1007, t1008, 
	    t1012, t1031, t1066, t1071, t1076, t1083, t1092, t1095, t1098, 
	    t1102, t1108, t1114, t1118, t1121, t1123, t1134, t1135, t1137, 
	    t1138, t1139, t1146, t1149, t1153, t1170, t1176, t1184, t1190, 
	    t1206, t1209, t1212, t1218, sigmaaa, sigmaab, sigmabb, t1227, 
	    t1230, t1232, t1234, t1240, t1242, t1245, t1248, t1252, t1273, 
	    t1286, t1325, t1355, t1362, t1365, t1378, t1397, t1400, t1403, 
	    t1406, t1411, t1467, t1479, t1482, t1497, t1500, t1503, t1506, 
	    t1511, t1565, t1568, t1571, t1575, t1609, t1621, t1624, t1627, 
	    t1631;


/*     F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N.C. Handy */
/*     Development and assessment of new exchange-correlation functionals */
/*     J. Chem. Phys. 109 (1998) 6264-6271 */


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
		    zk[i__] = t2 * -.9305257363491 * rhob * (1.0932 - t8 * 
			    .002976224 / t10 + t19 * 8.95872e-5 / t20 - t27 * 
			    4.3427136e-7 / t20 / t10 + t36 * 1.15035392e-9 / 
			    t37) - rhob * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (.222601 - t8 * .00677244 / t62 - t19 
			    * 5.0068e-4 / t66 - t27 * .006419968 / t66 / t62 
			    + t36 * .002486336 / t74);
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
		    zk[i__] = t2 * -.9305257363491 * rhoa * (1.0932 - t8 * 
			    .002976224 / t10 + t19 * 8.95872e-5 / t20 - t27 * 
			    4.3427136e-7 / t20 / t10 + t36 * 1.15035392e-9 / 
			    t37) - rhoa * .03109 * (t45 * .1274696188700087 + 
			    1.) * t60 * (.222601 - t8 * .00677244 / t62 - t19 
			    * 5.0068e-4 / t66 - t27 * .006419968 / t66 / t62 
			    + t36 * .002486336 / t74);
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
		    s1 = t4 * -.9305257363491 * rhoa * (1.0932 - t10 * 
			    .002976224 / t12 + t21 * 8.95872e-5 / t22 - t29 * 
			    4.3427136e-7 / t22 / t12 + t38 * 1.15035392e-9 / 
			    t39) - t50 * .03109 * t62 * (.222601 - t10 * 
			    .00677244 / t64 - t21 * 5.0068e-4 / t68 - t29 * 
			    .006419968 / t68 / t64 + t38 * .002486336 / t76);
		    zk[i__] = s1 - t84 * .9305257363491 * rhob * (1.0932 - 
			    t90 * .002976224 / t92 + t101 * 8.95872e-5 / t102 
			    - t109 * 4.3427136e-7 / t102 / t92 + t118 * 
			    1.15035392e-9 / t119) - t130 * .03109 * t142 * (
			    .222601 - t90 * .00677244 / t144 - t101 * 
			    5.0068e-4 / t148 - t109 * .006419968 / t148 / 
			    t144 + t118 * .002486336 / t156) + (t164 * (-t182 
			    + (t166 * .06901399211255825 + 1.) * 
			    .03799574853701528 * t193 * t205 * (1. - t211 * 
			    1.) + ((t166 * .1274696188700087 + 1.) * -.03109 *
			     t227 + t182) * 1.923661050931536 * t205 * t211) 
			    + t50 * .03109 * t62 + t130 * .03109 * t142) * (
			    t243 * .02011722 / t246 + .729974 - t250 * 
			    4.15548e-4 / t251 + t250 * 1.74649824e-6 * t243 / 
			    t251 / t246 - t260 * 5.80422672e-9 / t261);
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
		    t41 = 1.0932 - t8 * .002976224 * t11 + t19 * 8.95872e-5 * 
			    t21 - t27 * 4.3427136e-7 * t29 + t36 * 
			    1.15035392e-9 * t38;
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
		    t78 = .222601 - t8 * .00677244 * t63 - t19 * 5.0068e-4 * 
			    t67 - t27 * .006419968 * t71 + t36 * .002486336 * 
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
			    .9305257363491 * (t87 * .007936597333333333 * t11 
			    - t93 * 5.095447893333333e-4 * t21 + t98 * 
			    5.38536448e-6 * t29 - t104 * 2.616712533333333e-8 
			    * t38 + t111 * 4.908176725333333e-11 * t113) - 
			    t47 * .03109 * t60 * t78 + t44 * 
			    .001321010150222857 * t122 * t79 + t48 * 1. / 
			    t126 * (-1.853395810515781 / t130 / t49 * t133 - 
			    t122 * 1.28158207914907 * t133 - 
			    .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * .01805984 * t63 - t93 * 
			    9.416746666666667e-4 * t67 + t98 * 
			    .05082568533333333 * t71 - t104 * 
			    .03679286613333333 * t75 + t111 * 
			    .005304183466666667 * t159);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t168 = sigmabb * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.002976224 *
			     t11 + t168 * 1.91079296e-4 * t21 - t171 * 
			    2.01951168e-6 * t29 + t174 * 9.812672e-9 * t38 - 
			    t180 * 1.840566272e-11 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.00677244 * t63 + t168 * 3.53128e-4 *
			     t67 - t171 * .019059632 * t71 + t174 * 
			    .0137973248 * t75 - t180 * .0019890688 * t159);
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
		    t41 = 1.0932 - t8 * .002976224 * t11 + t19 * 8.95872e-5 * 
			    t21 - t27 * 4.3427136e-7 * t29 + t36 * 
			    1.15035392e-9 * t38;
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
		    t78 = .222601 - t8 * .00677244 * t63 - t19 * 5.0068e-4 * 
			    t67 - t27 * .006419968 * t71 + t36 * .002486336 * 
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
			    .9305257363491 * (t87 * .007936597333333333 * t11 
			    - t93 * 5.095447893333333e-4 * t21 + t98 * 
			    5.38536448e-6 * t29 - t104 * 2.616712533333333e-8 
			    * t38 + t111 * 4.908176725333333e-11 * t113) - 
			    t47 * .03109 * t60 * t78 + t44 * 
			    .001321010150222857 * t122 * t79 + t48 * 1. / 
			    t126 * (-1.853395810515781 / t130 / t49 * t133 - 
			    t122 * 1.28158207914907 * t133 - 
			    .8223668877838045 / t52 * t133 - 
			    .1603914194192128 / t45 * t133) / t59 * t78 - t48 
			    * .03109 * t60 * (t87 * .01805984 * t63 - t93 * 
			    9.416746666666667e-4 * t67 + t98 * 
			    .05082568533333333 * t71 - t104 * 
			    .03679286613333333 * t75 + t111 * 
			    .005304183466666667 * t159);
		    vrhob[i__] = 0.;
		    t168 = sigmaaa * t18;
		    t171 = t14 * t26;
		    t174 = t24 * t35;
		    t180 = t32 / t2 / t25 / t16;
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.002976224 *
			     t11 + t168 * 1.91079296e-4 * t21 - t171 * 
			    2.01951168e-6 * t29 + t174 * 9.812672e-9 * t38 - 
			    t180 * 1.840566272e-11 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.00677244 * t63 + t168 * 3.53128e-4 *
			     t67 - t171 * .019059632 * t71 + t174 * 
			    .0137973248 * t75 - t180 * .0019890688 * t159);
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
		    t43 = 1.0932 - t10 * .002976224 * t13 + t21 * 8.95872e-5 *
			     t23 - t29 * 4.3427136e-7 * t31 + t38 * 
			    1.15035392e-9 * t40;
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
		    t80 = .222601 - t10 * .00677244 * t65 - t21 * 5.0068e-4 * 
			    t69 - t29 * .006419968 * t73 + t38 * .002486336 * 
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
		    t123 = 1.0932 - t90 * .002976224 * t93 + t101 * 
			    8.95872e-5 * t103 - t109 * 4.3427136e-7 * t111 + 
			    t118 * 1.15035392e-9 * t120;
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
		    t160 = .222601 - t90 * .00677244 * t145 - t101 * 
			    5.0068e-4 * t149 - t109 * .006419968 * t153 + 
			    t118 * .002486336 * t157;
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
		    t265 = t243 * .02011722 * t247 + .729974 - t250 * 
			    4.15548e-4 * t252 + t255 * 1.74649824e-6 * t257 - 
			    t260 * 5.80422672e-9 * t262;
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
			    t272 * .007936597333333333 * t13 - t278 * 
			    5.095447893333333e-4 * t23 + t283 * 5.38536448e-6 
			    * t31 - t289 * 2.616712533333333e-8 * t40 + t296 *
			     4.908176725333333e-11 * t298) - t304 * .03109 * 
			    t80 + t308 * .001321010150222857 * t81;
		    vrhoa[i__] = s1 + t50 * 1. * t312 * t329 * t330 * t80 - 
			    t50 * .03109 * t62 * (t272 * .01805984 * t65 - 
			    t278 * 9.416746666666667e-4 * t69 + t283 * 
			    .05082568533333333 * t73 - t289 * 
			    .03679286613333333 * t77 + t296 * 
			    .005304183466666667 * t344) + (-t182 + t216 + 
			    t233 + t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t404 * t213 + t194 * 
			    .03799574853701528 * t205 * (-t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t404 * t211 + t442 * 
			    7.694644203726145 - t445) + t304 * .03109 - t308 *
			     .001321010150222857 * t62 - t50 * 1. * t312 * 
			    t329 * t330) * t265 + t240 * (t272 * -.02682296 * 
			    t247 + t459 * .00126906576 * t272 - t462 * 
			    1.363476096e-5 * t272 + t465 * 7.28718336e-8 * 
			    t272 - t470 * 1.8573525504e-10 * t272);
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
			    * (t480 * .007936597333333333 * t93 - t486 * 
			    5.095447893333333e-4 * t103 + t491 * 
			    5.38536448e-6 * t111 - t497 * 
			    2.616712533333333e-8 * t120 + t504 * 
			    4.908176725333333e-11 * t506) - t512 * .03109 * 
			    t160 + t516 * .001321010150222857 * t161;
		    vrhob[i__] = s1 + t130 * 1. * t520 * t537 * t538 * t160 - 
			    t130 * .03109 * t142 * (t480 * .01805984 * t145 - 
			    t486 * 9.416746666666667e-4 * t149 + t491 * 
			    .05082568533333333 * t153 - t497 * 
			    .03679286613333333 * t157 + t504 * 
			    .005304183466666667 * t552) + (-t182 + t216 + 
			    t233 + t164 * (t355 + t376 - t380 - t394 + t194 * 
			    .03799574853701528 * t565 * t213 + t194 * 
			    .03799574853701528 * t205 * (t410 + t414) + t438 
			    + t230 * 1.923661050931536 * t565 * t211 - t442 * 
			    7.694644203726145 - t445) + t512 * .03109 - t516 *
			     .001321010150222857 * t142 - t130 * 1. * t520 * 
			    t537 * t538) * t265 + t240 * (t480 * -.02682296 * 
			    t247 + t459 * .00126906576 * t480 - t462 * 
			    1.363476096e-5 * t480 + t465 * 7.28718336e-8 * 
			    t480 - t470 * 1.8573525504e-10 * t480);
		    t602 = sigmaaa * t20;
		    t605 = t16 * t28;
		    t608 = t26 * t37;
		    t614 = t34 / t4 / t27 / t18;
		    vsigmaaa[i__] = t5 * -.9305257363491 * (t9 * -.002976224 *
			     t13 + t602 * 1.91079296e-4 * t23 - t605 * 
			    2.01951168e-6 * t31 + t608 * 9.812672e-9 * t40 - 
			    t614 * 1.840566272e-11 * t298) - t50 * .03109 * 
			    t62 * (t9 * -.00677244 * t65 + t602 * 3.53128e-4 *
			     t69 - t605 * .019059632 * t73 + t608 * 
			    .0137973248 * t77 - t614 * .0019890688 * t344) + 
			    t240 * (t9 * .01005861 * t247 - t459 * 
			    4.7589966e-4 * t9 + t462 * 5.11303536e-6 * t9 - 
			    t465 * 2.73269376e-8 * t9 + t470 * 
			    6.965072064e-11 * t9);
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t660 = t114 / t84 / t107 / t98;
		    vsigmabb[i__] = t85 * -.9305257363491 * (t89 * 
			    -.002976224 * t93 + t648 * 1.91079296e-4 * t103 - 
			    t651 * 2.01951168e-6 * t111 + t654 * 9.812672e-9 *
			     t120 - t660 * 1.840566272e-11 * t506) - t130 * 
			    .03109 * t142 * (t89 * -.00677244 * t145 + t648 * 
			    3.53128e-4 * t149 - t651 * .019059632 * t153 + 
			    t654 * .0137973248 * t157 - t660 * .0019890688 * 
			    t552) + t240 * (t89 * .01005861 * t247 - t459 * 
			    4.7589966e-4 * t89 + t462 * 5.11303536e-6 * t89 - 
			    t465 * 2.73269376e-8 * t89 + t470 * 
			    6.965072064e-11 * t89);
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
		    t41 = 1.0932 - t8 * .002976224 * t11 + t19 * 8.95872e-5 * 
			    t21 - t27 * 4.3427136e-7 * t29 + t36 * 
			    1.15035392e-9 * t38;
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
		    t78 = .222601 - t8 * .00677244 * t63 - t19 * 5.0068e-4 * 
			    t67 - t27 * .006419968 * t71 + t36 * .002486336 * 
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
		    t116 = t87 * .007936597333333333 * t11 - t93 * 
			    5.095447893333333e-4 * t21 + t98 * 5.38536448e-6 *
			     t29 - t104 * 2.616712533333333e-8 * t38 + t111 * 
			    4.908176725333333e-11 * t113;
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
		    t162 = t87 * .01805984 * t63 - t93 * 9.416746666666667e-4 
			    * t67 + t98 * .05082568533333333 * t71 - t104 * 
			    .03679286613333333 * t75 + t111 * 
			    .005304183466666667 * t159;
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
		    vsigmabb[i__] = t3 * -.9305257363491 * (t7 * -.002976224 *
			     t11 + t168 * 1.91079296e-4 * t21 - t171 * 
			    2.01951168e-6 * t29 + t174 * 9.812672e-9 * t38 - 
			    t180 * 1.840566272e-11 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.00677244 * t63 + t168 * 3.53128e-4 *
			     t67 - t171 * .019059632 * t71 + t174 * 
			    .0137973248 * t75 - t180 * .0019890688 * t159);
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
			    -.02910085688888889 * t11 + t213 * 
			    .003311774037333333 * t21 - t217 * 
			    5.933856915911111e-5 * t29 + t223 * 
			    4.776147922488889e-7 * t38 - t229 * 
			    1.819969344853333e-9 * t113 + t236 * 
			    2.617694253511111e-12 * t238) + t47 * 2. * t127 * 
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
			    t207 * -.06621941333333333 * t63 + t213 * 
			    .01559585422222222 * t67 - t217 * 
			    .4584356209777778 * t71 + t223 * 
			    .5105712014222222 * t75 - t229 * .15451807744 * 
			    t159 + t236 * .01414448924444444 * t319);
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
			    2.02984192e-4 * t21 - t328 * 5.567657728e-6 * t29 
			    + t331 * 5.367215616e-8 * t38 - t334 * 
			    2.3062540288e-10 * t113 + t338 * 3.681132544e-13 *
			     t238) - t48 * .03109 * t60 * (t18 * .001707616 * 
			    t67 - t328 * .0382605152 * t71 + t331 * 
			    .0528277536 * t75 - t334 * .01899413504 * t159 + 
			    t338 * .0019890688 * t319);
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
		    t41 = 1.0932 - t8 * .002976224 * t11 + t19 * 8.95872e-5 * 
			    t21 - t27 * 4.3427136e-7 * t29 + t36 * 
			    1.15035392e-9 * t38;
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
		    t78 = .222601 - t8 * .00677244 * t63 - t19 * 5.0068e-4 * 
			    t67 - t27 * .006419968 * t71 + t36 * .002486336 * 
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
		    t116 = t87 * .007936597333333333 * t11 - t93 * 
			    5.095447893333333e-4 * t21 + t98 * 5.38536448e-6 *
			     t29 - t104 * 2.616712533333333e-8 * t38 + t111 * 
			    4.908176725333333e-11 * t113;
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
		    t162 = t87 * .01805984 * t63 - t93 * 9.416746666666667e-4 
			    * t67 + t98 * .05082568533333333 * t71 - t104 * 
			    .03679286613333333 * t75 + t111 * 
			    .005304183466666667 * t159;
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
		    vsigmaaa[i__] = t3 * -.9305257363491 * (t7 * -.002976224 *
			     t11 + t168 * 1.91079296e-4 * t21 - t171 * 
			    2.01951168e-6 * t29 + t174 * 9.812672e-9 * t38 - 
			    t180 * 1.840566272e-11 * t113) - t48 * .03109 * 
			    t60 * (t7 * -.00677244 * t63 + t168 * 3.53128e-4 *
			     t67 - t171 * .019059632 * t71 + t174 * 
			    .0137973248 * t75 - t180 * .0019890688 * t159);
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
			    -.02910085688888889 * t11 + t213 * 
			    .003311774037333333 * t21 - t217 * 
			    5.933856915911111e-5 * t29 + t223 * 
			    4.776147922488889e-7 * t38 - t229 * 
			    1.819969344853333e-9 * t113 + t236 * 
			    2.617694253511111e-12 * t238) + t47 * 2. * t127 * 
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
			    t207 * -.06621941333333333 * t63 + t213 * 
			    .01559585422222222 * t67 - t217 * 
			    .4584356209777778 * t71 + t223 * 
			    .5105712014222222 * t75 - t229 * .15451807744 * 
			    t159 + t236 * .01414448924444444 * t319);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t328 = sigmaaa * t26;
		    t331 = t14 * t35;
		    t334 = t24 * t179;
		    t338 = t32 / t233;
		    v2sigmaaa2[i__] = t3 * -.9305257363491 * (t18 * 
			    2.02984192e-4 * t21 - t328 * 5.567657728e-6 * t29 
			    + t331 * 5.367215616e-8 * t38 - t334 * 
			    2.3062540288e-10 * t113 + t338 * 3.681132544e-13 *
			     t238) - t48 * .03109 * t60 * (t18 * .001707616 * 
			    t67 - t328 * .0382605152 * t71 + t331 * 
			    .0528277536 * t75 - t334 * .01899413504 * t159 + 
			    t338 * .0019890688 * t319);
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
		    t43 = 1.0932 - t10 * .002976224 * t13 + t21 * 8.95872e-5 *
			     t23 - t29 * 4.3427136e-7 * t31 + t38 * 
			    1.15035392e-9 * t40;
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
		    t80 = .222601 - t10 * .00677244 * t65 - t21 * 5.0068e-4 * 
			    t69 - t29 * .006419968 * t73 + t38 * .002486336 * 
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
		    t123 = 1.0932 - t90 * .002976224 * t93 + t101 * 
			    8.95872e-5 * t103 - t109 * 4.3427136e-7 * t111 + 
			    t118 * 1.15035392e-9 * t120;
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
		    t160 = .222601 - t90 * .00677244 * t145 - t101 * 
			    5.0068e-4 * t149 - t109 * .006419968 * t153 + 
			    t118 * .002486336 * t157;
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
		    t265 = t243 * .02011722 * t247 + .729974 - t250 * 
			    4.15548e-4 * t252 + t255 * 1.74649824e-6 * t257 - 
			    t260 * 5.80422672e-9 * t262;
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
		    t301 = t272 * .007936597333333333 * t13 - t278 * 
			    5.095447893333333e-4 * t23 + t283 * 5.38536448e-6 
			    * t31 - t289 * 2.616712533333333e-8 * t40 + t296 *
			     4.908176725333333e-11 * t298;
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
		    t347 = t272 * .01805984 * t65 - t278 * 
			    9.416746666666667e-4 * t69 + t283 * 
			    .05082568533333333 * t73 - t289 * 
			    .03679286613333333 * t77 + t296 * 
			    .005304183466666667 * t344;
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
		    t473 = t272 * -.02682296 * t247 + t459 * .00126906576 * 
			    t272 - t462 * 1.363476096e-5 * t272 + t465 * 
			    7.28718336e-8 * t272 - t470 * 1.8573525504e-10 * 
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
		    t509 = t480 * .007936597333333333 * t93 - t486 * 
			    5.095447893333333e-4 * t103 + t491 * 
			    5.38536448e-6 * t111 - t497 * 
			    2.616712533333333e-8 * t120 + t504 * 
			    4.908176725333333e-11 * t506;
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
		    t555 = t480 * .01805984 * t145 - t486 * 
			    9.416746666666667e-4 * t149 + t491 * 
			    .05082568533333333 * t153 - t497 * 
			    .03679286613333333 * t157 + t504 * 
			    .005304183466666667 * t552;
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
		    t598 = t480 * -.02682296 * t247 + t459 * .00126906576 * 
			    t480 - t462 * 1.363476096e-5 * t480 + t465 * 
			    7.28718336e-8 * t480 - t470 * 1.8573525504e-10 * 
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
		    t617 = t9 * -.002976224 * t13 + t602 * 1.91079296e-4 * 
			    t23 - t605 * 2.01951168e-6 * t31 + t608 * 
			    9.812672e-9 * t40 - t614 * 1.840566272e-11 * t298;
		    t630 = t9 * -.00677244 * t65 + t602 * 3.53128e-4 * t69 - 
			    t605 * .019059632 * t73 + t608 * .0137973248 * 
			    t77 - t614 * .0019890688 * t344;
		    t631 = t62 * t630;
		    t644 = t9 * .01005861 * t247 - t459 * 4.7589966e-4 * t9 + 
			    t462 * 5.11303536e-6 * t9 - t465 * 2.73269376e-8 *
			     t9 + t470 * 6.965072064e-11 * t9;
		    vsigmaaa[i__] = t5 * -.9305257363491 * t617 - t50 * 
			    .03109 * t631 + t240 * t644;
		    vsigmaab[i__] = 0.;
		    t648 = sigmabb * t100;
		    t651 = t96 * t108;
		    t654 = t106 * t117;
		    t659 = 1 / t84 / t107 / t98;
		    t660 = t114 * t659;
		    t663 = t89 * -.002976224 * t93 + t648 * 1.91079296e-4 * 
			    t103 - t651 * 2.01951168e-6 * t111 + t654 * 
			    9.812672e-9 * t120 - t660 * 1.840566272e-11 * 
			    t506;
		    t676 = t89 * -.00677244 * t145 + t648 * 3.53128e-4 * t149 
			    - t651 * .019059632 * t153 + t654 * .0137973248 * 
			    t157 - t660 * .0019890688 * t552;
		    t677 = t142 * t676;
		    t690 = t89 * .01005861 * t247 - t459 * 4.7589966e-4 * t89 
			    + t462 * 5.11303536e-6 * t89 - t465 * 
			    2.73269376e-8 * t89 + t470 * 6.965072064e-11 * 
			    t89;
		    vsigmabb[i__] = t85 * -.9305257363491 * t663 - t130 * 
			    .03109 * t677 + t240 * t690;
		    t694 = sigmaaa / t7 / t17;
		    t697 = t17 * t269;
		    t700 = t16 / t4 / t697;
		    t704 = t26 / t35;
		    t710 = t34 / t7 / t27 / t17;
		    t716 = t292 / t4 / t27 / t697;
/* Computing 2nd power */
		    d__1 = t27;
		    t720 = d__1 * d__1;
		    t723 = t34 * t16 / t720 / t6;
		    t725 = 1 / t39 / t22;
		    t736 = t375 * 2.;
		    t739 = t444 * 15.38928840745229;
		    t740 = t437 * 3.847322101863073;
		    t741 = t379 * .001748158859896213;
		    t742 = t393 * 2.249999913366216;
		    t744 = t354 * .005495599555936838;
		    t751 = 1 / t361 / t165 * t210;
		    t754 = 1 / t208 / t164;
		    t755 = t362 * t754;
		    t759 = 1 / t174 / t165 * t210;
		    t761 = t351 * t754;
		    t765 = 1 / t172 / t165 * t210;
		    t767 = t366 * t754;
		    t771 = 1 / t166 / t165 * t210;
		    t773 = t369 * t754;
		    t777 = t358 * (t751 * -.8309097827459833 + t755 * 
			    1.99418347859036 - t759 * .4945709824779306 + 
			    t761 * 1.483712947433792 - t765 * 
			    .2001071587498409 + t767 * .8004286349993634 - 
			    t771 * .04215565168327908 + t773 * 
			    .2529339100996745) * t373;
		    t778 = t777 * 1.;
/* Computing 2nd power */
		    d__1 = t372;
		    t782 = d__1 * d__1;
		    t784 = t168 / t356 / t176 * t782 * t373;
		    t785 = t784 * 2.;
		    t787 = t761 * .001748158859896213 * t378;
		    t788 = t759 * t180;
		    t789 = t788 * .001831866518645613;
		    t790 = t761 * t180;
		    t791 = t790 * .005495599555936838;
		    t794 = t353 * t357 * t372 * t373;
		    t795 = t794 * .08837926660346786;
/* Computing 2nd power */
		    d__1 = t356;
		    t796 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t179;
		    t799 = d__1 * d__1;
		    t802 = t168 / t796 * t782 / t799;
		    t803 = t802 * 16.0818243221511;
/* Computing 2nd power */
		    d__1 = t428;
		    t815 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t421;
		    t831 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t226;
		    t834 = d__1 * d__1;
		    t845 = t759 * 8.806734334819047e-4 * t227 - t761 * 
			    .002642020300445714 * t227 - t353 * 
			    .08497974591333914 * t422 * t428 * t429 - t218 * 
			    2. / t421 / t223 * t815 * t429 + t423 * 1. * (
			    t751 * -1.544496508763151 + t755 * 
			    3.706791621031562 - t759 * .854388052766047 + 
			    t761 * 2.563164158298141 - t765 * 
			    .4111834438919023 + t767 * 1.644733775567609 - 
			    t771 * .05346380647307093 + t773 * 
			    .3207828388384256) * t429 + t218 * 
			    32.1646831778707 / t831 * t815 / t834 - t788 * 
			    .001831866518645613 + t790 * .005495599555936838 
			    + t794 * .08837926660346786 + t784 * 2. - t777 * 
			    1. - t802 * 16.0818243221511;
		    t848 = t845 * 1.923661050931536 * t205 * t211;
/* Computing 2nd power */
		    d__1 = t199;
		    t849 = d__1 * d__1;
		    t850 = 1 / t849;
/* Computing 2nd power */
		    d__1 = t397;
		    t851 = d__1 * d__1;
		    t854 = t352 * 2.;
		    t856 = t196 * 2. * t754;
		    t857 = -t854 + t856;
/* Computing 2nd power */
		    d__1 = t203;
		    t860 = d__1 * d__1;
		    t861 = 1 / t860;
/* Computing 2nd power */
		    d__1 = t401;
		    t862 = d__1 * d__1;
		    t868 = t850 * .4444444444444444 * t851 + t199 * 
			    1.333333333333333 * t857 + t861 * 
			    .4444444444444444 * t862 - t203 * 
			    1.333333333333333 * t857;
/* Computing 2nd power */
		    d__1 = t381;
		    t872 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t388;
		    t875 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t192;
		    t877 = d__1 * d__1;
		    t882 = t184 * 33.30964519106732 / t872 * t875 / t877 * 
			    t205 * t213;
		    t884 = t389 * t391 * t415;
		    t888 = t389 * t390 * t404 * t213;
		    t895 = t184 * 2.249999913366216 / t381 / t189 * t875 * 
			    t392;
		    t896 = t194 * .07599149707403056 * t404 * t415 + t778 - 
			    t785 + t787 + t789 - t791 - t795 + t803 + t848 + 
			    t194 * .03799574853701528 * t868 * t213 - t882 - 
			    t884 * 2.249999913366216 - t888 * 
			    2.249999913366216 + t895;
		    t908 = t383 * 1.124999956683108 * (t751 * 
			    -1.132974264373283 + t755 * 2.71913823449588 - 
			    t759 * .4994648585728036 + t761 * 
			    1.498394575718411 - t765 * .1075243117819161 + 
			    t767 * .4300972471276643 - t771 * 
			    .04247805766949639 + t773 * .2548683460169784) * 
			    t392;
		    t910 = t759 * 5.827196199654043e-4 * t378;
		    t913 = t353 * t193 * t404 * t213;
		    t916 = t353 * t377 * t415;
		    t922 = t353 * .05176049209143758 * t382 * t388 * t390 * 
			    t214;
		    t923 = t436 * t409;
		    t926 = t436 * 15.38928840745229 * t413;
		    t928 = t435 * t404 * t211;
		    t930 = t206 * t210;
		    t931 = t231 * t930;
		    t932 = t931 * 23.08393261117844;
		    t933 = t408 * t412;
		    t934 = t231 * t933;
		    t936 = t930 * 12.;
		    t937 = t933 * 32.;
		    t940 = t207 / t209 / t208;
		    t941 = t940 * 20.;
		    t946 = t439 * t413;
		    t948 = t439 * t409;
		    t954 = t231 * 38.47322101863073 * t940;
		    t955 = -t908 - t910 - t913 * .001748158859896213 - t916 * 
			    .001748158859896213 + t922 + t923 * 
			    15.38928840745229 - t926 + t928 * 
			    3.847322101863073 + t932 - t934 * 
			    61.55715362980916 + t194 * .03799574853701528 * 
			    t205 * (-t936 + t937 - t941) - t946 * 
			    15.38928840745229 + t948 * 15.38928840745229 + 
			    t230 * 1.923661050931536 * t868 * t211 + t954;
		    t960 = 1 / t17;
		    t963 = 1 / t269;
		    t967 = 1 / t56 / t46;
		    t984 = -1.544496508763151 / t316 / t46 * t960 + t317 * 
			    3.706791621031562 * t963 - t967 * 
			    .854388052766047 * t960 + t307 * 
			    2.563164158298141 * t963 - .4111834438919023 / 
			    t54 / t46 * t960 + t323 * 1.644733775567609 * 
			    t963 - .05346380647307093 / t47 / t46 * t960 + 
			    t326 * .3207828388384256 * t963;
		    t990 = 1 / t311 / t58;
/* Computing 2nd power */
		    d__1 = t329;
		    t991 = d__1 * d__1;
		    t996 = t963 * t967;
		    t1001 = t49 * t312;
/* Computing 2nd power */
		    d__1 = t311;
		    t1004 = d__1 * d__1;
		    t1005 = 1 / t1004;
/* Computing 2nd power */
		    d__1 = t61;
		    t1007 = d__1 * d__1;
		    t1008 = 1 / t1007;
		    t1012 = t736 + t440 * 3.847322101863073 + t442 * 
			    15.38928840745229 - t739 + t740 - t741 - t742 + 
			    t406 * .07599149707403056 + t744 + t417 * 
			    .07599149707403056 + t164 * (t896 + t955) - t50 * 
			    1. * t312 * t984 * t330 + t50 * 2. * t990 * t991 *
			     t330 - t996 * 8.806734334819047e-4 * t62 + t308 *
			     .08497974591333914 * t452 - t1001 * 2. * t331 - 
			    t50 * 32.1646831778707 * t1005 * t991 * t1008;
		    t1031 = 1 / t76 / t68;
		    t1066 = t243 * t257;
		    t1071 = t250 * t262;
		    t1076 = t255 * t469;
		    t1083 = t260 / t261 / t251;
		    s1 = t5 * -.9305257363491 * (t694 * -.02910085688888889 * 
			    t13 + t700 * .003311774037333333 * t23 - t704 * 
			    5.933856915911111e-5 * t31 + t710 * 
			    4.776147922488889e-7 * t40 - t716 * 
			    1.819969344853333e-9 * t298 + t723 * 
			    2.617694253511111e-12 * t725) - .4135669939329333 
			    / t7 * t43 - t4 * 2.4814019635976 * t301 + t1012 *
			     t265 + t308 * .002642020300445714 * t348 + t1001 
			    * 2. * t332 - t304 * .06218 * t347 - t50 * .03109 
			    * t62 * (t694 * -.06621941333333333 * t65 + t700 *
			     .01559585422222222 * t69 - t704 * 
			    .4584356209777778 * t73 + t710 * 
			    .5105712014222222 * t77 - t716 * .15451807744 * 
			    t344 + t723 * .01414448924444444 * t1031);
		    v2rhoa2[i__] = s1 + t313 * 1. * t984 * t330 * t80 - t50 * 
			    2. * t990 * t991 * t330 * t80 - t308 * 
			    .08497974591333914 * t312 * t332 + t996 * 
			    8.806734334819047e-4 * t81 + t50 * 
			    32.1646831778707 * t1005 * t991 * t1008 * t80 + 
			    t313 * 2. * t331 * t347 + t455 * 2. * t473 + t240 
			    * (t694 * .09835085333333333 * t247 - t700 * 
			    .00190667136 * t252 + t1066 * 5.666441472e-5 * 
			    t700 - t459 * .00465324112 * t694 - t1071 * 
			    6.1872159744e-7 * t700 + t462 * 4.999412352e-5 * 
			    t694 + t1076 * 3.32248670208e-9 * t700 - t465 * 
			    2.671967232e-7 * t694 - t1083 * 7.4294102016e-12 *
			     t700 + t470 * 6.8102926848e-10 * t694);
		    t1092 = sigmabb / t87 / t97;
		    t1095 = t97 * t477;
		    t1098 = t96 / t84 / t1095;
		    t1102 = t106 / t115;
		    t1108 = t114 / t87 / t107 / t97;
		    t1114 = t500 / t84 / t107 / t1095;
/* Computing 2nd power */
		    d__1 = t107;
		    t1118 = d__1 * d__1;
		    t1121 = t114 * t96 / t1118 / t86;
		    t1123 = 1 / t119 / t102;
/* Computing 2nd power */
		    d__1 = t519;
		    t1134 = d__1 * d__1;
		    t1135 = 1 / t1134;
/* Computing 2nd power */
		    d__1 = t537;
		    t1137 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t141;
		    t1138 = d__1 * d__1;
		    t1139 = 1 / t1138;
		    t1146 = 1 / t97;
		    t1149 = 1 / t477;
		    t1153 = 1 / t136 / t126;
		    t1170 = -1.544496508763151 / t524 / t126 * t1146 + t525 * 
			    3.706791621031562 * t1149 - t1153 * 
			    .854388052766047 * t1146 + t515 * 
			    2.563164158298141 * t1149 - .4111834438919023 / 
			    t134 / t126 * t1146 + t531 * 1.644733775567609 * 
			    t1149 - .05346380647307093 / t127 / t126 * t1146 
			    + t534 * .3207828388384256 * t1149;
		    t1176 = 1 / t519 / t138;
		    t1184 = t1149 * t1153;
		    t1190 = t129 * t520;
/* Computing 2nd power */
		    d__1 = t559;
		    t1206 = d__1 * d__1;
		    t1209 = t854 + t856;
/* Computing 2nd power */
		    d__1 = t562;
		    t1212 = d__1 * d__1;
		    t1218 = t850 * .4444444444444444 * t1206 + t199 * 
			    1.333333333333333 * t1209 + t861 * 
			    .4444444444444444 * t1212 - t203 * 
			    1.333333333333333 * t1209;
		    t1227 = t353 * t193 * t565 * t213;
		    t1230 = t353 * t377 * t569;
		    t1232 = t573 * t409;
		    t1234 = t573 * t413;
		    t1240 = t435 * t565 * t211;
		    t1242 = t194 * .03799574853701528 * t205 * (-t936 - t937 
			    - t941) + t194 * .03799574853701528 * t1218 * 
			    t213 + t194 * .07599149707403056 * t565 * t569 - 
			    t1227 * .001748158859896213 - t1230 * 
			    .001748158859896213 - t1232 * 15.38928840745229 - 
			    t1234 * 15.38928840745229 + t230 * 
			    1.923661050931536 * t1218 * t211 + t1240 * 
			    3.847322101863073 + t778 - t785 + t787 + t789 - 
			    t791;
		    t1245 = t389 * t390 * t565 * t213;
		    t1248 = t389 * t391 * t569;
		    t1252 = -t795 + t803 + t848 - t1245 * 2.249999913366216 - 
			    t882 + t895 - t908 - t910 - t1248 * 
			    2.249999913366216 + t922 - t923 * 
			    15.38928840745229 - t926 + t932 + t934 * 
			    61.55715362980916 + t954;
		    t1273 = t736 - t442 * 15.38928840745229 - t739 + t740 - 
			    t741 - t742 + t744 + t574 * 3.847322101863073 + 
			    t571 * .07599149707403056 + t567 * 
			    .07599149707403056 + t164 * (t1242 + t1252) - 
			    t1190 * 2. * t539 - t130 * 32.1646831778707 * 
			    t1135 * t1137 * t1139 - t130 * 1. * t520 * t1170 *
			     t538 + t130 * 2. * t1176 * t1137 * t538 - t1184 *
			     8.806734334819047e-4 * t142 + t516 * 
			    .08497974591333914 * t583;
		    t1286 = 1 / t156 / t148;
		    s1 = t85 * -.9305257363491 * (t1092 * -.02910085688888889 
			    * t93 + t1098 * .003311774037333333 * t103 - 
			    t1102 * 5.933856915911111e-5 * t111 + t1108 * 
			    4.776147922488889e-7 * t120 - t1114 * 
			    1.819969344853333e-9 * t506 + t1121 * 
			    2.617694253511111e-12 * t1123) - 
			    .4135669939329333 / t87 * t123 - t84 * 
			    2.4814019635976 * t509 + t130 * 32.1646831778707 *
			     t1135 * t1137 * t1139 * t160 + t521 * 1. * t1170 
			    * t538 * t160 - t130 * 2. * t1176 * t1137 * t538 *
			     t160 + t516 * .002642020300445714 * t556 + t1184 
			    * 8.806734334819047e-4 * t161;
		    v2rhob2[i__] = s1 - t516 * .08497974591333914 * t520 * 
			    t540 + t1190 * 2. * t540 + t521 * 2. * t539 * 
			    t555 - t512 * .06218 * t555 + t1273 * t265 - t130 
			    * .03109 * t142 * (t1092 * -.06621941333333333 * 
			    t145 + t1098 * .01559585422222222 * t149 - t1102 *
			     .4584356209777778 * t153 + t1108 * 
			    .5105712014222222 * t157 - t1114 * .15451807744 * 
			    t552 + t1121 * .01414448924444444 * t1286) + t240 
			    * (t1092 * .09835085333333333 * t247 - t1098 * 
			    .00190667136 * t252 + t1066 * 5.666441472e-5 * 
			    t1098 - t459 * .00465324112 * t1092 - t1071 * 
			    6.1872159744e-7 * t1098 + t462 * 4.999412352e-5 * 
			    t1092 + t1076 * 3.32248670208e-9 * t1098 - t465 * 
			    2.671967232e-7 * t1092 - t1083 * 7.4294102016e-12 
			    * t1098 + t470 * 6.8102926848e-10 * t1092) + t586 
			    * 2. * t598;
		    t1325 = t1227 * -8.740794299481065e-4 - t1230 * 
			    8.740794299481065e-4 + t1232 * 7.694644203726145 
			    - t1234 * 7.694644203726145 + t1240 * 
			    1.923661050931536 + t778 - t785 + t787 + t789 - 
			    t791 - t795 + t803 + t848 - t1245 * 
			    1.124999956683108 - t882 - t884 * 
			    1.124999956683108 - t888 * 1.124999956683108;
		    t1355 = t850 * .4444444444444444 * t397 * t559 + t199 * 
			    2.666666666666667 * t196 * t754 + t861 * 
			    .4444444444444444 * t401 * t562 - t203 * 
			    2.666666666666667 * t196 * t754;
		    t1362 = t895 - t908 - t910 - t1248 * 1.124999956683108 - 
			    t913 * 8.740794299481065e-4 - t916 * 
			    8.740794299481065e-4 + t922 - t926 + t928 * 
			    1.923661050931536 - t931 * 23.08393261117844 - 
			    t946 * 7.694644203726145 - t948 * 
			    7.694644203726145 + t954 + t194 * 
			    .03799574853701528 * t404 * t569 + t194 * 
			    .03799574853701528 * t205 * (t936 - t941) + t194 *
			     .03799574853701528 * t565 * t415 + t194 * 
			    .03799574853701528 * t1355 * t213 + t230 * 
			    1.923661050931536 * t1355 * t211;
		    t1365 = t744 + t736 - t741 - t742 + t568 + t572 + t740 + 
			    t575 - t739 + t407 + t418 + t441 + t164 * (t1325 
			    + t1362);
		    t1378 = t271 * sigmabb * t479;
		    v2rhoab[i__] = t1365 * 1. * t265 + t455 * 1. * t598 + 
			    t586 * 1. * t473 + t240 * 1. * (t272 * 
			    -.00190667136 * t252 * sigmabb * t479 + t1066 * 
			    5.666441472e-5 * sigmaaa * t1378 - t1071 * 
			    6.1872159744e-7 * sigmaaa * t1378 + t1076 * 
			    3.32248670208e-9 * sigmaaa * t1378 - t1083 * 
			    7.4294102016e-12 * sigmaaa * t1378);
		    t1397 = sigmaaa * t277;
		    t1400 = t16 * t282;
		    t1403 = t26 * t288;
		    t1406 = t34 * t295;
		    t1411 = t292 / t720 / rhoa;
		    s1 = t4 * -1.2407009817988 * t617 - t5 * .9305257363491 * 
			    (t271 * .007936597333333333 * t13 - t1397 * 
			    .001050835968 * t23 + t1400 * 
			    2.023245175466667e-5 * t31 - t1403 * 
			    1.692928750933333e-7 * t40 + t1406 * 
			    6.640828416e-10 * t298 - t1411 * 
			    9.816353450666667e-13 * t725) - t304 * .03109 * 
			    t630 + t308 * .001321010150222857 * t631;
		    v2rhoasigmaaa[i__] = s1 + t313 * 1. * t331 * t630 - t50 * 
			    .03109 * t62 * (t271 * .01805984 * t65 - t1397 * 
			    .005495317333333333 * t69 + t1400 * 
			    .1528537258666667 * t73 - t1403 * 
			    .1776668757333333 * t77 + t1406 * .05595521024 * 
			    t344 - t1411 * .005304183466666667 * t1031) + 
			    t455 * t644 + t240 * (t271 * -.02682296 * t247 + 
			    t1397 * 7.1500176e-4 * t252 - t1066 * 
			    2.124915552e-5 * t1397 + t459 * .00126906576 * 
			    t271 + t1071 * 2.3202059904e-7 * t1397 - t462 * 
			    1.363476096e-5 * t271 - t1076 * 1.24593251328e-9 *
			     t1397 + t465 * 7.28718336e-8 * t271 + t1083 * 
			    2.7860288256e-12 * t1397 - t470 * 
			    1.8573525504e-10 * t271);
		    v2rhoasigmaab[i__] = 0.;
		    t1467 = t272 * t89;
		    v2rhoasigmabb[i__] = t455 * t690 + t240 * (t272 * 
			    7.1500176e-4 * t252 * t89 - t1066 * 
			    2.124915552e-5 * t1467 + t1071 * 2.3202059904e-7 *
			     t1467 - t1076 * 1.24593251328e-9 * t1467 + t1083 
			    * 2.7860288256e-12 * t1467);
		    t1479 = t252 * t9;
		    t1482 = t480 * t9;
		    v2rhobsigmaaa[i__] = t586 * t644 + t240 * (t480 * 
			    7.1500176e-4 * t1479 - t1066 * 2.124915552e-5 * 
			    t1482 + t1071 * 2.3202059904e-7 * t1482 - t1076 * 
			    1.24593251328e-9 * t1482 + t1083 * 
			    2.7860288256e-12 * t1482);
		    v2rhobsigmaab[i__] = 0.;
		    t1497 = sigmabb * t485;
		    t1500 = t96 * t490;
		    t1503 = t106 * t496;
		    t1506 = t114 * t503;
		    t1511 = t500 / t1118 / rhob;
		    s1 = t84 * -1.2407009817988 * t663 - t85 * .9305257363491 
			    * (t479 * .007936597333333333 * t93 - t1497 * 
			    .001050835968 * t103 + t1500 * 
			    2.023245175466667e-5 * t111 - t1503 * 
			    1.692928750933333e-7 * t120 + t1506 * 
			    6.640828416e-10 * t506 - t1511 * 
			    9.816353450666667e-13 * t1123) - t512 * .03109 * 
			    t676 + t516 * .001321010150222857 * t677;
		    v2rhobsigmabb[i__] = s1 + t521 * 1. * t539 * t676 - t130 *
			     .03109 * t142 * (t479 * .01805984 * t145 - t1497 
			    * .005495317333333333 * t149 + t1500 * 
			    .1528537258666667 * t153 - t1503 * 
			    .1776668757333333 * t157 + t1506 * .05595521024 * 
			    t552 - t1511 * .005304183466666667 * t1286) + 
			    t586 * t690 + t240 * (t479 * -.02682296 * t247 + 
			    t1497 * 7.1500176e-4 * t252 - t1066 * 
			    2.124915552e-5 * t1497 + t459 * .00126906576 * 
			    t479 + t1071 * 2.3202059904e-7 * t1497 - t462 * 
			    1.363476096e-5 * t479 - t1076 * 1.24593251328e-9 *
			     t1497 + t465 * 7.28718336e-8 * t479 + t1083 * 
			    2.7860288256e-12 * t1497 - t470 * 
			    1.8573525504e-10 * t479);
		    t1565 = sigmaaa * t28;
		    t1568 = t16 * t37;
		    t1571 = t26 * t613;
		    t1575 = t34 / t720;
		    v2sigmaaa2[i__] = t5 * -.9305257363491 * (t20 * 
			    2.02984192e-4 * t23 - t1565 * 5.567657728e-6 * 
			    t31 + t1568 * 5.367215616e-8 * t40 - t1571 * 
			    2.3062540288e-10 * t298 + t1575 * 3.681132544e-13 
			    * t725) - t50 * .03109 * t62 * (t20 * .001707616 *
			     t69 - t1565 * .0382605152 * t73 + t1568 * 
			    .0528277536 * t77 - t1571 * .01899413504 * t344 + 
			    t1575 * .0019890688 * t1031) + t240 * (t20 * 
			    -2.6812566e-4 * t252 + t1066 * 7.96843332e-6 * 
			    t20 - t1071 * 8.700772464e-8 * t20 + t1076 * 
			    4.6722469248e-10 * t20 - t1083 * 1.0447608096e-12 
			    * t20);
		    v2sigmaaaab[i__] = 0.;
		    t1609 = t9 * t89;
		    v2sigmaaabb[i__] = t240 * (t1479 * -2.6812566e-4 * t89 + 
			    t1066 * 7.96843332e-6 * t1609 - t1071 * 
			    8.700772464e-8 * t1609 + t1076 * 4.6722469248e-10 
			    * t1609 - t1083 * 1.0447608096e-12 * t1609);
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
		    t1621 = sigmabb * t108;
		    t1624 = t96 * t117;
		    t1627 = t106 * t659;
		    t1631 = t114 / t1118;
		    v2sigmabb2[i__] = t85 * -.9305257363491 * (t100 * 
			    2.02984192e-4 * t103 - t1621 * 5.567657728e-6 * 
			    t111 + t1624 * 5.367215616e-8 * t120 - t1627 * 
			    2.3062540288e-10 * t506 + t1631 * 3.681132544e-13 
			    * t1123) - t130 * .03109 * t142 * (t100 * 
			    .001707616 * t149 - t1621 * .0382605152 * t153 + 
			    t1624 * .0528277536 * t157 - t1627 * .01899413504 
			    * t552 + t1631 * .0019890688 * t1286) + t240 * (
			    t100 * -2.6812566e-4 * t252 + t1066 * 
			    7.96843332e-6 * t100 - t1071 * 8.700772464e-8 * 
			    t100 + t1076 * 4.6722469248e-10 * t100 - t1083 * 
			    1.0447608096e-12 * t100);
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
} /* uks_xc_hcth__ */

/* Subroutine */ int rks_xc_hcth__(integer *ideriv, integer *npt, doublereal *
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
	    t240, t241, t256, rho, t260, t273, t274, t287, t292, t295, t298, 
	    t299, t302, t303, t308, t309, t314, t315, t318, t321, t323, t325, 
	    t348, t356, t370, t377, t378, t380, t381, t382, t389, t390, t392, 
	    t393, t396, t397, t399, t403, t405, t409, t411, t413, t423, t425, 
	    t429, t433, t436, t448, t449, t452, t456, t467, t469, t473, t474, 
	    t477, t481, t525, t528, t531, t534, t539, t572, t574, t576, t578, 
	    t581, t592, t595, t598, t602, sigma;


/*     F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N.C. Handy */
/*     Development and assessment of new exchange-correlation functionals */
/*     J. Chem. Phys. 109 (1998) 6264-6271 */


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
		zk[i__] = t2 * -.7385587663820224 * rho * (1.0932 - t8 * 
			.004724461108493003 / t10 + t19 * 2.25745598162284e-4 
			/ t20 - t27 * 1.73708544e-6 / t20 / t10 + t36 * 
			7.304292090974968e-9 / t37) - t48 * .03109 * t60 * (
			.222601 - t8 * .01075057838039151 / t62 - t19 * 
			.00126163454252273 / t66 - t27 * .025679872 / t66 / 
			t62 + t36 * .01578724952778562 / t74) + (rho * 
			-.062182 * (t45 * .1325688999052018 + 1.) * t93 + t48 
			* .03109 * t60) * (t8 * .0319340961906757 / t100 + 
			.729974 - t19 * .00104711534488343 / t104 + t27 * 
			6.98599296e-6 / t104 / t100 - t36 * 
			3.685454240475973e-8 / t112);
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
		t41 = 1.0932 - t8 * .004724461108493003 * t11 + t19 * 
			2.25745598162284e-4 * t21 - t27 * 1.73708544e-6 * t29 
			+ t36 * 7.304292090974968e-9 * t38;
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
		t78 = .222601 - t8 * .01075057838039151 * t63 - t19 * 
			.00126163454252273 * t67 - t27 * .025679872 * t71 + 
			t36 * .01578724952778562 * t75;
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
		t116 = t8 * .0319340961906757 * t101 + .729974 - t19 * 
			.00104711534488343 * t105 + t27 * 6.98599296e-6 * 
			t109 - t36 * 3.685454240475973e-8 * t113;
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
			(t123 * .02519712591196268 * t11 - t129 * 
			.002567944823781261 * t21 + t134 * 4.308291584e-5 * 
			t29 - t140 * 3.323017782489364e-7 * t38 + t147 * 
			9.894264276562486e-10 * t149) - t155 * .03109 * t78 + 
			t159 * .001664368495390566 * t79;
		vrhoa[i__] = s1 + t48 * .5 * t163 * t180 * t181 * t78 - t48 * 
			.015545 * t60 * (t123 * .05733641802875474 * t63 - 
			t129 * .004745742938744286 * t67 + t134 * 
			.4066054826666667 * t71 - t140 * .4672402752398278 * 
			t75 + t147 * .1069256384345231 * t195) + (t83 * 
			-.062182 * t93 + rho * (t172 * .002747799777968419 * 
			t93 + t83 * 1. / t206 * (t170 * -.99709173929518 - 
			t172 * .7418564737168958 - t175 * .4002143174996817 - 
			t178 * .1264669550498372) / t92) + t155 * .03109 - 
			t159 * .001664368495390566 * t60 - t48 * .5 * t163 * 
			t180 * t181) * t116 + t98 * (t123 * 
			-.08515758984180187 * t101 + t129 * 
			.006395690658899341 * t105 - t134 * 1.0907808768e-4 * 
			t109 + t140 * 9.254146025239327e-7 * t113 - t147 * 
			3.744188120519821e-9 * t238);
		t245 = sigma * t18;
		t248 = t14 * t26;
		t251 = t24 * t35;
		t257 = t32 / t2 / t25 / t16;
		vsigmaaa[i__] = t3 * -.7385587663820224 * (t7 * 
			-.01889784443397201 * t11 + t245 * 
			.001925958617835946 * t21 - t248 * 3.231218688e-5 * 
			t29 + t251 * 2.492263336867023e-7 * t38 - t257 * 
			7.420698207421865e-10 * t149) - t48 * .03109 * t60 * (
			t7 * -.04300231352156605 * t63 + t245 * 
			.003559307204058214 * t67 - t248 * .304954112 * t71 + 
			t251 * .3504302064298708 * t75 - t257 * 
			.08019422882589234 * t195) + t98 * 2. * (t7 * 
			.0638681923813514 * t101 - t245 * .004796767994174505 
			* t105 + t248 * 8.180856576e-5 * t109 - t251 * 
			6.940609518929495e-7 * t113 + t257 * 
			2.808141090389866e-9 * t238);
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
		t41 = 1.0932 - t8 * .004724461108493003 * t11 + t19 * 
			2.25745598162284e-4 * t21 - t27 * 1.73708544e-6 * t29 
			+ t36 * 7.304292090974968e-9 * t38;
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
		t78 = .222601 - t8 * .01075057838039151 * t63 - t19 * 
			.00126163454252273 * t67 - t27 * .025679872 * t71 + 
			t36 * .01578724952778562 * t75;
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
		t116 = t8 * .0319340961906757 * t101 + .729974 - t19 * 
			.00104711534488343 * t105 + t27 * 6.98599296e-6 * 
			t109 - t36 * 3.685454240475973e-8 * t113;
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
		t152 = t123 * .02519712591196268 * t11 - t129 * 
			.002567944823781261 * t21 + t134 * 4.308291584e-5 * 
			t29 - t140 * 3.323017782489364e-7 * t38 + t147 * 
			9.894264276562486e-10 * t149;
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
		t198 = t123 * .05733641802875474 * t63 - t129 * 
			.004745742938744286 * t67 + t134 * .4066054826666667 *
			 t71 - t140 * .4672402752398278 * t75 + t147 * 
			.1069256384345231 * t195;
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
		t230 = t123 * .08515758984180187 * t101;
		t232 = t129 * .006395690658899341 * t105;
		t234 = t134 * 1.0907808768e-4 * t109;
		t235 = t140 * t113;
		t238 = 1 / t112 / t100;
		t240 = t147 * 3.744188120519821e-9 * t238;
		t241 = -t230 + t232 - t234 + t235 * 9.254146025239327e-7 - 
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
		t260 = t7 * -.01889784443397201 * t11 + t245 * 
			.001925958617835946 * t21 - t248 * 3.231218688e-5 * 
			t29 + t251 * 2.492263336867023e-7 * t38 - t257 * 
			7.420698207421865e-10 * t149;
		t273 = t7 * -.04300231352156605 * t63 + t245 * 
			.003559307204058214 * t67 - t248 * .304954112 * t71 + 
			t251 * .3504302064298708 * t75 - t257 * 
			.08019422882589234 * t195;
		t274 = t60 * t273;
		t287 = t7 * .0638681923813514 * t101 - t245 * 
			.004796767994174505 * t105 + t248 * 8.180856576e-5 * 
			t109 - t251 * 6.940609518929495e-7 * t113 + t257 * 
			2.808141090389866e-9 * t238;
		vsigmaaa[i__] = t3 * -.7385587663820224 * t260 - t48 * .03109 
			* t274 + t98 * 2. * t287;
		t292 = sigma / t5 / t15;
		t295 = t15 * t120;
		t298 = t14 / t2 / t295;
		t299 = t298 * t105;
		t302 = t24 / t33;
		t303 = t302 * t109;
		t308 = t32 / t5 / t25 / t15;
		t309 = t308 * t113;
		t314 = t143 / t2 / t25 / t295;
		t315 = t314 * t238;
/* Computing 2nd power */
		d__1 = t25;
		t318 = d__1 * d__1;
		t321 = t32 * t14 / t318 / t4;
		t323 = 1 / t112 / t104;
		t325 = t321 * 4.754822529024e-10 * t323;
		t348 = 1 / t37 / t20;
		t356 = t47 * t163;
		t370 = 1 / t74 / t66;
/* Computing 2nd power */
		d__1 = t162;
		t377 = d__1 * d__1;
		t378 = 1 / t377;
/* Computing 2nd power */
		d__1 = t180;
		t380 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t59;
		t381 = d__1 * d__1;
		t382 = 1 / t381;
		t389 = 1 / t15;
		t390 = 1 / t167 / t44 * t389;
		t392 = 1 / t120;
		t393 = t168 * t392;
		t396 = 1 / t54 / t44;
		t397 = t396 * t389;
		t399 = t158 * t392;
		t403 = 1 / t52 / t44 * t389;
		t405 = t174 * t392;
		t409 = 1 / t45 / t44 * t389;
		t411 = t177 * t392;
		t413 = t390 * -6.934554859331846 + t393 * 16.64293166239643 - 
			t397 * 4.305845969834537 + t399 * 12.91753790950361 - 
			t403 * 2.326004811900819 + t405 * 9.304019247603276 - 
			t409 * .3394740105503081 + t411 * 2.036844063301848;
		t423 = t397 * .001831866518645613 * t93;
		t425 = t399 * .005495599555936838 * t93;
		t429 = t172 * .08837926660346786 * t207 * t213 * t214;
/* Computing 2nd power */
		d__1 = t213;
		t433 = d__1 * d__1;
		t436 = t83 * 2. / t206 / t89 * t433 * t214;
		t448 = t208 * 1. * (t390 * -.8309097827459833 + t393 * 
			1.99418347859036 - t397 * .4945709824779306 + t399 * 
			1.483712947433792 - t403 * .2001071587498409 + t405 * 
			.8004286349993634 - t409 * .04215565168327908 + t411 *
			 .2529339100996745) * t214;
/* Computing 2nd power */
		d__1 = t206;
		t449 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t92;
		t452 = d__1 * d__1;
		t456 = t83 * 16.0818243221511 / t449 * t433 / t452;
		t467 = log(29.60857464321668 / (t49 * 8.157414703487641 + t45 
			* 2.247591863577616 + t52 * .4300972471276643 + t54 * 
			.1911512595127338) + 1.);
		t469 = (t45 * .06901399211255825 + 1.) * t467 * t169;
		t473 = t216 * 2.;
		t474 = t204 * .005495599555936838;
		t477 = t392 * t396;
		t481 = 1 / t162 / t56;
		s1 = t98 * (t292 * .6244889921732137 * t101 - t299 * 
			.06611977455216065 + t303 * .00170653661184 - t309 * 
			2.250092278626939e-5 + t315 * 1.614116494367631e-7 - 
			t325) + t227 * 3. * t241 + t159 * .003328736990781132 
			* t199 - .6564966812284644 / t5 * t41 - t2 * 
			1.969490043685393 * t152 - t3 * .3692793831910112 * (
			t292 * -.184778923354393 * t11 + t298 * 
			.03338059057705277 * t21 - t302 * 
			9.494171065457778e-4 * t29 + t308 * 
			1.213065957842335e-5 * t38 - t314 * 
			7.337656600781108e-8 * t149 + t321 * 
			1.675324322247111e-10 * t348) - t155 * .06218 * t198 
			+ t356 * 2. * t183 - t48 * .015545 * t60 * (t292 * 
			-.4204670655442014 * t63 + t298 * .1571963602053569 * 
			t67 - t302 * 7.334969935644444 * t71 + t308 * 
			12.96770019587685 * t75 - t314 * 6.229778507390148 * 
			t195 + t321 * .9052473116444444 * t370);
		s2 = s1 + t48 * 16.08234158893535 * t378 * t380 * t382 * t78 
			+ (t48 * -.5 * t163 * t413 * t181 - t48 * 
			16.08234158893535 * t378 * t380 * t382 + rho * (t423 
			- t425 - t429 - t436 + t448 + t456 + t469 * 
			.03377399869956914) + t473 + t474 - t356 * 2. * t182 
			- t477 * .002219157993854088 * t60 + t48 * 1. * t481 *
			 t380 * t181 + t159 * .1070677706909338 * t224) * 
			t116 + t164 * 1. * t182 * t198 + t477 * 
			.002219157993854088 * t79;
		v2rhoa2[i__] = s2 - t159 * .1070677706909338 * t163 * t183 - 
			t48 * 1. * t481 * t380 * t181 * t78 + t164 * .5 * 
			t413 * t181 * t78 + (t474 + t473 + rho * (t423 - t425 
			- t429 - t436 + t448 + t456 - t469 * 
			.03377399869956914)) * t116 + t227 * (-t230 + t232 - 
			t234 + t235 * 9.254146025239327e-7 - t240) + t98 * (
			t299 * -.01921804305356549 + t303 * 9.0663063552e-4 - 
			t309 * 1.571454903442721e-5 + t315 * 
			1.339542698862844e-7 - t325);
		t525 = sigma * t128;
		t528 = t14 * t133;
		t531 = t24 * t139;
		t534 = t32 * t146;
		t539 = t143 / t318 / rho;
		t572 = t525 * t105;
		t574 = t528 * t109;
		t576 = t531 * t113;
		t578 = t534 * t238;
		t581 = t539 * 3.566116896768e-10 * t323;
		s1 = t2 * -.9847450218426965 * t260 - t3 * .3692793831910112 *
			 (t122 * .1007885036478507 * t11 - t525 * 
			.02118352569711769 * t21 + t528 * 
			6.474384561493333e-4 * t29 - t531 * 
			8.599542016444107e-6 * t38 + t534 * 
			5.354828486437394e-8 * t149 - t539 * 
			1.256493241685333e-10 * t348) - t155 * .03109 * t273 
			+ t159 * .001664368495390566 * t274;
		v2rhoasigmaaa[i__] = s1 + t164 * .5 * t182 * t273 - t48 * 
			.015545 * t60 * (t122 * .2293456721150189 * t63 - 
			t525 * .1107786557459012 * t67 + t528 * 
			4.891319227733333 * t71 - t531 * 9.024914734047895 * 
			t75 + t534 * 4.511945422890826 * t195 - t539 * 
			.6789354837333333 * t370) + t227 * 2. * t287 + t98 * (
			t122 * -.3406303593672075 * t101 + t572 * 
			.03999629492577148 - t574 * .00111628532736 + t576 * 
			1.548757018591614e-5 - t578 * 1.154424548967926e-7 + 
			t581) + t98 * (t572 * .01441353229017411 - t574 * 
			6.7997297664e-4 + t576 * 1.178591177582041e-5 - t578 *
			 1.004657024147133e-7 + t581);
		t592 = sigma * t26;
		t595 = t14 * t35;
		t598 = t24 * t256;
		t602 = t32 / t318;
		v2sigmaaa2[i__] = t3 * -.7385587663820224 * (t18 * 
			.00818380980149448 * t21 - t592 * 3.56330094592e-4 * 
			t29 + t595 * 5.452751177586271e-6 * t38 - t598 * 
			3.719293436531171e-8 * t149 + t602 * 
			9.42369931264e-11 * t348) - t48 * .03109 * t60 * (t18 
			* .06884676299319308 * t67 - t592 * 2.4486729728 * 
			t71 + t595 * 5.366965224816438 * t75 - t598 * 
			3.06318215186455 * t195 + t602 * .5092016128 * t370) 
			+ t98 * 4. * (t18 * -.01081014921763059 * t105 + t592 
			* 5.0997973248e-4 * t109 - t595 * 
			8.839433831865308e-6 * t113 + t598 * 
			7.534927681103499e-8 * t238 - t602 * 
			2.674587672576e-10 * t323);
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
} /* rks_xc_hcth__ */


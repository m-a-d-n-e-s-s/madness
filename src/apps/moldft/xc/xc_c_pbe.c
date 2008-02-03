/* xc_c_pbe.f -- translated by f2c (version 20050501).
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

/* :C_PBEsubrstart */
/*    Generated: Wed Sep  3 12:41:37 GMT 2003 */
/* Subroutine */ int uks_c_pbe__(integer *ideriv, integer *npt, doublereal *
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
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t4, t5, t6, t7, t8, t9, s3, t11, t20, 
	    t21, t23, t31, t33, t17, t18, t26, t27, t35, t37, t38, t49, t13, 
	    t19, t32, t36, t41, t42, t44, t45, t46, t47, t48, t50, t55, t66, 
	    t72, t73, t75, t77, t78, t79, t80, t82, t84, t89, t90, t95, t98, 
	    t16, t24, t28, t29, t34, t100, t40, t102, t43, t104, t54, t60, 
	    t116, t62, t65, t67, t83, t86, t105, t109, t114, t15, t52, t57, 
	    t69, t70, t81, t85, t91, t92, t93, t96, t99, t101, t103, t106, 
	    t107, t110, t111, t115, t118, t124, t127, t128, t131, t132, t135, 
	    t137, t138, t139, t140, t143, t146, t150, t151, t155, t156, t160, 
	    t163, t167, t168, t171, t177, t178, t180, t184, t191, t201, t207, 
	    t208, t209, t210, t211, rho, t214, t215, t218, t222, t225, t230, 
	    t231, t235, t236, t249, t256, t259, t261, t262, t265, t268, t285, 
	    t288, t291, t298, t306, t316, t320, t323, t326, t352, t353, t354, 
	    t357, t358, t361, t362, t363, t374, t378, t381, t423, t61, t76, 
	    t87, t88, t120, t121, t126, t130, t145, t149, t161, t172, t173, 
	    t174, t190, t192, t193, t199, t204, t206, t226, t244, t245, t246, 
	    t248, t250, t254, t287, t290, t293, t297, t309, t322, t327, t338, 
	    t341, t359, t373, t387, t39, t53, t129, t133, t134, t136, t147, 
	    t148, t152, t157, t164, t165, t179, t182, t183, t185, t186, rhoa, 
	    rhob, t195, t196, t200, t202, t203, t212, t213, t216, t217, t221, 
	    t224, t227, t228, t232, t233, t234, t237, t238, t243, t247, t253, 
	    t258, t260, t264, t269, t270, t275, t276, t282, t283, t284, t286, 
	    t289, t292, t295, t299, t301, t304, t311, t317, t318, t324, t332, 
	    t333, t336, t340, t344, t345, t348, t349, t350, t356, t360, t364, 
	    t365, t366, t369, t375, t376, t382, t388, t391, t394, t395, t398, 
	    t399, t400, t401, t404, t407, t412, t417, t418, t421, t427, t428, 
	    t429, t432, t433, t437, t439, t440, t441, t443, t448, t454, t457, 
	    t460, t461, t465, t467, t468, t471, t472, t473, t476, t477, sigma,
	     t478, t480, t482, t483, t484, t485, t487, t488, t489, t491, t492,
	     t498, t504, t506, t507, t508, t510, t511, t514, t515, t518, t522,
	     t525, t526, t528, t532, t533, t535, t536, t539, t540, t541, t542,
	     t543, t546, t547, t548, t551, t552, t555, t557, t558, t560, t561,
	     t563, t565, t566, t570, t571, t574, t575, t578, t579, t582, t585,
	     t589, t590, t593, t599, t602, t606, t607, t610, t613, t617, t618,
	     t621, t624, t625, t628, t635, t636, t639, t641, t642, t645, t648,
	     t652, t653, t656, t657, t661, t665, t666, t669, t670, t671, t677,
	     t680, t683, t684, t687, t688, t691, t692, t693, t695, t696, t697,
	     t706, t708, t721, t722, t723, t726, t729, t738, t739, t743, t746,
	     t747, t750, t753, t762, t766, t780, t786, t787, t795, t799, t801,
	     t802, t805, t809, t819, t823, t827, t828, t829, t831, t836, t837,
	     t838, t839, t840, t841, t842, t843, t844, t848, t849, t850, t851,
	     t858, t860, t861, t864, t866, t869, t873, t875, t879, t881, t885,
	     t886, t887, t890, t893, t894, t898, t900, t901, t904, t906, t911,
	     t914, t917, t918, t919, t920, t923, t927, t939, t955, t958, t969,
	     t972, t983, t986, t987, t989, t992, t994, t1006, t1008, t1012, 
	    sigmaaa, sigmaab, sigmabb, t1017, t1023, t1024, t1029, t1031, 
	    t1034, t1038, t1054, t1073, t1077, t1078, t1081, t1082, t1085, 
	    t1089, t1090, t1091, t1101, t1105, t1108, t1111, t1114, t1117, 
	    t1120, t1129, t1150, t1165, t1169, t1177, t1179, t1180, t1184, 
	    t1185, t1186, t1187, t1188, t1191, t1192, t1198, t1206, t1207, 
	    t1210, t1211, t1214, t1218, t1223, t1227, t1230, t1233, t1239, 
	    t1243, t1247, t1264, t1272, t1276, t1285, t1286, t1290, t1291, 
	    t1292, t1295, t1296, t1300, t1301, t1304, t1305, t1308, t1309, 
	    t1311, t1312, t1313, t1318, t1347, t1350, t1353, t1354, t1365, 
	    t1369, t1370, t1373, t1376, t1379, t1382, t1386, t1389, t1438, 
	    t1442, t1446, t1449, t1450, t1455, t1458, t1459, t1460, t1461, 
	    t1464, t1465, t1466, t1467, t1468, t1473, t1496, t1499, t1511, 
	    t1515, t1518, t1519, t1522, t1525, t1528, t1561, t1571, t1575, 
	    t1581, t1584, t1586, t1587, t1593, t1594, t1604, t1606, t1607, 
	    t1609, t1620, t1621, t1629, t1642, t1646, t1648, t1649, t1672, 
	    t1674, t1676, t1704, t1711, t1713, t1714, t1729, t1732, t1765, 
	    t1774, t1785, t1787, t1789;


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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b2);
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t17 = log(32.16395899738507 / (t6 * 11.12037486309468 + 
			    t3 * 3.844746237447211 + t9 * 1.644733775567609 + 
			    t11 * .2405871291288192) + 1.);
		    t18 = (t3 * .1274696188700087 + 1.) * t17;
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhob, &c_b2);
		    t23 = 1 / t21 / t20;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t31 = .2162211495206379 / t27 * sigmabb * t23;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t35 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t49 = log(sigmabb * .2162211495206379 * t23 * (t31 + 1.) /
			     (t31 + 1. + .04675158550002605 / t33 * t35 / t38 
			    / t37) + 1.);
		    zk[i__] = rhob * (t18 * -.0310907 + t49 * 
			    .01554534543482745);
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b2);
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t17 = log(32.16395899738507 / (t6 * 11.12037486309468 + 
			    t3 * 3.844746237447211 + t9 * 1.644733775567609 + 
			    t11 * .2405871291288192) + 1.);
		    t18 = (t3 * .1274696188700087 + 1.) * t17;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhoa, &c_b2);
		    t23 = 1 / t21 / t20;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t31 = .2162211495206379 / t27 * sigmaaa * t23;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t35 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t49 = log(sigmaaa * .2162211495206379 * t23 * (t31 + 1.) /
			     (t31 + 1. + .04675158550002605 / t33 * t35 / t38 
			    / t37) + 1.);
		    zk[i__] = rhoa * (t18 * -.0310907 + t49 * 
			    .01554534543482745);
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b2);
		    t8 = pow_dd(&t4, &c_b3);
		    t11 = sqrt(t4);
/* Computing 2nd power */
		    d__1 = t5;
		    t13 = d__1 * d__1;
		    t19 = log(16.08197949869254 / (t8 * 5.98255043577108 + t5 
			    * 2.225569421150687 + t11 * .8004286349993634 + 
			    t13 * .1897004325747559) + 1.);
		    t21 = (t5 * .1325688999052018 + 1.) * .0621814 * t19;
		    t32 = log(29.60874997779344 / (t8 * 8.157414703487641 + 
			    t5 * 2.247591863577616 + t11 * .4300972471276643 
			    + t13 * .1911512595127338) + 1.);
		    t35 = rhoa - rhob * 1.;
		    t36 = t35 * t4;
		    t37 = t36 + 1.;
		    t38 = pow_dd(&t37, &c_b2);
		    t41 = 1. - t36 * 1.;
		    t42 = pow_dd(&t41, &c_b2);
		    t44 = t38 * t37 + t42 * t41 - 2.;
/* Computing 2nd power */
		    d__1 = t35;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rho;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t48 = d__1 * d__1;
		    t50 = t46 / t48;
		    t55 = (t5 * .06901399211255825 + 1.) * .037995525 * t32 * 
			    t44 * (1. - t50 * 1.);
		    t66 = log(32.16395899738507 / (t8 * 11.12037486309468 + 
			    t5 * 3.844746237447211 + t11 * 1.644733775567609 
			    + t13 * .2405871291288192) + 1.);
		    t72 = ((t5 * .1274696188700087 + 1.) * -.0310907 * t66 + 
			    t21) * 1.923661050931536 * t44 * t50;
/* Computing 2nd power */
		    d__1 = t38;
		    t73 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t42;
		    t75 = d__1 * d__1;
		    t77 = t73 * .5 + t75 * .5;
/* Computing 2nd power */
		    d__1 = t77;
		    t78 = d__1 * d__1;
		    t79 = t78 * t77;
		    t80 = 1 / t78;
		    t82 = pow_dd(&rho, &c_b2);
		    t84 = 1 / t82 / t47;
		    t89 = exp((-t21 + t55 + t72) * -32.16396844291482 / t79);
		    t90 = t89 - 1.;
		    t95 = .1362107888567592 / t90 * sigma * t80 * t84;
/* Computing 2nd power */
		    d__1 = t90;
		    t98 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = sigma;
		    t100 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t78;
		    t102 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t82;
		    t104 = d__1 * d__1;
		    t116 = log(sigma * .1362107888567592 * t80 * t84 * (t95 + 
			    1.) / (t95 + 1. + .01855337900098064 / t98 * t100 
			    / t102 / t104 / t48) + 1.);
		    zk[i__] = rho * (-t21 + t55 + t72 + t79 * 
			    .0310906908696549 * t116);
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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1274696188700087 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.16395899738507 / t13 + 1.;
		    t17 = log(t16);
		    t18 = t5 * t17;
		    t19 = t18 * .0310907;
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhob, &c_b2);
		    t23 = 1 / t21 / t20;
		    t24 = sigmabb * t23;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t28 = 1 / t27;
		    t29 = t28 * sigmabb;
		    t31 = t29 * .2162211495206379 * t23;
		    t32 = t31 + 1.;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
		    t34 = 1 / t33;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t35 = d__1 * d__1;
		    t36 = t34 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t40 = 1 / t38 / t37;
		    t41 = t36 * t40;
		    t43 = t31 + 1. + t41 * .04675158550002605;
		    t44 = 1 / t43;
		    t45 = t32 * t44;
		    t48 = t24 * .2162211495206379 * t45 + 1.;
		    t49 = log(t48);
		    t50 = t49 * .01554534543482745;
		    zk[i__] = rhob * (-t19 + t50);
		    vrhoa[i__] = 0.;
		    t54 = 1 / t21 / t20 / rhob;
		    t60 = 1 / t38 / t37 / rhob;
		    t62 = t28 * t44;
/* Computing 2nd power */
		    d__1 = t43;
		    t65 = d__1 * d__1;
		    t67 = t32 / t65;
		    t78 = 1 / t48;
		    t82 = 1 / t20;
		    t83 = 1 / t11 * t82;
		    t84 = t83 * t17;
/* Computing 2nd power */
		    d__1 = t13;
		    t86 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t89 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t89;
		    t90 = d__1 * d__1;
		    t105 = t5 / t86 * (-1.853395810515781 / t90 / t6 * t82 - 
			    t83 * 1.28158207914907 - .8223668877838045 / t9 * 
			    t82 - .1603914194192128 / t3 * t82) / t16;
		    t109 = t84 * -.08497977086918237 - t105 * 
			    64.32793688582964;
		    t114 = t34 * sigmabb;
		    vrhob[i__] = -t19 + t50 + rhob * .01554534543482745 * (
			    sigmabb * -.5045160155481551 * t54 * t45 - t35 * 
			    .1090870328333941 * t60 * t62 - t24 * 
			    .2162211495206379 * t67 * (t29 * 
			    -.5045160155481551 * t54 - t36 * 
			    .2181740656667882 * t60)) * t78 + rhob * (t84 * 
			    .001321039893133927 + t105 * 1. + (t41 * 
			    -.04675158550002605 * t109 * t26 * t44 - t24 * 
			    .2162211495206379 * t67 * (t114 * 
			    -.2162211495206379 * t23 * t109 * t26 - 
			    .0935031710000521 / t33 / t27 * t35 * t40 * t109 *
			     t26)) * .01554534543482745 * t78);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = rhob * .01554534543482745 * (t23 * 
			    .2162211495206379 * t32 * t44 + sigmabb * 
			    .04675158550002605 * t40 * t62 - t24 * 
			    .2162211495206379 * t67 * (t28 * 
			    .2162211495206379 * t23 + t114 * 
			    .0935031710000521 * t40)) * t78;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1274696188700087 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.16395899738507 / t13 + 1.;
		    t17 = log(t16);
		    t18 = t5 * t17;
		    t19 = t18 * .0310907;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhoa, &c_b2);
		    t23 = 1 / t21 / t20;
		    t24 = sigmaaa * t23;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t28 = 1 / t27;
		    t29 = t28 * sigmaaa;
		    t31 = t29 * .2162211495206379 * t23;
		    t32 = t31 + 1.;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
		    t34 = 1 / t33;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t35 = d__1 * d__1;
		    t36 = t34 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t40 = 1 / t38 / t37;
		    t41 = t36 * t40;
		    t43 = t31 + 1. + t41 * .04675158550002605;
		    t44 = 1 / t43;
		    t45 = t32 * t44;
		    t48 = t24 * .2162211495206379 * t45 + 1.;
		    t49 = log(t48);
		    t50 = t49 * .01554534543482745;
		    zk[i__] = rhoa * (-t19 + t50);
		    t54 = 1 / t21 / t20 / rhoa;
		    t60 = 1 / t38 / t37 / rhoa;
		    t62 = t28 * t44;
/* Computing 2nd power */
		    d__1 = t43;
		    t65 = d__1 * d__1;
		    t67 = t32 / t65;
		    t78 = 1 / t48;
		    t82 = 1 / t20;
		    t83 = 1 / t11 * t82;
		    t84 = t83 * t17;
/* Computing 2nd power */
		    d__1 = t13;
		    t86 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t6;
		    t89 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t89;
		    t90 = d__1 * d__1;
		    t105 = t5 / t86 * (-1.853395810515781 / t90 / t6 * t82 - 
			    t83 * 1.28158207914907 - .8223668877838045 / t9 * 
			    t82 - .1603914194192128 / t3 * t82) / t16;
		    t109 = t84 * -.08497977086918237 - t105 * 
			    64.32793688582964;
		    t114 = t34 * sigmaaa;
		    vrhoa[i__] = -t19 + t50 + rhoa * .01554534543482745 * (
			    sigmaaa * -.5045160155481551 * t54 * t45 - t35 * 
			    .1090870328333941 * t60 * t62 - t24 * 
			    .2162211495206379 * t67 * (t29 * 
			    -.5045160155481551 * t54 - t36 * 
			    .2181740656667882 * t60)) * t78 + rhoa * (t84 * 
			    .001321039893133927 + t105 * 1. + (t41 * 
			    -.04675158550002605 * t109 * t26 * t44 - t24 * 
			    .2162211495206379 * t67 * (t114 * 
			    -.2162211495206379 * t23 * t109 * t26 - 
			    .0935031710000521 / t33 / t27 * t35 * t40 * t109 *
			     t26)) * .01554534543482745 * t78);
		    vrhob[i__] = 0.;
		    vsigmaaa[i__] = rhoa * .01554534543482745 * (t23 * 
			    .2162211495206379 * t32 * t44 + sigmaaa * 
			    .04675158550002605 * t40 * t62 - t24 * 
			    .2162211495206379 * t67 * (t28 * 
			    .2162211495206379 * t23 + t114 * 
			    .0935031710000521 * t40)) * t78;
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b2);
		    t7 = t5 * .1325688999052018 + 1.;
		    t8 = pow_dd(&t4, &c_b3);
		    t11 = sqrt(t4);
/* Computing 2nd power */
		    d__1 = t5;
		    t13 = d__1 * d__1;
		    t15 = t8 * 5.98255043577108 + t5 * 2.225569421150687 + 
			    t11 * .8004286349993634 + t13 * .1897004325747559;
		    t18 = 16.08197949869254 / t15 + 1.;
		    t19 = log(t18);
		    t21 = t7 * .0621814 * t19;
		    t23 = t5 * .06901399211255825 + 1.;
		    t28 = t8 * 8.157414703487641 + t5 * 2.247591863577616 + 
			    t11 * .4300972471276643 + t13 * .1911512595127338;
		    t31 = 29.60874997779344 / t28 + 1.;
		    t32 = log(t31);
		    t33 = t23 * t32;
		    t35 = rhoa - rhob * 1.;
		    t36 = t35 * t4;
		    t37 = t36 + 1.;
		    t38 = pow_dd(&t37, &c_b2);
		    t41 = 1. - t36 * 1.;
		    t42 = pow_dd(&t41, &c_b2);
		    t44 = t38 * t37 + t42 * t41 - 2.;
/* Computing 2nd power */
		    d__1 = t35;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rho;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t48 = d__1 * d__1;
		    t49 = 1 / t48;
		    t50 = t46 * t49;
		    t52 = 1. - t50 * 1.;
		    t55 = t33 * .037995525 * t44 * t52;
		    t57 = t5 * .1274696188700087 + 1.;
		    t62 = t8 * 11.12037486309468 + t5 * 3.844746237447211 + 
			    t11 * 1.644733775567609 + t13 * .2405871291288192;
		    t65 = 32.16395899738507 / t62 + 1.;
		    t66 = log(t65);
		    t69 = t57 * -.0310907 * t66 + t21;
		    t70 = t69 * t44;
		    t72 = t70 * 1.923661050931536 * t50;
/* Computing 2nd power */
		    d__1 = t38;
		    t73 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t42;
		    t75 = d__1 * d__1;
		    t77 = t73 * .5 + t75 * .5;
/* Computing 2nd power */
		    d__1 = t77;
		    t78 = d__1 * d__1;
		    t79 = t78 * t77;
		    t80 = 1 / t78;
		    t81 = sigma * t80;
		    t82 = pow_dd(&rho, &c_b2);
		    t84 = 1 / t82 / t47;
		    t85 = -t21 + t55 + t72;
		    t86 = 1 / t79;
		    t89 = exp(t85 * -32.16396844291482 * t86);
		    t90 = t89 - 1.;
		    t91 = 1 / t90;
		    t92 = t91 * sigma;
		    t93 = t80 * t84;
		    t95 = t92 * .1362107888567592 * t93;
		    t96 = t95 + 1.;
/* Computing 2nd power */
		    d__1 = t90;
		    t98 = d__1 * d__1;
		    t99 = 1 / t98;
/* Computing 2nd power */
		    d__1 = sigma;
		    t100 = d__1 * d__1;
		    t101 = t99 * t100;
/* Computing 2nd power */
		    d__1 = t78;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
/* Computing 2nd power */
		    d__1 = t82;
		    t104 = d__1 * d__1;
		    t106 = 1 / t104 / t48;
		    t107 = t103 * t106;
		    t110 = t95 + 1. + t101 * .01855337900098064 * t107;
		    t111 = 1 / t110;
		    t115 = t81 * .1362107888567592 * t84 * t96 * t111 + 1.;
		    t116 = log(t115);
		    t118 = t79 * .0310906908696549 * t116;
		    zk[i__] = rho * (-t21 + t55 + t72 + t118);
		    t124 = t38 * 1.333333333333333 * t4 - t42 * 
			    1.333333333333333 * t4;
		    t127 = t33 * .037995525 * t124 * t52;
		    t128 = t45 * t35;
		    t131 = t33 * t44 * t128 * t49;
		    t132 = t131 * .1519821;
		    t135 = t69 * 1.923661050931536 * t124 * t50;
		    t137 = t70 * t128 * t49;
		    t138 = t137 * 7.694644203726145;
		    t139 = t78 * t116;
		    t140 = 1 / t38;
		    t143 = 1 / t42;
		    t146 = t140 * .3333333333333333 * t4 - t143 * 
			    .3333333333333333 * t4;
		    t150 = sigma * t86 * t84;
		    t151 = t96 * t111;
		    t155 = t99 * sigma;
		    t156 = t155 * t80;
		    t160 = t85 * t103;
		    t163 = (t127 - t132 + t135 + t138) * -32.16396844291482 * 
			    t86 + t160 * 96.49190532874446 * t146;
		    t167 = t156 * .1362107888567592 * t84 * t163 * t89;
		    t168 = t86 * t84;
		    t171 = t92 * .2724215777135184 * t168 * t146;
		    t177 = t81 * t84;
/* Computing 2nd power */
		    d__1 = t110;
		    t178 = d__1 * d__1;
		    t180 = t96 / t178;
		    t184 = 1 / t98 / t90 * t100 * t103;
		    t191 = 1 / t102 / t77 * t106;
		    t201 = 1 / t115;
		    t207 = 1 / t47;
		    t208 = 1 / t13 * t207;
		    t209 = t208 * t19;
		    t210 = t209 * .002747773264188438;
/* Computing 2nd power */
		    d__1 = t15;
		    t211 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t8;
		    t214 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t214;
		    t215 = d__1 * d__1;
		    t218 = 1 / t215 / t8 * t207;
		    t222 = 1 / t11 * t207;
		    t225 = 1 / t5 * t207;
		    t230 = t7 / t211 * (t218 * -.99709173929518 - t208 * 
			    .7418564737168958 - t222 * .4002143174996817 - 
			    t225 * .1264669550498372) / t18;
		    t231 = t230 * 1.;
		    t235 = t208 * 8.7407428755417e-4 * t32 * t44 * t52;
/* Computing 2nd power */
		    d__1 = t28;
		    t236 = d__1 * d__1;
		    t249 = t23 * 1.125 / t236 * (t218 * -1.35956911724794 - 
			    t208 * .7491972878592054 - t222 * 
			    .2150486235638321 - t225 * .1274341730084892) / 
			    t31 * t44 * t52;
		    t256 = t38 * -1.333333333333333 * t35 * t207 + t42 * 
			    1.333333333333333 * t35 * t207;
		    t259 = t33 * .037995525 * t256 * t52;
		    t261 = t48 * rho;
		    t262 = 1 / t261;
		    t265 = t33 * .1519821 * t44 * t46 * t262;
/* Computing 2nd power */
		    d__1 = t62;
		    t268 = d__1 * d__1;
		    t285 = (t208 * .001321039893133927 * t66 + t57 * 1. / 
			    t268 * (t218 * -1.853395810515781 - t208 * 
			    1.28158207914907 - t222 * .8223668877838045 - 
			    t225 * .1603914194192128) / t65 - t209 * 
			    .002747773264188438 - t230 * 1.) * 
			    1.923661050931536 * t44 * t50;
		    t288 = t69 * 1.923661050931536 * t256 * t50;
		    t291 = t70 * 7.694644203726145 * t46 * t262;
		    t298 = t140 * -.3333333333333333 * t35 * t207 + t143 * 
			    .3333333333333333 * t35 * t207;
		    t306 = 1 / t82 / t47 / rho;
		    t316 = (t210 + t231 - t235 - t249 + t259 + t265 + t285 + 
			    t288 - t291) * -32.16396844291482 * t86 + t160 * 
			    96.49190532874446 * t298;
		    t320 = t156 * .1362107888567592 * t84 * t316 * t89;
		    t323 = t92 * .2724215777135184 * t168 * t298;
		    t326 = t92 * .3178251739991049 * t80 * t306;
		    t352 = t210 + t231 - t235 - t249 + t259 + t265 + t285 + 
			    t288 - t291 + t139 * .09327207260896469 * t298 + 
			    t79 * .0310906908696549 * (t150 * 
			    -.2724215777135184 * t151 * t298 - t81 * 
			    .3178251739991049 * t306 * t96 * t111 + t81 * 
			    .1362107888567592 * t84 * (-t320 - t323 - t326) * 
			    t111 - t177 * .1362107888567592 * t180 * (-t320 - 
			    t323 - t326 - t184 * .03710675800196129 * t106 * 
			    t316 * t89 - t101 * .07421351600392257 * t191 * 
			    t298 - t101 * .08658243533790967 * t103 / t104 / 
			    t261)) * t201;
		    t353 = rho * t352;
		    vrhoa[i__] = rho * (t127 - t132 + t135 + t138 + t139 * 
			    .09327207260896469 * t146 + t79 * 
			    .0310906908696549 * (t150 * -.2724215777135184 * 
			    t151 * t146 + t81 * .1362107888567592 * t84 * (
			    -t167 - t171) * t111 - t177 * .1362107888567592 * 
			    t180 * (-t167 - t171 - t184 * .03710675800196129 *
			     t106 * t163 * t89 - t101 * .07421351600392257 * 
			    t191 * t146)) * t201) - t21 + t55 + t72 + t118 + 
			    t353;
		    t354 = -t124;
		    t357 = t33 * .037995525 * t354 * t52;
		    t358 = t131 * .1519821;
		    t361 = t69 * 1.923661050931536 * t354 * t50;
		    t362 = t137 * 7.694644203726145;
		    t363 = -t146;
		    t374 = (t357 + t358 + t361 - t362) * -32.16396844291482 * 
			    t86 + t160 * 96.49190532874446 * t363;
		    t378 = t156 * .1362107888567592 * t84 * t374 * t89;
		    t381 = t92 * .2724215777135184 * t168 * t363;
		    vrhob[i__] = rho * (t357 + t358 + t361 - t362 + t139 * 
			    .09327207260896469 * t363 + t79 * 
			    .0310906908696549 * (t150 * -.2724215777135184 * 
			    t151 * t363 + t81 * .1362107888567592 * t84 * (
			    -t378 - t381) * t111 - t177 * .1362107888567592 * 
			    t180 * (-t378 - t381 - t184 * .03710675800196129 *
			     t106 * t374 * t89 - t101 * .07421351600392257 * 
			    t191 * t363)) * t201) - t21 + t55 + t72 + t118 + 
			    t353;
		    t423 = rho * t79 * (t93 * .1362107888567592 * t151 + 
			    sigma * .01855337900098064 * t103 * t106 * t91 * 
			    t111 - t177 * .1362107888567592 * t180 * (t91 * 
			    .1362107888567592 * t80 * t84 + t155 * 
			    .03710675800196129 * t107)) * t201;
		    vsigmaaa[i__] = t423 * .0310906908696549;
		    vsigmaab[i__] = t423 * .06218138173930979;
		    vsigmabb[i__] = vsigmaaa[i__];
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
		    t2 = 1 / rhob;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1274696188700087 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.16395899738507 / t13 + 1.;
		    t17 = log(t16);
		    t18 = t5 * t17;
		    t19 = t18 * .0310907;
/* Computing 2nd power */
		    d__1 = rhob;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhob, &c_b2);
		    t23 = 1 / t21 / t20;
		    t24 = sigmabb * t23;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t28 = 1 / t27;
		    t29 = t28 * sigmabb;
		    t31 = t29 * .2162211495206379 * t23;
		    t32 = t31 + 1.;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
		    t34 = 1 / t33;
/* Computing 2nd power */
		    d__1 = sigmabb;
		    t35 = d__1 * d__1;
		    t36 = t34 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t40 = 1 / t38 / t37;
		    t41 = t36 * t40;
		    t43 = t31 + 1. + t41 * .04675158550002605;
		    t44 = 1 / t43;
		    t45 = t32 * t44;
		    t48 = t24 * .2162211495206379 * t45 + 1.;
		    t49 = log(t48);
		    t50 = t49 * .01554534543482745;
		    zk[i__] = rhob * (-t19 + t50);
		    vrhoa[i__] = 0.;
		    t52 = t20 * rhob;
		    t54 = 1 / t21 / t52;
		    t55 = sigmabb * t54;
		    t60 = 1 / t38 / t37 / rhob;
		    t61 = t35 * t60;
		    t62 = t28 * t44;
/* Computing 2nd power */
		    d__1 = t43;
		    t65 = d__1 * d__1;
		    t66 = 1 / t65;
		    t67 = t32 * t66;
		    t70 = t36 * t60;
		    t72 = t29 * -.5045160155481551 * t54 - t70 * 
			    .2181740656667882;
		    t73 = t67 * t72;
		    t76 = t55 * -.5045160155481551 * t45 - t61 * 
			    .1090870328333941 * t62 - t24 * .2162211495206379 
			    * t73;
		    t77 = rhob * t76;
		    t78 = 1 / t48;
		    t81 = 1 / t11;
		    t82 = 1 / t20;
		    t83 = t81 * t82;
		    t84 = t83 * t17;
/* Computing 2nd power */
		    d__1 = t13;
		    t86 = d__1 * d__1;
		    t87 = 1 / t86;
		    t88 = t5 * t87;
/* Computing 2nd power */
		    d__1 = t6;
		    t89 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t89;
		    t90 = d__1 * d__1;
		    t91 = t90 * t6;
		    t92 = 1 / t91;
		    t96 = 1 / t9;
		    t99 = 1 / t3;
		    t102 = t92 * -1.853395810515781 * t82 - t83 * 
			    1.28158207914907 - t96 * .8223668877838045 * t82 
			    - t99 * .1603914194192128 * t82;
		    t103 = 1 / t16;
		    t105 = t88 * t102 * t103;
		    t109 = t84 * -.08497977086918237 - t105 * 
			    64.32793688582964;
		    t110 = t109 * t26;
		    t111 = t110 * t44;
		    t114 = t34 * sigmabb;
		    t120 = 1 / t33 / t27;
		    t121 = t120 * t35;
		    t126 = t114 * -.2162211495206379 * t23 * t109 * t26 - 
			    t121 * .0935031710000521 * t40 * t109 * t26;
		    t127 = t67 * t126;
		    t130 = t41 * -.04675158550002605 * t111 - t24 * 
			    .2162211495206379 * t127;
		    t131 = t130 * t78;
		    vrhob[i__] = -t19 + t50 + t77 * .01554534543482745 * t78 
			    + rhob * (t84 * .001321039893133927 + t105 * 1. + 
			    t131 * .01554534543482745);
		    vsigmaaa[i__] = 0.;
		    vsigmaab[i__] = 0.;
		    t135 = t23 * t32;
		    t138 = sigmabb * t40;
		    t145 = t28 * .2162211495206379 * t23 + t114 * 
			    .0935031710000521 * t40;
		    t149 = t135 * .2162211495206379 * t44 + t138 * 
			    .04675158550002605 * t62 - t24 * 
			    .2162211495206379 * t67 * t145;
		    vsigmabb[i__] = rhob * .01554534543482745 * t149 * t78;
		    v2rhoa2[i__] = 0.;
		    v2rhoab[i__] = 0.;
		    t155 = 1 / t21 / t37;
		    t161 = 1 / t38 / t37 / t20;
		    t167 = t28 * t66;
		    t172 = 1 / t65 / t43;
		    t173 = t32 * t172;
/* Computing 2nd power */
		    d__1 = t72;
		    t174 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t76;
		    t190 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t48;
		    t192 = d__1 * d__1;
		    t193 = 1 / t192;
		    t199 = t70 * t111;
		    t204 = t41 * .04675158550002605 * t110 * t66 * t72;
		    t206 = t55 * .5045160155481551 * t127;
		    t209 = t61 * .1090870328333941 * t167 * t126;
		    t214 = t24 * .4324422990412758 * t32 * t172 * t126 * t72;
		    t226 = t24 * .2162211495206379 * t67 * (t114 * 
			    .5045160155481551 * t54 * t109 * t26 + t121 * 
			    .4363481313335765 * t60 * t109 * t26);
		    t230 = t130 * t193;
		    t244 = 1 / t37;
		    t245 = 1 / t11 / t2 * t244;
		    t246 = t245 * t17;
		    t248 = 1 / t52;
		    t249 = t81 * t248;
		    t250 = t249 * t17;
		    t254 = t83 * t87 * t102 * t103;
/* Computing 2nd power */
		    d__1 = t102;
		    t259 = d__1 * d__1;
		    t261 = t5 / t86 / t13 * t259 * t103;
		    t285 = t88 * (-1.544496508763151 / t91 / t2 * t244 + t92 *
			     3.706791621031562 * t248 - t245 * 
			    .854388052766047 + t249 * 2.563164158298141 - 
			    .4111834438919023 / t9 / t2 * t244 + t96 * 
			    1.644733775567609 * t248 - .05346380647307093 / 
			    t3 / t2 * t244 + t99 * .3207828388384256 * t248) *
			     t103;
/* Computing 2nd power */
		    d__1 = t86;
		    t287 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t16;
		    t290 = d__1 * d__1;
		    t293 = t5 / t287 * t259 / t290;
/* Computing 2nd power */
		    d__1 = t109;
		    t297 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t26;
		    t298 = d__1 * d__1;
		    t309 = t246 * -.05665318057945491 + t250 * 
			    .1699595417383647 + t254 * 5.46657173168712 + 
			    t261 * 128.6558737716593 - t285 * 
			    64.32793688582964 - t293 * 2069.041124382199;
/* Computing 2nd power */
		    d__1 = t126;
		    t322 = d__1 * d__1;
		    t327 = t23 * t297;
/* Computing 2nd power */
		    d__1 = t33;
		    t338 = d__1 * d__1;
		    t341 = t40 * t297;
/* Computing 2nd power */
		    d__1 = t130;
		    t359 = d__1 * d__1;
		    s1 = t76 * .0310906908696549 * t78 + rhob * 
			    .01554534543482745 * (sigmabb * 1.681720051827184 
			    * t155 * t45 + t35 * .8726962626671529 * t161 * 
			    t62 + t55 * 1.00903203109631 * t73 + t61 * 
			    .2181740656667882 * t167 * t72 + t24 * 
			    .4324422990412758 * t173 * t174 - t24 * 
			    .2162211495206379 * t67 * (t29 * 
			    1.681720051827184 * t155 + t36 * 
			    1.236319705445133 * t161)) * t78 - rhob * 
			    .01554534543482745 * t190 * t193 + t84 * 
			    .002642079786267853 + t105 * 2.;
		    s2 = s1 + t131 * .0310906908696549 + rhob * ((t199 * 
			    .2181740656667882 + t204 + t206 + t209 + t214 - 
			    t226) * .01554534543482745 * t78 - t230 * 
			    .01554534543482745 * t76);
		    s3 = s2 + rhob * .01554534543482745 * (t199 * 
			    .2181740656667882 + t206 + t209 + t204 + t214 - 
			    t226) * t78;
		    v2rhob2[i__] = s3 - t77 * .01554534543482745 * t230 + 
			    rhob * (t246 * 8.806932620892844e-4 - t250 * 
			    .002642079786267853 - t254 * .08497974591333914 - 
			    t261 * 2. + t285 * 1. + t293 * 32.16395899738507 
			    + (t35 * .0935031710000521 * t40 * t120 * t297 * 
			    t298 * t44 - t41 * .04675158550002605 * t309 * 
			    t26 * t44 - t41 * .04675158550002605 * t297 * t26 
			    * t44 + t41 * .0935031710000521 * t110 * t66 * 
			    t126 + t24 * .4324422990412758 * t173 * t322 - 
			    t24 * .2162211495206379 * t67 * (t120 * 
			    .4324422990412758 * sigmabb * t327 * t298 - t114 *
			     .2162211495206379 * t23 * t309 * t26 - t114 * 
			    .2162211495206379 * t327 * t26 + 
			    .2805095130001563 / t338 * t35 * t341 * t298 - 
			    t121 * .0935031710000521 * t40 * t309 * t26 - 
			    t121 * .0935031710000521 * t341 * t26)) * 
			    .01554534543482745 * t78 - t359 * 
			    .01554534543482745 * t193);
		    v2sigmaaa2[i__] = 0.;
		    v2sigmaaaab[i__] = 0.;
		    v2sigmaaabb[i__] = 0.;
		    v2sigmaab2[i__] = 0.;
		    v2sigmaabbb[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t145;
		    t373 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t149;
		    t387 = d__1 * d__1;
		    v2sigmabb2[i__] = rhob * .01554534543482745 * (t40 * 
			    .0935031710000521 * t28 * t44 - t135 * 
			    .4324422990412758 * t66 * t145 - t138 * 
			    .0935031710000521 * t167 * t145 + t24 * 
			    .4324422990412758 * t173 * t373 - sigmabb * 
			    .02021736311745604 / t37 / t52 * t67 * t34) * t78 
			    - rhob * .01554534543482745 * t387 * t193;
		} else if (rhob < 1e-20) {
		    rho = rhoa;
/* Computing MAX */
		    d__1 = 0., d__2 = sigmaaa1[i__];
		    sigmaaa = max(d__1,d__2);
		    sigma = sigmaaa;
		    t2 = 1 / rhoa;
		    t3 = pow_dd(&t2, &c_b2);
		    t5 = t3 * .1274696188700087 + 1.;
		    t6 = pow_dd(&t2, &c_b3);
		    t9 = sqrt(t2);
/* Computing 2nd power */
		    d__1 = t3;
		    t11 = d__1 * d__1;
		    t13 = t6 * 11.12037486309468 + t3 * 3.844746237447211 + 
			    t9 * 1.644733775567609 + t11 * .2405871291288192;
		    t16 = 32.16395899738507 / t13 + 1.;
		    t17 = log(t16);
		    t18 = t5 * t17;
		    t19 = t18 * .0310907;
/* Computing 2nd power */
		    d__1 = rhoa;
		    t20 = d__1 * d__1;
		    t21 = pow_dd(&rhoa, &c_b2);
		    t23 = 1 / t21 / t20;
		    t24 = sigmaaa * t23;
		    t26 = exp(t18 * 2.000000587336264);
		    t27 = t26 - 1.;
		    t28 = 1 / t27;
		    t29 = t28 * sigmaaa;
		    t31 = t29 * .2162211495206379 * t23;
		    t32 = t31 + 1.;
/* Computing 2nd power */
		    d__1 = t27;
		    t33 = d__1 * d__1;
		    t34 = 1 / t33;
/* Computing 2nd power */
		    d__1 = sigmaaa;
		    t35 = d__1 * d__1;
		    t36 = t34 * t35;
/* Computing 2nd power */
		    d__1 = t20;
		    t37 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t21;
		    t38 = d__1 * d__1;
		    t40 = 1 / t38 / t37;
		    t41 = t36 * t40;
		    t43 = t31 + 1. + t41 * .04675158550002605;
		    t44 = 1 / t43;
		    t45 = t32 * t44;
		    t48 = t24 * .2162211495206379 * t45 + 1.;
		    t49 = log(t48);
		    t50 = t49 * .01554534543482745;
		    zk[i__] = rhoa * (-t19 + t50);
		    t52 = t20 * rhoa;
		    t54 = 1 / t21 / t52;
		    t55 = sigmaaa * t54;
		    t60 = 1 / t38 / t37 / rhoa;
		    t61 = t35 * t60;
		    t62 = t28 * t44;
/* Computing 2nd power */
		    d__1 = t43;
		    t65 = d__1 * d__1;
		    t66 = 1 / t65;
		    t67 = t32 * t66;
		    t70 = t36 * t60;
		    t72 = t29 * -.5045160155481551 * t54 - t70 * 
			    .2181740656667882;
		    t73 = t67 * t72;
		    t76 = t55 * -.5045160155481551 * t45 - t61 * 
			    .1090870328333941 * t62 - t24 * .2162211495206379 
			    * t73;
		    t77 = rhoa * t76;
		    t78 = 1 / t48;
		    t81 = 1 / t11;
		    t82 = 1 / t20;
		    t83 = t81 * t82;
		    t84 = t83 * t17;
/* Computing 2nd power */
		    d__1 = t13;
		    t86 = d__1 * d__1;
		    t87 = 1 / t86;
		    t88 = t5 * t87;
/* Computing 2nd power */
		    d__1 = t6;
		    t89 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t89;
		    t90 = d__1 * d__1;
		    t91 = t90 * t6;
		    t92 = 1 / t91;
		    t96 = 1 / t9;
		    t99 = 1 / t3;
		    t102 = t92 * -1.853395810515781 * t82 - t83 * 
			    1.28158207914907 - t96 * .8223668877838045 * t82 
			    - t99 * .1603914194192128 * t82;
		    t103 = 1 / t16;
		    t105 = t88 * t102 * t103;
		    t109 = t84 * -.08497977086918237 - t105 * 
			    64.32793688582964;
		    t110 = t109 * t26;
		    t111 = t110 * t44;
		    t114 = t34 * sigmaaa;
		    t120 = 1 / t33 / t27;
		    t121 = t120 * t35;
		    t126 = t114 * -.2162211495206379 * t23 * t109 * t26 - 
			    t121 * .0935031710000521 * t40 * t109 * t26;
		    t127 = t67 * t126;
		    t130 = t41 * -.04675158550002605 * t111 - t24 * 
			    .2162211495206379 * t127;
		    t131 = t130 * t78;
		    vrhoa[i__] = -t19 + t50 + t77 * .01554534543482745 * t78 
			    + rhoa * (t84 * .001321039893133927 + t105 * 1. + 
			    t131 * .01554534543482745);
		    vrhob[i__] = 0.;
		    t135 = t23 * t32;
		    t138 = sigmaaa * t40;
		    t145 = t28 * .2162211495206379 * t23 + t114 * 
			    .0935031710000521 * t40;
		    t149 = t135 * .2162211495206379 * t44 + t138 * 
			    .04675158550002605 * t62 - t24 * 
			    .2162211495206379 * t67 * t145;
		    vsigmaaa[i__] = rhoa * .01554534543482745 * t149 * t78;
		    vsigmaab[i__] = 0.;
		    vsigmabb[i__] = 0.;
		    t155 = 1 / t21 / t37;
		    t161 = 1 / t38 / t37 / t20;
		    t167 = t28 * t66;
		    t172 = 1 / t65 / t43;
		    t173 = t32 * t172;
/* Computing 2nd power */
		    d__1 = t72;
		    t174 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t76;
		    t190 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t48;
		    t192 = d__1 * d__1;
		    t193 = 1 / t192;
		    t199 = t70 * t111;
		    t204 = t41 * .04675158550002605 * t110 * t66 * t72;
		    t206 = t55 * .5045160155481551 * t127;
		    t209 = t61 * .1090870328333941 * t167 * t126;
		    t214 = t24 * .4324422990412758 * t32 * t172 * t126 * t72;
		    t226 = t24 * .2162211495206379 * t67 * (t114 * 
			    .5045160155481551 * t54 * t109 * t26 + t121 * 
			    .4363481313335765 * t60 * t109 * t26);
		    t230 = t130 * t193;
		    t244 = 1 / t37;
		    t245 = 1 / t11 / t2 * t244;
		    t246 = t245 * t17;
		    t248 = 1 / t52;
		    t249 = t81 * t248;
		    t250 = t249 * t17;
		    t254 = t83 * t87 * t102 * t103;
/* Computing 2nd power */
		    d__1 = t102;
		    t259 = d__1 * d__1;
		    t261 = t5 / t86 / t13 * t259 * t103;
		    t285 = t88 * (-1.544496508763151 / t91 / t2 * t244 + t92 *
			     3.706791621031562 * t248 - t245 * 
			    .854388052766047 + t249 * 2.563164158298141 - 
			    .4111834438919023 / t9 / t2 * t244 + t96 * 
			    1.644733775567609 * t248 - .05346380647307093 / 
			    t3 / t2 * t244 + t99 * .3207828388384256 * t248) *
			     t103;
/* Computing 2nd power */
		    d__1 = t86;
		    t287 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t16;
		    t290 = d__1 * d__1;
		    t293 = t5 / t287 * t259 / t290;
/* Computing 2nd power */
		    d__1 = t109;
		    t297 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t26;
		    t298 = d__1 * d__1;
		    t309 = t246 * -.05665318057945491 + t250 * 
			    .1699595417383647 + t254 * 5.46657173168712 + 
			    t261 * 128.6558737716593 - t285 * 
			    64.32793688582964 - t293 * 2069.041124382199;
/* Computing 2nd power */
		    d__1 = t126;
		    t322 = d__1 * d__1;
		    t327 = t23 * t297;
/* Computing 2nd power */
		    d__1 = t33;
		    t338 = d__1 * d__1;
		    t341 = t40 * t297;
/* Computing 2nd power */
		    d__1 = t130;
		    t359 = d__1 * d__1;
		    s1 = t76 * .0310906908696549 * t78 + rhoa * 
			    .01554534543482745 * (sigmaaa * 1.681720051827184 
			    * t155 * t45 + t35 * .8726962626671529 * t161 * 
			    t62 + t55 * 1.00903203109631 * t73 + t61 * 
			    .2181740656667882 * t167 * t72 + t24 * 
			    .4324422990412758 * t173 * t174 - t24 * 
			    .2162211495206379 * t67 * (t29 * 
			    1.681720051827184 * t155 + t36 * 
			    1.236319705445133 * t161)) * t78 - rhoa * 
			    .01554534543482745 * t190 * t193 + t84 * 
			    .002642079786267853 + t105 * 2.;
		    s2 = s1 + t131 * .0310906908696549 + rhoa * ((t199 * 
			    .2181740656667882 + t204 + t206 + t209 + t214 - 
			    t226) * .01554534543482745 * t78 - t230 * 
			    .01554534543482745 * t76);
		    s3 = s2 + rhoa * .01554534543482745 * (t199 * 
			    .2181740656667882 + t206 + t209 + t204 + t214 - 
			    t226) * t78;
		    v2rhoa2[i__] = s3 - t77 * .01554534543482745 * t230 + 
			    rhoa * (t246 * 8.806932620892844e-4 - t250 * 
			    .002642079786267853 - t254 * .08497974591333914 - 
			    t261 * 2. + t285 * 1. + t293 * 32.16395899738507 
			    + (t35 * .0935031710000521 * t40 * t120 * t297 * 
			    t298 * t44 - t41 * .04675158550002605 * t309 * 
			    t26 * t44 - t41 * .04675158550002605 * t297 * t26 
			    * t44 + t41 * .0935031710000521 * t110 * t66 * 
			    t126 + t24 * .4324422990412758 * t173 * t322 - 
			    t24 * .2162211495206379 * t67 * (t120 * 
			    .4324422990412758 * sigmaaa * t327 * t298 - t114 *
			     .2162211495206379 * t23 * t309 * t26 - t114 * 
			    .2162211495206379 * t327 * t26 + 
			    .2805095130001563 / t338 * t35 * t341 * t298 - 
			    t121 * .0935031710000521 * t40 * t309 * t26 - 
			    t121 * .0935031710000521 * t341 * t26)) * 
			    .01554534543482745 * t78 - t359 * 
			    .01554534543482745 * t193);
		    v2rhob2[i__] = 0.;
		    v2rhoab[i__] = 0.;
/* Computing 2nd power */
		    d__1 = t145;
		    t373 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t149;
		    t387 = d__1 * d__1;
		    v2sigmaaa2[i__] = rhoa * .01554534543482745 * (t40 * 
			    .0935031710000521 * t28 * t44 - t135 * 
			    .4324422990412758 * t66 * t145 - t138 * 
			    .0935031710000521 * t167 * t145 + t24 * 
			    .4324422990412758 * t173 * t373 - sigmaaa * 
			    .02021736311745604 / t37 / t52 * t67 * t34) * t78 
			    - rhoa * .01554534543482745 * t387 * t193;
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
		    t4 = 1 / rho;
		    t5 = pow_dd(&t4, &c_b2);
		    t7 = t5 * .1325688999052018 + 1.;
		    t8 = pow_dd(&t4, &c_b3);
		    t11 = sqrt(t4);
/* Computing 2nd power */
		    d__1 = t5;
		    t13 = d__1 * d__1;
		    t15 = t8 * 5.98255043577108 + t5 * 2.225569421150687 + 
			    t11 * .8004286349993634 + t13 * .1897004325747559;
		    t18 = 16.08197949869254 / t15 + 1.;
		    t19 = log(t18);
		    t21 = t7 * .0621814 * t19;
		    t23 = t5 * .06901399211255825 + 1.;
		    t28 = t8 * 8.157414703487641 + t5 * 2.247591863577616 + 
			    t11 * .4300972471276643 + t13 * .1911512595127338;
		    t31 = 29.60874997779344 / t28 + 1.;
		    t32 = log(t31);
		    t33 = t23 * t32;
		    t35 = rhoa - rhob * 1.;
		    t36 = t35 * t4;
		    t37 = t36 + 1.;
		    t38 = pow_dd(&t37, &c_b2);
		    t39 = t38 * t37;
		    t41 = 1. - t36 * 1.;
		    t42 = pow_dd(&t41, &c_b2);
		    t43 = t42 * t41;
		    t44 = t39 + t43 - 2.;
/* Computing 2nd power */
		    d__1 = t35;
		    t45 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t45;
		    t46 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = rho;
		    t47 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t47;
		    t48 = d__1 * d__1;
		    t49 = 1 / t48;
		    t50 = t46 * t49;
		    t52 = 1. - t50 * 1.;
		    t53 = t44 * t52;
		    t55 = t33 * .037995525 * t53;
		    t57 = t5 * .1274696188700087 + 1.;
		    t62 = t8 * 11.12037486309468 + t5 * 3.844746237447211 + 
			    t11 * 1.644733775567609 + t13 * .2405871291288192;
		    t65 = 32.16395899738507 / t62 + 1.;
		    t66 = log(t65);
		    t69 = t57 * -.0310907 * t66 + t21;
		    t70 = t69 * t44;
		    t72 = t70 * 1.923661050931536 * t50;
/* Computing 2nd power */
		    d__1 = t38;
		    t73 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t42;
		    t75 = d__1 * d__1;
		    t77 = t73 * .5 + t75 * .5;
/* Computing 2nd power */
		    d__1 = t77;
		    t78 = d__1 * d__1;
		    t79 = t78 * t77;
		    t80 = 1 / t78;
		    t81 = sigma * t80;
		    t82 = pow_dd(&rho, &c_b2);
		    t84 = 1 / t82 / t47;
		    t85 = -t21 + t55 + t72;
		    t86 = 1 / t79;
		    t89 = exp(t85 * -32.16396844291482 * t86);
		    t90 = t89 - 1.;
		    t91 = 1 / t90;
		    t92 = t91 * sigma;
		    t93 = t80 * t84;
		    t95 = t92 * .1362107888567592 * t93;
		    t96 = t95 + 1.;
/* Computing 2nd power */
		    d__1 = t90;
		    t98 = d__1 * d__1;
		    t99 = 1 / t98;
/* Computing 2nd power */
		    d__1 = sigma;
		    t100 = d__1 * d__1;
		    t101 = t99 * t100;
/* Computing 2nd power */
		    d__1 = t78;
		    t102 = d__1 * d__1;
		    t103 = 1 / t102;
/* Computing 2nd power */
		    d__1 = t82;
		    t104 = d__1 * d__1;
		    t106 = 1 / t104 / t48;
		    t107 = t103 * t106;
		    t110 = t95 + 1. + t101 * .01855337900098064 * t107;
		    t111 = 1 / t110;
		    t115 = t81 * .1362107888567592 * t84 * t96 * t111 + 1.;
		    t116 = log(t115);
		    t118 = t79 * .0310906908696549 * t116;
		    zk[i__] = rho * (-t21 + t55 + t72 + t118);
		    t124 = t38 * 1.333333333333333 * t4 - t42 * 
			    1.333333333333333 * t4;
		    t126 = t33 * t124 * t52;
		    t127 = t126 * .037995525;
		    t128 = t45 * t35;
		    t129 = t44 * t128;
		    t131 = t33 * t129 * t49;
		    t132 = t131 * .1519821;
		    t133 = t69 * t124;
		    t134 = t133 * t50;
		    t135 = t134 * 1.923661050931536;
		    t136 = t128 * t49;
		    t137 = t70 * t136;
		    t138 = t137 * 7.694644203726145;
		    t139 = t78 * t116;
		    t140 = 1 / t38;
		    t143 = 1 / t42;
		    t146 = t140 * .3333333333333333 * t4 - t143 * 
			    .3333333333333333 * t4;
		    t147 = t139 * t146;
		    t148 = t147 * .09327207260896469;
		    t149 = sigma * t86;
		    t150 = t149 * t84;
		    t151 = t96 * t111;
		    t152 = t151 * t146;
		    t155 = t99 * sigma;
		    t156 = t155 * t80;
		    t157 = t127 - t132 + t135 + t138;
		    t160 = t85 * t103;
		    t163 = t157 * -32.16396844291482 * t86 + t160 * 
			    96.49190532874446 * t146;
		    t164 = t84 * t163;
		    t165 = t164 * t89;
		    t167 = t156 * .1362107888567592 * t165;
		    t168 = t86 * t84;
		    t171 = t92 * .2724215777135184 * t168 * t146;
		    t172 = -t167 - t171;
		    t177 = t81 * t84;
/* Computing 2nd power */
		    d__1 = t110;
		    t178 = d__1 * d__1;
		    t179 = 1 / t178;
		    t180 = t96 * t179;
		    t182 = 1 / t98 / t90;
		    t183 = t182 * t100;
		    t184 = t183 * t103;
		    t185 = t106 * t163;
		    t186 = t185 * t89;
		    t190 = 1 / t102 / t77;
		    t191 = t190 * t106;
		    t192 = t191 * t146;
		    t195 = -t167 - t171 - t184 * .03710675800196129 * t186 - 
			    t101 * .07421351600392257 * t192;
		    t196 = t180 * t195;
		    t199 = t150 * -.2724215777135184 * t152 + t81 * 
			    .1362107888567592 * t84 * t172 * t111 - t177 * 
			    .1362107888567592 * t196;
		    t200 = t79 * t199;
		    t201 = 1 / t115;
		    t202 = t200 * t201;
		    t203 = t202 * .0310906908696549;
		    t206 = 1 / t13;
		    t207 = 1 / t47;
		    t208 = t206 * t207;
		    t209 = t208 * t19;
		    t210 = t209 * .002747773264188438;
/* Computing 2nd power */
		    d__1 = t15;
		    t211 = d__1 * d__1;
		    t212 = 1 / t211;
		    t213 = t7 * t212;
/* Computing 2nd power */
		    d__1 = t8;
		    t214 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t214;
		    t215 = d__1 * d__1;
		    t216 = t215 * t8;
		    t217 = 1 / t216;
		    t218 = t217 * t207;
		    t221 = 1 / t11;
		    t222 = t221 * t207;
		    t224 = 1 / t5;
		    t225 = t224 * t207;
		    t227 = t218 * -.99709173929518 - t208 * .7418564737168958 
			    - t222 * .4002143174996817 - t225 * 
			    .1264669550498372;
		    t228 = 1 / t18;
		    t230 = t213 * t227 * t228;
		    t231 = t230 * 1.;
		    t232 = t32 * t44;
		    t233 = t232 * t52;
		    t234 = t208 * t233;
		    t235 = t234 * 8.7407428755417e-4;
/* Computing 2nd power */
		    d__1 = t28;
		    t236 = d__1 * d__1;
		    t237 = 1 / t236;
		    t238 = t23 * t237;
		    t243 = t218 * -1.35956911724794 - t208 * 
			    .7491972878592054 - t222 * .2150486235638321 - 
			    t225 * .1274341730084892;
		    t244 = t238 * t243;
		    t245 = 1 / t31;
		    t246 = t245 * t44;
		    t247 = t246 * t52;
		    t248 = t244 * t247;
		    t249 = t248 * 1.125;
		    t250 = t38 * t35;
		    t253 = t42 * t35;
		    t256 = t250 * -1.333333333333333 * t207 + t253 * 
			    1.333333333333333 * t207;
		    t258 = t33 * t256 * t52;
		    t259 = t258 * .037995525;
		    t260 = t44 * t46;
		    t261 = t48 * rho;
		    t262 = 1 / t261;
		    t264 = t33 * t260 * t262;
		    t265 = t264 * .1519821;
/* Computing 2nd power */
		    d__1 = t62;
		    t268 = d__1 * d__1;
		    t269 = 1 / t268;
		    t270 = t57 * t269;
		    t275 = t218 * -1.853395810515781 - t208 * 
			    1.28158207914907 - t222 * .8223668877838045 - 
			    t225 * .1603914194192128;
		    t276 = 1 / t65;
		    t282 = t208 * .001321039893133927 * t66 + t270 * 1. * 
			    t275 * t276 - t209 * .002747773264188438 - t230 * 
			    1.;
		    t283 = t282 * t44;
		    t284 = t283 * t50;
		    t285 = t284 * 1.923661050931536;
		    t286 = t69 * t256;
		    t287 = t286 * t50;
		    t288 = t287 * 1.923661050931536;
		    t289 = t46 * t262;
		    t290 = t70 * t289;
		    t291 = t290 * 7.694644203726145;
		    t292 = t140 * t35;
		    t295 = t143 * t35;
		    t298 = t292 * -.3333333333333333 * t207 + t295 * 
			    .3333333333333333 * t207;
		    t299 = t139 * t298;
		    t301 = t151 * t298;
		    t304 = t47 * rho;
		    t306 = 1 / t82 / t304;
		    t311 = t210 + t231 - t235 - t249 + t259 + t265 + t285 + 
			    t288 - t291;
		    t316 = t311 * -32.16396844291482 * t86 + t160 * 
			    96.49190532874446 * t298;
		    t317 = t84 * t316;
		    t318 = t317 * t89;
		    t320 = t156 * .1362107888567592 * t318;
		    t323 = t92 * .2724215777135184 * t168 * t298;
		    t324 = t80 * t306;
		    t326 = t92 * .3178251739991049 * t324;
		    t327 = -t320 - t323 - t326;
		    t332 = t106 * t316;
		    t333 = t332 * t89;
		    t336 = t191 * t298;
		    t340 = 1 / t104 / t261;
		    t341 = t103 * t340;
		    t344 = -t320 - t323 - t326 - t184 * .03710675800196129 * 
			    t333 - t101 * .07421351600392257 * t336 - t101 * 
			    .08658243533790967 * t341;
		    t345 = t180 * t344;
		    t348 = t150 * -.2724215777135184 * t301 - t81 * 
			    .3178251739991049 * t306 * t96 * t111 + t81 * 
			    .1362107888567592 * t84 * t327 * t111 - t177 * 
			    .1362107888567592 * t345;
		    t349 = t79 * t348;
		    t350 = t349 * t201;
		    t352 = t210 + t231 - t235 - t249 + t259 + t265 + t285 + 
			    t288 - t291 + t299 * .09327207260896469 + t350 * 
			    .0310906908696549;
		    t353 = rho * t352;
		    vrhoa[i__] = rho * (t127 - t132 + t135 + t138 + t148 + 
			    t203) - t21 + t55 + t72 + t118 + t353;
		    t354 = -t124;
		    t356 = t33 * t354 * t52;
		    t357 = t356 * .037995525;
		    t358 = t131 * .1519821;
		    t359 = t69 * t354;
		    t360 = t359 * t50;
		    t361 = t360 * 1.923661050931536;
		    t362 = t137 * 7.694644203726145;
		    t363 = -t146;
		    t364 = t139 * t363;
		    t365 = t364 * .09327207260896469;
		    t366 = t151 * t363;
		    t369 = t357 + t358 + t361 - t362;
		    t374 = t369 * -32.16396844291482 * t86 + t160 * 
			    96.49190532874446 * t363;
		    t375 = t84 * t374;
		    t376 = t375 * t89;
		    t378 = t156 * .1362107888567592 * t376;
		    t381 = t92 * .2724215777135184 * t168 * t363;
		    t382 = -t378 - t381;
		    t387 = t106 * t374;
		    t388 = t387 * t89;
		    t391 = t191 * t363;
		    t394 = -t378 - t381 - t184 * .03710675800196129 * t388 - 
			    t101 * .07421351600392257 * t391;
		    t395 = t180 * t394;
		    t398 = t150 * -.2724215777135184 * t366 + t81 * 
			    .1362107888567592 * t84 * t382 * t111 - t177 * 
			    .1362107888567592 * t395;
		    t399 = t79 * t398;
		    t400 = t399 * t201;
		    t401 = t400 * .0310906908696549;
		    vrhob[i__] = rho * (t357 + t358 + t361 - t362 + t365 + 
			    t401) - t21 + t55 + t72 + t118 + t353;
		    t404 = rho * t79;
		    t407 = sigma * t103;
		    t412 = t91 * t80;
		    t417 = t412 * .1362107888567592 * t84 + t155 * 
			    .03710675800196129 * t107;
		    t418 = t180 * t417;
		    t421 = t93 * .1362107888567592 * t151 + t407 * 
			    .01855337900098064 * t106 * t91 * t111 - t177 * 
			    .1362107888567592 * t418;
		    t423 = t404 * t421 * t201;
		    vsigmaaa[i__] = t423 * .0310906908696549;
		    vsigmaab[i__] = t423 * .06218138173930979;
		    vsigmabb[i__] = vsigmaaa[i__];
		    t427 = t208 * 8.7407428755417e-4 * t32 * t124 * t52;
		    t428 = t48 * t47;
		    t429 = 1 / t428;
		    t432 = t206 * t429 * t232 * t128;
		    t433 = t432 * .00349629715021668;
		    t437 = t244 * 1.125 * t245 * t124 * t52;
		    t439 = t244 * t246 * t136;
		    t440 = t439 * 4.5;
		    t441 = 1 / t73;
		    t443 = 1 / t304;
		    t448 = 1 / t75;
		    t454 = t441 * -.4444444444444444 * t35 * t443 - t38 * 
			    1.333333333333333 * t207 - t448 * 
			    .4444444444444444 * t35 * t443 + t42 * 
			    1.333333333333333 * t207;
		    t457 = t33 * .037995525 * t454 * t52;
		    t460 = t33 * t256 * t128 * t49;
		    t461 = t460 * .1519821;
		    t465 = t33 * .1519821 * t124 * t46 * t262;
		    t467 = t33 * t129 * t262;
		    t468 = t467 * .6079284;
		    t471 = t282 * 1.923661050931536 * t124 * t50;
		    t472 = t283 * t136;
		    t473 = t472 * 7.694644203726145;
		    t476 = t69 * 1.923661050931536 * t454 * t50;
		    t477 = t286 * t136;
		    t478 = t477 * 7.694644203726145;
		    t480 = t133 * 7.694644203726145 * t289;
		    t482 = t70 * t128 * t262;
		    t483 = t482 * 30.77857681490458;
		    t484 = t77 * t116;
		    t485 = t298 * t146;
		    t487 = t484 * .1865441452179294 * t485;
		    t488 = t78 * t199;
		    t489 = t201 * t298;
		    t491 = t488 * .09327207260896469 * t489;
		    t492 = 1 / t39;
		    t498 = 1 / t43;
		    t504 = t492 * .1111111111111111 * t35 * t443 - t140 * 
			    .3333333333333333 * t207 + t498 * 
			    .1111111111111111 * t35 * t443 + t143 * 
			    .3333333333333333 * t207;
		    t506 = t139 * .09327207260896469 * t504;
		    t507 = t78 * t348;
		    t508 = t201 * t146;
		    t510 = t507 * .09327207260896469 * t508;
		    t511 = t407 * t84;
		    t514 = t511 * .8172647331405553 * t151 * t485;
		    t515 = t172 * t111;
		    t518 = t150 * .2724215777135184 * t515 * t298;
		    t522 = t150 * .2724215777135184 * t180 * t298 * t195;
		    t525 = t150 * .2724215777135184 * t151 * t504;
		    t526 = t149 * t306;
		    t528 = t526 * .6356503479982097 * t152;
		    t532 = t81 * .3178251739991049 * t306 * t172 * t111;
		    t533 = t81 * t306;
		    t535 = t533 * .3178251739991049 * t196;
		    t536 = t327 * t111;
		    t539 = t150 * .2724215777135184 * t536 * t146;
		    t540 = t182 * sigma;
		    t541 = t540 * t80;
/* Computing 2nd power */
		    d__1 = t89;
		    t542 = d__1 * d__1;
		    t543 = t542 * t163;
		    t546 = t541 * .2724215777135184 * t317 * t543;
		    t547 = t155 * t86;
		    t548 = t89 * t146;
		    t551 = t547 * .2724215777135184 * t317 * t548;
		    t552 = -t427 + t433 - t437 + t440 + t457 - t461 + t465 + 
			    t468 + t471 + t473 + t476 + t478 - t480 - t483;
		    t555 = t311 * t103;
		    t557 = t555 * 96.49190532874446 * t146;
		    t558 = t157 * t103;
		    t560 = t558 * 96.49190532874446 * t298;
		    t561 = t85 * t190;
		    t563 = t561 * 385.9676213149779 * t485;
		    t565 = t160 * 96.49190532874446 * t504;
		    t566 = t552 * -32.16396844291482 * t86 + t557 + t560 - 
			    t563 + t565;
		    t570 = t156 * .1362107888567592 * t84 * t566 * t89;
		    t571 = t163 * t89;
		    t574 = t156 * .1362107888567592 * t317 * t571;
		    t575 = t84 * t298;
		    t578 = t547 * .2724215777135184 * t575 * t571;
		    t579 = t92 * t103;
		    t582 = t579 * .8172647331405553 * t575 * t146;
		    t585 = t92 * .2724215777135184 * t168 * t504;
		    t589 = t156 * .3178251739991049 * t306 * t163 * t89;
		    t590 = t86 * t306;
		    t593 = t92 * .6356503479982097 * t590 * t146;
		    t599 = t327 * t179;
		    t602 = t177 * .1362107888567592 * t599 * t195;
		    t606 = t150 * .2724215777135184 * t180 * t344 * t146;
		    t607 = t172 * t179;
		    t610 = t177 * .1362107888567592 * t607 * t344;
		    t613 = t96 / t178 / t110;
		    t617 = t177 * .2724215777135184 * t613 * t344 * t195;
/* Computing 2nd power */
		    d__1 = t98;
		    t618 = d__1 * d__1;
		    t621 = 1 / t618 * t100 * t103;
		    t624 = t621 * .1113202740058839 * t332 * t543;
		    t625 = t183 * t190;
		    t628 = t625 * .1484270320078451 * t332 * t548;
		    t635 = t184 * .03710675800196129 * t332 * t571;
		    t636 = t106 * t298;
		    t639 = t625 * .1484270320078451 * t636 * t571;
		    t641 = 1 / t102 / t78;
		    t642 = t101 * t641;
		    t645 = t642 * .3710675800196129 * t636 * t146;
		    t648 = t101 * .07421351600392257 * t191 * t504;
		    t652 = t184 * .1731648706758193 * t340 * t163 * t89;
		    t653 = t190 * t340;
		    t656 = t101 * .3463297413516387 * t653 * t146;
		    t657 = t546 + t551 - t570 - t574 + t578 + t582 - t585 + 
			    t589 + t593 + t624 + t628 - t184 * 
			    .03710675800196129 * t106 * t566 * t89 - t635 + 
			    t639 + t645 - t648 + t652 + t656;
		    t661 = t514 - t518 + t522 - t525 + t528 - t532 + t535 - 
			    t539 + t81 * .1362107888567592 * t84 * (t546 + 
			    t551 - t570 - t574 + t578 + t582 - t585 + t589 + 
			    t593) * t111 - t602 + t606 - t610 + t617 - t177 * 
			    .1362107888567592 * t180 * t657;
/* Computing 2nd power */
		    d__1 = t115;
		    t665 = d__1 * d__1;
		    t666 = 1 / t665;
		    t669 = t349 * .0310906908696549 * t666 * t199;
		    t670 = -t427 + t433 - t437 + t440 + t457 - t461 + t465 + 
			    t468 + t471 + t473 + t476 + t478 - t480 - t483 + 
			    t487 + t491 + t506 + t510 + t79 * 
			    .0310906908696549 * t661 * t201 - t669;
		    t671 = rho * t670;
		    t677 = t441 * .4444444444444444 * t207 + t448 * 
			    .4444444444444444 * t207;
		    t680 = t33 * .037995525 * t677 * t52;
		    t683 = t33 * t124 * t128 * t49;
		    t684 = t683 * .3039642;
		    t687 = t33 * t44 * t45 * t49;
		    t688 = t687 * .4559463;
		    t691 = t69 * 1.923661050931536 * t677 * t50;
		    t692 = t133 * t136;
		    t693 = t692 * 15.38928840745229;
		    t695 = t70 * t45 * t49;
		    t696 = t695 * 23.08393261117844;
/* Computing 2nd power */
		    d__1 = t146;
		    t697 = d__1 * d__1;
		    t706 = t492 * -.1111111111111111 * t207 - t498 * 
			    .1111111111111111 * t207;
		    t708 = t139 * .09327207260896469 * t706;
		    t721 = t150 * .2724215777135184 * t151 * t706;
/* Computing 2nd power */
		    d__1 = t163;
		    t722 = d__1 * d__1;
		    t723 = t84 * t722;
		    t726 = t541 * .2724215777135184 * t723 * t542;
		    t729 = t547 * .5448431554270369 * t164 * t548;
		    t738 = t160 * 96.49190532874446 * t706;
		    t739 = (t680 - t684 - t688 + t691 + t693 + t696) * 
			    -32.16396844291482 * t86 + t558 * 
			    192.9838106574889 * t146 - t561 * 
			    385.9676213149779 * t697 + t738;
		    t743 = t156 * .1362107888567592 * t84 * t739 * t89;
		    t746 = t156 * .1362107888567592 * t723 * t89;
		    t747 = t103 * t84;
		    t750 = t92 * .8172647331405553 * t747 * t697;
		    t753 = t92 * .2724215777135184 * t168 * t706;
/* Computing 2nd power */
		    d__1 = t195;
		    t762 = d__1 * d__1;
		    t766 = t106 * t722;
		    t780 = t641 * t106;
		    t786 = t101 * .07421351600392257 * t191 * t706;
		    t787 = t726 + t729 - t743 - t746 + t750 - t753 + t621 * 
			    .1113202740058839 * t766 * t542 + t625 * 
			    .2968540640156903 * t185 * t548 - t184 * 
			    .03710675800196129 * t106 * t739 * t89 - t184 * 
			    .03710675800196129 * t766 * t89 + t101 * 
			    .3710675800196129 * t780 * t697 - t786;
/* Computing 2nd power */
		    d__1 = t199;
		    t795 = d__1 * d__1;
		    t799 = t680 - t684 - t688 + t691 + t693 + t696 + t484 * 
			    .1865441452179294 * t697 + t488 * 
			    .1865441452179294 * t508 + t708 + t79 * 
			    .0310906908696549 * (t511 * .8172647331405553 * 
			    t151 * t697 - t150 * .5448431554270369 * t515 * 
			    t146 + t150 * .5448431554270369 * t180 * t146 * 
			    t195 - t721 + t81 * .1362107888567592 * t84 * (
			    t726 + t729 - t743 - t746 + t750 - t753) * t111 - 
			    t177 * .2724215777135184 * t607 * t195 + t177 * 
			    .2724215777135184 * t613 * t762 - t177 * 
			    .1362107888567592 * t180 * t787) * t201 - t79 * 
			    .0310906908696549 * t795 * t666;
		    t801 = t439 * 4.5;
		    t802 = -t427 - t437 + t457 + t465 + t433 + t801 - t461 + 
			    t468 + t471 + t476 - t480 + t473 + t478 - t483;
		    t805 = t802 * -32.16396844291482 * t86 + t560 + t557 - 
			    t563 + t565;
		    t809 = t156 * .1362107888567592 * t84 * t805 * t89;
		    t819 = t546 + t578 + t589 - t809 - t574 + t551 + t582 + 
			    t593 - t585 + t624 + t639 + t652 - t184 * 
			    .03710675800196129 * t106 * t805 * t89 - t635 + 
			    t628 + t645 + t656 - t648;
		    t823 = t514 + t528 - t539 + t606 - t525 - t518 - t532 + 
			    t81 * .1362107888567592 * t84 * (t546 + t578 + 
			    t589 - t809 - t574 + t551 + t582 + t593 - t585) * 
			    t111 - t610 + t522 + t535 - t602 + t617 - t177 * 
			    .1362107888567592 * t180 * t819;
		    t827 = -t427 - t437 + t457 + t465 + t433 + t801 - t461 + 
			    t468 + t471 + t476 - t480 + t473 + t478 - t483 + 
			    t487 + t510 + t506 + t491 + t79 * 
			    .0310906908696549 * t823 * t201 - t669;
		    t828 = rho * t827;
		    t829 = t209 * .005495546528376876;
		    t831 = t230 * 2.;
		    t836 = t290 * 15.38928840745229;
		    t837 = t299 * .1865441452179294;
		    t838 = t287 * 3.847322101863073;
		    t839 = t258 * .07599105;
		    t840 = t234 * .00174814857510834;
		    t841 = t248 * 2.25;
		    t842 = t264 * .3039642;
		    t843 = t284 * 3.847322101863073;
		    t844 = t350 * .06218138173930979;
		    t848 = 1 / t13 / t4 * t49;
		    t849 = t848 * t19;
		    t850 = t849 * .001831848842792292;
/* Computing 2nd power */
		    d__1 = t348;
		    t851 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t227;
		    t858 = d__1 * d__1;
		    t860 = t7 / t211 / t15 * t858 * t228;
		    t861 = t860 * 2.;
		    t864 = 1 / t216 / t4 * t49;
		    t866 = t217 * t443;
		    t869 = t206 * t443;
		    t873 = 1 / t11 / t4 * t49;
		    t875 = t221 * t443;
		    t879 = 1 / t5 / t4 * t49;
		    t881 = t224 * t443;
		    t885 = t213 * (t864 * -.8309097827459833 + t866 * 
			    1.99418347859036 - t848 * .4945709824779306 + 
			    t869 * 1.483712947433792 - t873 * 
			    .2001071587498409 + t875 * .8004286349993634 - 
			    t879 * .04215565168327908 + t881 * 
			    .2529339100996745) * t228;
		    t886 = t885 * 1.;
/* Computing 2nd power */
		    d__1 = t211;
		    t887 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t18;
		    t890 = d__1 * d__1;
		    t893 = t7 / t887 * t858 / t890;
		    t894 = t893 * 16.08197949869254;
		    t898 = t244 * 2.25 * t245 * t256 * t52;
		    t900 = t869 * .00174814857510834 * t233;
/* Computing 2nd power */
		    d__1 = t236;
		    t901 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t243;
		    t904 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t31;
		    t906 = d__1 * d__1;
		    t911 = t23 * 33.30984372501762 / t901 * t904 / t906 * t44 
			    * t52;
		    t914 = t244 * 9. * t246 * t289;
		    t917 = t208 * t212 * t227 * t228;
		    t918 = t917 * .08837926660346786;
		    t919 = t869 * t19;
		    t920 = t919 * .005495546528376876;
		    t923 = t33 * .7599105 * t260 * t429;
		    t927 = t33 * .3039642 * t256 * t46 * t262;
/* Computing 2nd power */
		    d__1 = t275;
		    t939 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t268;
		    t955 = d__1 * d__1;
/* Computing 2nd power */
		    d__1 = t65;
		    t958 = d__1 * d__1;
		    t969 = t848 * 8.806932620892844e-4 * t66 - t869 * 
			    .002642079786267853 * t66 - t208 * 
			    .08497974591333914 * t269 * t275 * t276 - t57 * 
			    2. / t268 / t62 * t939 * t276 + t270 * 1. * (t864 
			    * -1.544496508763151 + t866 * 3.706791621031562 - 
			    t848 * .854388052766047 + t869 * 
			    2.563164158298141 - t873 * .4111834438919023 + 
			    t875 * 1.644733775567609 - t879 * 
			    .05346380647307093 + t881 * .3207828388384256) * 
			    t276 + t57 * 32.16395899738507 / t955 * t939 / 
			    t958 - t849 * .001831848842792292 + t919 * 
			    .005495546528376876 + t917 * .08837926660346786 + 
			    t860 * 2. - t885 * 1. - t893 * 16.08197949869254;
		    t972 = t969 * 1.923661050931536 * t44 * t50;
		    t983 = t441 * .4444444444444444 * t45 * t49 + t250 * 
			    2.666666666666667 * t443 + t448 * 
			    .4444444444444444 * t45 * t49 - t253 * 
			    2.666666666666667 * t443;
		    t986 = t33 * .037995525 * t983 * t52;
		    t987 = t850 - t79 * .0310906908696549 * t851 * t666 - 
			    t861 + t886 + t894 - t898 + t900 - t911 - t914 - 
			    t918 - t920 - t923 + t927 + t972 + t986;
		    t989 = t286 * 15.38928840745229 * t289;
		    t992 = t69 * 1.923661050931536 * t983 * t50;
		    t994 = t848 * 5.8271619170278e-4 * t233;
		    t1006 = t238 * 1.125 * (t864 * -1.132974264373283 + t866 *
			     2.71913823449588 - t848 * .4994648585728036 + 
			    t869 * 1.498394575718411 - t873 * 
			    .1075243117819161 + t875 * .4300972471276643 - 
			    t879 * .04247805766949639 + t881 * 
			    .2548683460169784) * t247;
		    t1008 = 1 / t48 / t304;
		    t1012 = t206 * .00699259430043336 * t1008 * t232 * t46;
		    t1017 = t208 * .05176049408441869 * t237 * t243 * t245 * 
			    t53;
		    t1023 = t23 * 2.25 / t236 / t28 * t904 * t247;
/* Computing 2nd power */
		    d__1 = t298;
		    t1024 = d__1 * d__1;
		    t1029 = t282 * 3.847322101863073 * t256 * t50;
		    t1031 = t283 * 15.38928840745229 * t289;
		    t1034 = t70 * 38.47322101863073 * t46 * t429;
		    t1038 = t208 * .00174814857510834 * t32 * t256 * t52;
		    t1054 = 1 / t82 / t48;
		    t1073 = t492 * -.1111111111111111 * t45 * t49 + t292 * 
			    .6666666666666667 * t443 - t498 * 
			    .1111111111111111 * t45 * t49 - t295 * 
			    .6666666666666667 * t443;
/* Computing 2nd power */
		    d__1 = t316;
		    t1077 = d__1 * d__1;
		    t1078 = t84 * t1077;
		    t1081 = t541 * .2724215777135184 * t1078 * t542;
		    t1082 = t89 * t298;
		    t1085 = t547 * .5448431554270369 * t317 * t1082;
		    t1089 = t156 * .6356503479982097 * t306 * t316 * t89;
		    t1090 = t850 - t861 + t886 + t894 - t898 + t900 - t911 - 
			    t914 - t918 - t920 - t923 + t927;
		    t1091 = t972 + t986 - t989 + t992 - t994 - t1006 - t1012 
			    + t1017 + t1023 + t1029 - t1031 + t1034 - t1038;
		    t1101 = (t1090 + t1091) * -32.16396844291482 * t86 + t555 
			    * 192.9838106574889 * t298 - t561 * 
			    385.9676213149779 * t1024 + t160 * 
			    96.49190532874446 * t1073;
		    t1105 = t156 * .1362107888567592 * t84 * t1101 * t89;
		    t1108 = t156 * .1362107888567592 * t1078 * t89;
		    t1111 = t92 * .8172647331405553 * t747 * t1024;
		    t1114 = t92 * 1.271300695996419 * t590 * t298;
		    t1117 = t92 * .2724215777135184 * t168 * t1073;
		    t1120 = t92 * 1.059417246663683 * t80 * t1054;
/* Computing 2nd power */
		    d__1 = t344;
		    t1129 = d__1 * d__1;
		    t1150 = t106 * t1077;
		    t1165 = t1111 + t1114 - t1117 + t1120 + t1089 + t1081 + 
			    t1085 - t1108 - t1105 + t101 * .3710675800196129 *
			     t780 * t1024 + t101 * .6926594827032773 * t653 * 
			    t298 + t184 * .3463297413516387 * t340 * t316 * 
			    t89 - t101 * .07421351600392257 * t191 * t1073 - 
			    t184 * .03710675800196129 * t106 * t1101 * t89 + 
			    t621 * .1113202740058839 * t1150 * t542 + t101 * 
			    .4906338002481548 * t103 / t104 / t428 + t625 * 
			    .2968540640156903 * t332 * t1082 - t184 * 
			    .03710675800196129 * t1150 * t89;
		    t1169 = t511 * .8172647331405553 * t151 * t1024 + t526 * 
			    1.271300695996419 * t301 - t150 * 
			    .5448431554270369 * t536 * t298 + t150 * 
			    .5448431554270369 * t180 * t298 * t344 + t533 * 
			    .6356503479982097 * t345 + t81 * 
			    1.059417246663683 * t1054 * t96 * t111 - t81 * 
			    .6356503479982097 * t306 * t327 * t111 - t150 * 
			    .2724215777135184 * t151 * t1073 + t81 * 
			    .1362107888567592 * t84 * (t1081 + t1085 + t1089 
			    - t1105 - t1108 + t1111 + t1114 - t1117 + t1120) *
			     t111 - t177 * .2724215777135184 * t599 * t344 + 
			    t177 * .2724215777135184 * t613 * t1129 - t177 * 
			    .1362107888567592 * t180 * t1165;
		    t1177 = -t989 + t992 - t994 - t1006 - t1012 + t1017 + 
			    t1023 + t484 * .1865441452179294 * t1024 + t1029 
			    - t1031 + t1034 - t1038 + t79 * .0310906908696549 
			    * t1169 * t201 + t507 * .1865441452179294 * t489 
			    + t139 * .09327207260896469 * t1073;
		    t1179 = rho * (t987 + t1177);
		    t1180 = -t836 + t837 + t838 + t839 - t840 - t841 + t842 + 
			    t843 + t844 + t134 * 3.847322101863073 + t1179;
		    v2rhoa2[i__] = t671 - t131 * .3039642 + rho * t799 + t828 
			    + t829 + t202 * .06218138173930979 + t831 + t137 *
			     15.38928840745229 + t147 * .1865441452179294 + 
			    t126 * .07599105 + t1180;
		    t1184 = t33 * t354 * t128 * t49;
		    t1185 = t1184 * .3039642;
		    t1186 = t359 * t136;
		    t1187 = t1186 * 15.38928840745229;
/* Computing 2nd power */
		    d__1 = t363;
		    t1188 = d__1 * d__1;
		    t1191 = t78 * t398;
		    t1192 = t201 * t363;
		    t1198 = t382 * t111;
/* Computing 2nd power */
		    d__1 = t374;
		    t1206 = d__1 * d__1;
		    t1207 = t84 * t1206;
		    t1210 = t541 * .2724215777135184 * t1207 * t542;
		    t1211 = t89 * t363;
		    t1214 = t547 * .5448431554270369 * t375 * t1211;
		    t1218 = t369 * t103;
		    t1223 = (t680 + t1185 - t688 + t691 - t1187 + t696) * 
			    -32.16396844291482 * t86 + t1218 * 
			    192.9838106574889 * t363 - t561 * 
			    385.9676213149779 * t1188 + t738;
		    t1227 = t156 * .1362107888567592 * t84 * t1223 * t89;
		    t1230 = t156 * .1362107888567592 * t1207 * t89;
		    t1233 = t92 * .8172647331405553 * t747 * t1188;
		    t1239 = t382 * t179;
/* Computing 2nd power */
		    d__1 = t394;
		    t1243 = d__1 * d__1;
		    t1247 = t106 * t1206;
		    t1264 = t1210 + t1214 - t1227 - t1230 + t1233 - t753 + 
			    t621 * .1113202740058839 * t1247 * t542 + t625 * 
			    .2968540640156903 * t387 * t1211 - t184 * 
			    .03710675800196129 * t106 * t1223 * t89 - t184 * 
			    .03710675800196129 * t1247 * t89 + t101 * 
			    .3710675800196129 * t780 * t1188 - t786;
/* Computing 2nd power */
		    d__1 = t398;
		    t1272 = d__1 * d__1;
		    t1276 = t680 + t1185 - t688 + t691 - t1187 + t696 + t484 *
			     .1865441452179294 * t1188 + t1191 * 
			    .1865441452179294 * t1192 + t708 + t79 * 
			    .0310906908696549 * (t511 * .8172647331405553 * 
			    t151 * t1188 - t150 * .5448431554270369 * t1198 * 
			    t363 + t150 * .5448431554270369 * t180 * t363 * 
			    t394 - t721 + t81 * .1362107888567592 * t84 * (
			    t1210 + t1214 - t1227 - t1230 + t1233 - t753) * 
			    t111 - t177 * .2724215777135184 * t1239 * t394 + 
			    t177 * .2724215777135184 * t613 * t1243 - t177 * 
			    .1362107888567592 * t180 * t1264) * t201 - t79 * 
			    .0310906908696549 * t1272 * t666;
		    t1285 = t208 * 8.7407428755417e-4 * t32 * t354 * t52;
		    t1286 = t432 * .00349629715021668;
		    t1290 = t244 * 1.125 * t245 * t354 * t52;
		    t1291 = t439 * 4.5;
		    t1292 = -t454;
		    t1295 = t33 * .037995525 * t1292 * t52;
		    t1296 = t460 * .1519821;
		    t1300 = t33 * .1519821 * t354 * t46 * t262;
		    t1301 = t467 * .6079284;
		    t1304 = t282 * 1.923661050931536 * t354 * t50;
		    t1305 = t472 * 7.694644203726145;
		    t1308 = t69 * 1.923661050931536 * t1292 * t50;
		    t1309 = t477 * 7.694644203726145;
		    t1311 = t359 * 7.694644203726145 * t289;
		    t1312 = t482 * 30.77857681490458;
		    t1313 = t298 * t363;
		    t1318 = -t504;
		    t1347 = t542 * t374;
		    t1350 = t541 * .2724215777135184 * t317 * t1347;
		    t1353 = t547 * .2724215777135184 * t317 * t1211;
		    t1354 = -t1285 - t1286 - t1290 - t1291 + t1295 + t1296 + 
			    t1300 - t1301 + t1304 - t1305 + t1308 - t1309 - 
			    t1311 + t1312;
		    t1365 = t1354 * -32.16396844291482 * t86 + t555 * 
			    96.49190532874446 * t363 + t1218 * 
			    96.49190532874446 * t298 - t561 * 
			    385.9676213149779 * t1313 + t160 * 
			    96.49190532874446 * t1318;
		    t1369 = t156 * .1362107888567592 * t84 * t1365 * t89;
		    t1370 = t374 * t89;
		    t1373 = t156 * .1362107888567592 * t317 * t1370;
		    t1376 = t547 * .2724215777135184 * t575 * t1370;
		    t1379 = t579 * .8172647331405553 * t575 * t363;
		    t1382 = t92 * .2724215777135184 * t168 * t1318;
		    t1386 = t156 * .3178251739991049 * t306 * t374 * t89;
		    t1389 = t92 * .6356503479982097 * t590 * t363;
		    t1438 = t1350 + t1353 - t1369 - t1373 + t1376 + t1379 - 
			    t1382 + t1386 + t1389 + t621 * .1113202740058839 *
			     t332 * t1347 + t625 * .1484270320078451 * t332 * 
			    t1211 - t184 * .03710675800196129 * t106 * t1365 *
			     t89 - t184 * .03710675800196129 * t332 * t1370 + 
			    t625 * .1484270320078451 * t636 * t1370 + t642 * 
			    .3710675800196129 * t636 * t363 - t101 * 
			    .07421351600392257 * t191 * t1318 + t184 * 
			    .1731648706758193 * t340 * t374 * t89 + t101 * 
			    .3463297413516387 * t653 * t363;
		    t1442 = t511 * .8172647331405553 * t151 * t1313 - t150 * 
			    .2724215777135184 * t1198 * t298 + t150 * 
			    .2724215777135184 * t180 * t298 * t394 - t150 * 
			    .2724215777135184 * t151 * t1318 + t526 * 
			    .6356503479982097 * t366 - t81 * 
			    .3178251739991049 * t306 * t382 * t111 + t533 * 
			    .3178251739991049 * t395 - t150 * 
			    .2724215777135184 * t536 * t363 + t81 * 
			    .1362107888567592 * t84 * (t1350 + t1353 - t1369 
			    - t1373 + t1376 + t1379 - t1382 + t1386 + t1389) *
			     t111 - t177 * .1362107888567592 * t599 * t394 + 
			    t150 * .2724215777135184 * t180 * t344 * t363 - 
			    t177 * .1362107888567592 * t1239 * t344 + t177 * 
			    .2724215777135184 * t613 * t344 * t394 - t177 * 
			    .1362107888567592 * t180 * t1438;
		    t1446 = t666 * t398;
		    t1449 = -t1285 - t1286 - t1290 - t1291 + t1295 + t1296 + 
			    t1300 - t1301 + t1304 - t1305 + t1308 - t1309 - 
			    t1311 + t1312 + t484 * .1865441452179294 * t1313 
			    + t1191 * .09327207260896469 * t489 + t139 * 
			    .09327207260896469 * t1318 + t507 * 
			    .09327207260896469 * t1192 + t79 * 
			    .0310906908696549 * t1442 * t201 - t349 * 
			    .0310906908696549 * t1446;
		    t1450 = rho * t1449;
		    v2rhob2[i__] = t131 * .3039642 + t829 + rho * t1276 + 
			    t831 - t137 * 15.38928840745229 - t836 + t837 + 
			    t838 + t839 - t840 - t841 + t842 + t843 + t844 + 
			    t360 * 3.847322101863073 + t356 * .07599105 + 
			    t1450 * 2 + t1179 + t400 * .06218138173930979 + 
			    t364 * .1865441452179294;
		    t1455 = -t677;
		    t1458 = t33 * .037995525 * t1455 * t52;
		    t1459 = t683 * .1519821;
		    t1460 = t1184 * .1519821;
		    t1461 = t687 * .4559463;
		    t1464 = t69 * 1.923661050931536 * t1455 * t50;
		    t1465 = t692 * 7.694644203726145;
		    t1466 = t1186 * 7.694644203726145;
		    t1467 = t695 * 23.08393261117844;
		    t1468 = t146 * t363;
		    t1473 = -t706;
		    t1496 = t541 * .2724215777135184 * t164 * t1347;
		    t1499 = t547 * .2724215777135184 * t164 * t1211;
		    t1511 = (t1458 + t1459 - t1460 + t1461 + t1464 - t1465 + 
			    t1466 - t1467) * -32.16396844291482 * t86 + t558 *
			     96.49190532874446 * t363 + t1218 * 
			    96.49190532874446 * t146 - t561 * 
			    385.9676213149779 * t1468 + t160 * 
			    96.49190532874446 * t1473;
		    t1515 = t156 * .1362107888567592 * t84 * t1511 * t89;
		    t1518 = t156 * .1362107888567592 * t164 * t1370;
		    t1519 = t84 * t146;
		    t1522 = t547 * .2724215777135184 * t1519 * t1370;
		    t1525 = t579 * .8172647331405553 * t1519 * t363;
		    t1528 = t92 * .2724215777135184 * t168 * t1473;
		    t1561 = t106 * t146;
		    t1571 = t1496 + t1499 - t1515 - t1518 + t1522 + t1525 - 
			    t1528 + t621 * .1113202740058839 * t185 * t1347 + 
			    t625 * .1484270320078451 * t185 * t1211 - t184 * 
			    .03710675800196129 * t106 * t1511 * t89 - t184 * 
			    .03710675800196129 * t185 * t1370 + t625 * 
			    .1484270320078451 * t1561 * t1370 + t642 * 
			    .3710675800196129 * t1561 * t363 - t101 * 
			    .07421351600392257 * t191 * t1473;
		    t1575 = t511 * .8172647331405553 * t151 * t1468 - t150 * 
			    .2724215777135184 * t1198 * t146 + t150 * 
			    .2724215777135184 * t180 * t146 * t394 - t150 * 
			    .2724215777135184 * t151 * t1473 - t150 * 
			    .2724215777135184 * t515 * t363 + t81 * 
			    .1362107888567592 * t84 * (t1496 + t1499 - t1515 
			    - t1518 + t1522 + t1525 - t1528) * t111 - t177 * 
			    .1362107888567592 * t607 * t394 + t150 * 
			    .2724215777135184 * t180 * t195 * t363 - t177 * 
			    .1362107888567592 * t1239 * t195 + t177 * 
			    .2724215777135184 * t613 * t195 * t394 - t177 * 
			    .1362107888567592 * t180 * t1571;
		    t1581 = t1458 + t1459 - t1460 + t1461 + t1464 - t1465 + 
			    t1466 - t1467 + t484 * .1865441452179294 * t1468 
			    + t1191 * .09327207260896469 * t508 + t139 * 
			    .09327207260896469 * t1473 + t488 * 
			    .09327207260896469 * t1192 + t79 * 
			    .0310906908696549 * t1575 * t201 - t200 * 
			    .0310906908696549 * t1446;
		    t1584 = t829 + t671 * .5 - t840 + t148 - t841 + rho * 
			    t1581 + t842 + t1179 + t831 + t839 + t843 + t828 *
			     .5;
		    t1586 = -t836 + t838 + t203 + t844 + t357 + t361 + t365 + 
			    t401 + t1450 * 1. + t127 + t135 + t837;
		    v2rhoab[i__] = t1584 + t1586;
		    t1587 = t78 * t421;
		    t1593 = sigma * t190 * t106;
		    t1594 = t91 * t111;
		    t1604 = t99 * t80;
		    t1606 = t1604 * .1362107888567592 * t165;
		    t1607 = t91 * t86;
		    t1609 = t1607 * .2724215777135184 * t1519;
		    t1620 = t407 * t106;
		    t1621 = t91 * t179;
		    t1629 = t540 * t103;
		    t1642 = t666 * t421;
		    t1646 = rho * (t1587 * .09327207260896469 * t508 + t79 * 
			    .0310906908696549 * (t168 * -.2724215777135184 * 
			    t152 - t1593 * .03710675800196129 * t1594 * t146 
			    + t150 * .2724215777135184 * t180 * t146 * t417 + 
			    t93 * .1362107888567592 * t515 + t81 * 
			    .1362107888567592 * t84 * (-t1606 - t1609) * t111 
			    - t177 * .1362107888567592 * t607 * t417 - t93 * 
			    .1362107888567592 * t196 - t1620 * 
			    .01855337900098064 * t1621 * t195 + t177 * 
			    .2724215777135184 * t613 * t195 * t417 - t177 * 
			    .1362107888567592 * t180 * (-t1606 - t1609 - 
			    t1629 * .07421351600392257 * t186 - t155 * 
			    .1484270320078451 * t192)) * t201 - t200 * 
			    .0310906908696549 * t1642);
		    t1648 = t79 * t421 * t201;
		    t1649 = t1648 * .0310906908696549;
		    t1672 = t1604 * .1362107888567592 * t318;
		    t1674 = t1607 * .2724215777135184 * t575;
		    t1676 = t412 * .3178251739991049 * t306;
		    t1704 = t168 * -.2724215777135184 * t301 - t1593 * 
			    .03710675800196129 * t1594 * t298 + t150 * 
			    .2724215777135184 * t180 * t298 * t417 - t324 * 
			    .3178251739991049 * t151 - t407 * 
			    .04329121766895483 * t340 * t91 * t111 + t533 * 
			    .3178251739991049 * t418 + t93 * 
			    .1362107888567592 * t536 + t81 * 
			    .1362107888567592 * t84 * (-t1672 - t1674 - t1676)
			     * t111 - t177 * .1362107888567592 * t599 * t417 
			    - t93 * .1362107888567592 * t345 - t1620 * 
			    .01855337900098064 * t1621 * t344 + t177 * 
			    .2724215777135184 * t613 * t344 * t417 - t177 * 
			    .1362107888567592 * t180 * (-t1672 - t1674 - 
			    t1676 - t1629 * .07421351600392257 * t333 - t155 *
			     .1484270320078451 * t336 - t155 * 
			    .1731648706758193 * t341);
		    t1711 = rho * (t1587 * .09327207260896469 * t489 + t79 * 
			    .0310906908696549 * t1704 * t201 - t349 * 
			    .0310906908696549 * t1642);
		    v2rhoasigmaaa[i__] = t1646 + t1649 + t1711;
		    t1713 = t1648 * .06218138173930979;
		    t1714 = t1711 * 2.;
		    v2rhoasigmaab[i__] = t1646 * 2. + t1713 + t1714;
		    v2rhoasigmabb[i__] = v2rhoasigmaaa[i__];
		    t1729 = t1604 * .1362107888567592 * t376;
		    t1732 = t1607 * .2724215777135184 * t84 * t363;
		    t1765 = rho * (t1587 * .09327207260896469 * t1192 + t79 * 
			    .0310906908696549 * (t168 * -.2724215777135184 * 
			    t366 - t1593 * .03710675800196129 * t1594 * t363 
			    + t150 * .2724215777135184 * t180 * t363 * t417 + 
			    t93 * .1362107888567592 * t1198 + t81 * 
			    .1362107888567592 * t84 * (-t1729 - t1732) * t111 
			    - t177 * .1362107888567592 * t1239 * t417 - t93 * 
			    .1362107888567592 * t395 - t1620 * 
			    .01855337900098064 * t1621 * t394 + t177 * 
			    .2724215777135184 * t613 * t394 * t417 - t177 * 
			    .1362107888567592 * t180 * (-t1729 - t1732 - 
			    t1629 * .07421351600392257 * t388 - t155 * 
			    .1484270320078451 * t391)) * t201 - t399 * 
			    .0310906908696549 * t1642);
		    v2rhobsigmaaa[i__] = t1765 + t1649 + t1711;
		    v2rhobsigmaab[i__] = t1765 * 2. + t1713 + t1714;
		    v2rhobsigmabb[i__] = v2rhobsigmaaa[i__];
/* Computing 2nd power */
		    d__1 = t417;
		    t1774 = d__1 * d__1;
		    t1785 = t404 * (t107 * .03710675800196129 * t1594 - t93 * 
			    .2724215777135184 * t418 - t1620 * 
			    .03710675800196129 * t1621 * t417 + t177 * 
			    .2724215777135184 * t613 * t1774 - sigma * 
			    .005054340779364009 * t641 * t1008 * t180 * t99) *
			     t201;
/* Computing 2nd power */
		    d__1 = t421;
		    t1787 = d__1 * d__1;
		    t1789 = t404 * t1787 * t666;
		    v2sigmaaa2[i__] = t1785 * .0310906908696549 - t1789 * 
			    .0310906908696549;
		    v2sigmaaaab[i__] = t1785 * .06218138173930979 - t1789 * 
			    .06218138173930979;
		    v2sigmaaabb[i__] = v2sigmaaa2[i__];
		    v2sigmaab2[i__] = t1785 * .1243627634786196 - t1789 * 
			    .1243627634786196;
		    v2sigmaabbb[i__] = v2sigmaaaab[i__];
		    v2sigmabb2[i__] = v2sigmaaabb[i__];
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
} /* uks_c_pbe__ */

/* Subroutine */ int rks_c_pbe__(integer *ideriv, integer *npt, doublereal *
	rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa, 
	doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa, 
	doublereal *v2sigmaaa2)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, s2, t2, t3, t5, t6, t9, t11, t20, t21, t23, t31, 
	    t33, t17, t18, t26, t27, t35, t37, t38, t49, t13, t16, t19, t24, 
	    t28, t29, t32, t34, t36, t40, t43, t44, t45, t48, t50, t53, t54, 
	    t56, t57, t60, t61, t77, t80, t84, t85, t89, t91, t96, t98, t52, 
	    t55, t58, t59, t62, t63, t67, t70, t73, t74, t116, t76, t78, t81, 
	    t92, t93, t97, t100, t101, t103, t108, t111, t112, t115, t117, 
	    t121, t124, t125, t132, t133, t136, t151, t153, t155, t158, t160, 
	    t165, t169, t170, t171, t183, t197, t201, t202, t223, t226, t231, 
	    t233, t250, t251, t254, t255, t259, t261, t265, t271, t273, t274, 
	    t275, t276, t282, t287, t288, t289, t293, t297, t298, t299, t302, 
	    t311, t314, rho, t315, t325, t326, t329, t330, t365, t367, t378, 
	    t413, t427, sigma;


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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b2);
		t6 = pow_dd(&t2, &c_b3);
		t9 = sqrt(t2);
/* Computing 2nd power */
		d__1 = t3;
		t11 = d__1 * d__1;
		t17 = log(16.08197949869254 / (t6 * 5.98255043577108 + t3 * 
			2.225569421150687 + t9 * .8004286349993634 + t11 * 
			.1897004325747559) + 1.);
		t18 = (t3 * .1325688999052018 + 1.) * t17;
/* Computing 2nd power */
		d__1 = rho;
		t20 = d__1 * d__1;
		t21 = pow_dd(&rho, &c_b2);
		t23 = 1 / t21 / t20;
		t26 = exp(t18 * 2.000000587336264);
		t27 = t26 - 1.;
		t31 = .1362107888567592 / t27 * sigma * t23;
/* Computing 2nd power */
		d__1 = t27;
		t33 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = sigma;
		t35 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t21;
		t38 = d__1 * d__1;
		t49 = log(sigma * .1362107888567592 * t23 * (t31 + 1.) / (t31 
			+ 1. + .01855337900098064 / t33 * t35 / t38 / t37) + 
			1.);
		zk[i__] = rho * (t18 * -.0621814 + t49 * .0310906908696549);
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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b2);
		t5 = t3 * .1325688999052018 + 1.;
		t6 = pow_dd(&t2, &c_b3);
		t9 = sqrt(t2);
/* Computing 2nd power */
		d__1 = t3;
		t11 = d__1 * d__1;
		t13 = t6 * 5.98255043577108 + t3 * 2.225569421150687 + t9 * 
			.8004286349993634 + t11 * .1897004325747559;
		t16 = 16.08197949869254 / t13 + 1.;
		t17 = log(t16);
		t18 = t5 * t17;
		t19 = t18 * .0621814;
/* Computing 2nd power */
		d__1 = rho;
		t20 = d__1 * d__1;
		t21 = pow_dd(&rho, &c_b2);
		t23 = 1 / t21 / t20;
		t24 = sigma * t23;
		t26 = exp(t18 * 2.000000587336264);
		t27 = t26 - 1.;
		t28 = 1 / t27;
		t29 = t28 * sigma;
		t31 = t29 * .1362107888567592 * t23;
		t32 = t31 + 1.;
/* Computing 2nd power */
		d__1 = t27;
		t33 = d__1 * d__1;
		t34 = 1 / t33;
/* Computing 2nd power */
		d__1 = sigma;
		t35 = d__1 * d__1;
		t36 = t34 * t35;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t21;
		t38 = d__1 * d__1;
		t40 = 1 / t38 / t37;
		t43 = t31 + 1. + t36 * .01855337900098064 * t40;
		t44 = 1 / t43;
		t45 = t32 * t44;
		t48 = t24 * .1362107888567592 * t45 + 1.;
		t49 = log(t48);
		t50 = t49 * .0310906908696549;
		zk[i__] = rho * (-t19 + t50);
		t53 = 1 / t20;
		t54 = 1 / t11 * t53;
		t56 = t54 * .002747773264188438 * t17;
/* Computing 2nd power */
		d__1 = t13;
		t57 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t6;
		t60 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t60;
		t61 = d__1 * d__1;
		t77 = t5 * 1. / t57 * (-.99709173929518 / t61 / t6 * t53 - 
			t54 * .7418564737168958 - .4002143174996817 / t9 * 
			t53 - .1264669550498372 / t3 * t53) / t16;
		t80 = 1 / t21 / t20 / rho;
		t84 = t34 * sigma;
		t85 = t56 + t77;
		t89 = t84 * 4.381079514373337 * t23 * t85 * t26;
		t91 = t29 * .3178251739991049 * t80;
/* Computing 2nd power */
		d__1 = t43;
		t96 = d__1 * d__1;
		t98 = t32 / t96;
		t116 = 1 / t48;
		vrhoa[i__] = -t19 + t50 + rho * (t56 + t77 + (sigma * 
			-.3178251739991049 * t80 * t45 + t24 * 
			.1362107888567592 * (t89 - t91) * t44 - t24 * 
			.1362107888567592 * t98 * (t89 - t91 + 
			1.19350059339396 / t33 / t27 * t35 * t40 * t85 * t26 
			- t36 * .08658243533790967 / t38 / t37 / rho)) * 
			.0310906908696549 * t116);
		vsigmaaa[i__] = rho * .1243627634786196 * (t23 * 
			.1362107888567592 * t32 * t44 + sigma * 
			.01855337900098064 * t40 * t28 * t44 - t24 * 
			.1362107888567592 * t98 * (t28 * .1362107888567592 * 
			t23 + t84 * .03710675800196129 * t40)) * t116;
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
		t2 = 1 / rho;
		t3 = pow_dd(&t2, &c_b2);
		t5 = t3 * .1325688999052018 + 1.;
		t6 = pow_dd(&t2, &c_b3);
		t9 = sqrt(t2);
/* Computing 2nd power */
		d__1 = t3;
		t11 = d__1 * d__1;
		t13 = t6 * 5.98255043577108 + t3 * 2.225569421150687 + t9 * 
			.8004286349993634 + t11 * .1897004325747559;
		t16 = 16.08197949869254 / t13 + 1.;
		t17 = log(t16);
		t18 = t5 * t17;
		t19 = t18 * .0621814;
/* Computing 2nd power */
		d__1 = rho;
		t20 = d__1 * d__1;
		t21 = pow_dd(&rho, &c_b2);
		t23 = 1 / t21 / t20;
		t24 = sigma * t23;
		t26 = exp(t18 * 2.000000587336264);
		t27 = t26 - 1.;
		t28 = 1 / t27;
		t29 = t28 * sigma;
		t31 = t29 * .1362107888567592 * t23;
		t32 = t31 + 1.;
/* Computing 2nd power */
		d__1 = t27;
		t33 = d__1 * d__1;
		t34 = 1 / t33;
/* Computing 2nd power */
		d__1 = sigma;
		t35 = d__1 * d__1;
		t36 = t34 * t35;
/* Computing 2nd power */
		d__1 = t20;
		t37 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t21;
		t38 = d__1 * d__1;
		t40 = 1 / t38 / t37;
		t43 = t31 + 1. + t36 * .01855337900098064 * t40;
		t44 = 1 / t43;
		t45 = t32 * t44;
		t48 = t24 * .1362107888567592 * t45 + 1.;
		t49 = log(t48);
		t50 = t49 * .0310906908696549;
		zk[i__] = rho * (-t19 + t50);
		t52 = 1 / t11;
		t53 = 1 / t20;
		t54 = t52 * t53;
		t55 = t54 * t17;
		t56 = t55 * .002747773264188438;
/* Computing 2nd power */
		d__1 = t13;
		t57 = d__1 * d__1;
		t58 = 1 / t57;
		t59 = t5 * t58;
/* Computing 2nd power */
		d__1 = t6;
		t60 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t60;
		t61 = d__1 * d__1;
		t62 = t61 * t6;
		t63 = 1 / t62;
		t67 = 1 / t9;
		t70 = 1 / t3;
		t73 = t63 * -.99709173929518 * t53 - t54 * .7418564737168958 
			- t67 * .4002143174996817 * t53 - t70 * 
			.1264669550498372 * t53;
		t74 = 1 / t16;
		t76 = t59 * t73 * t74;
		t77 = t76 * 1.;
		t78 = t20 * rho;
		t80 = 1 / t21 / t78;
		t81 = sigma * t80;
		t84 = t34 * sigma;
		t85 = t56 + t77;
		t89 = t84 * 4.381079514373337 * t23 * t85 * t26;
		t91 = t29 * .3178251739991049 * t80;
		t92 = t89 - t91;
		t93 = t92 * t44;
/* Computing 2nd power */
		d__1 = t43;
		t96 = d__1 * d__1;
		t97 = 1 / t96;
		t98 = t32 * t97;
		t100 = 1 / t33 / t27;
		t101 = t100 * t35;
		t103 = t40 * t85 * t26;
		t108 = 1 / t38 / t37 / rho;
		t111 = t89 - t91 + t101 * 1.19350059339396 * t103 - t36 * 
			.08658243533790967 * t108;
		t112 = t98 * t111;
		t115 = t81 * -.3178251739991049 * t45 + t24 * 
			.1362107888567592 * t93 - t24 * .1362107888567592 * 
			t112;
		t116 = 1 / t48;
		t117 = t115 * t116;
		vrhoa[i__] = -t19 + t50 + rho * (t56 + t77 + t117 * 
			.0310906908696549);
		t121 = t23 * t32;
		t124 = sigma * t40;
		t125 = t28 * t44;
		t132 = t28 * .1362107888567592 * t23 + t84 * 
			.03710675800196129 * t40;
		t133 = t98 * t132;
		t136 = t121 * .1362107888567592 * t44 + t124 * 
			.01855337900098064 * t125 - t24 * .1362107888567592 * 
			t133;
		vsigmaaa[i__] = rho * .1243627634786196 * t136 * t116;
		t151 = log(29.60874997779344 / (t6 * 8.157414703487641 + t3 * 
			2.247591863577616 + t9 * .4300972471276643 + t11 * 
			.1911512595127338) + 1.);
		t153 = (t3 * .06901399211255825 + 1.) * t151 * t53;
		t155 = t49 * t53;
		t158 = 1 / t21 / t37;
		t160 = sigma * t158 * t45;
		t165 = t153 * -1.086299437397317 + t18 * 1.333333724890842 * 
			t53;
		t169 = t84 * .1362107888567592 * t23 * t165 * t26;
		t170 = t29 * t158;
		t171 = t170 * .06053812838078188;
		t183 = t36 / t38 / t37 / t20;
		t197 = -t165;
		t201 = t84 * .1362107888567592 * t23 * t197 * t26;
		t202 = t170 * .06053812838078188;
		t223 = 1 / t37;
		t226 = 1 / t78;
		t231 = 1 / t11 / t2 * t223;
		t233 = t52 * t226;
		t250 = t59 * 1. * (-.8309097827459833 / t62 / t2 * t223 + t63 
			* 1.99418347859036 * t226 - t231 * .4945709824779306 
			+ t233 * 1.483712947433792 - .2001071587498409 / t9 / 
			t2 * t223 + t67 * .8004286349993634 * t226 - 
			.04215565168327908 / t3 / t2 * t223 + t70 * 
			.2529339100996745 * t226) * t74;
/* Computing 2nd power */
		d__1 = t57;
		t251 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t73;
		t254 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t16;
		t255 = d__1 * d__1;
		t259 = t5 * 16.08197949869254 / t251 * t254 / t255;
		t261 = t231 * .001831848842792292 * t17;
		t265 = t54 * .08837926660346786 * t58 * t73 * t74;
		t271 = t5 * 2. / t57 / t13 * t254 * t74;
		t273 = t233 * .005495546528376876 * t17;
/* Computing 2nd power */
		d__1 = t115;
		t274 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t48;
		t275 = d__1 * d__1;
		t276 = 1 / t275;
		t282 = t92 * t97;
		t287 = 1 / t96 / t43;
		t288 = t32 * t287;
/* Computing 2nd power */
		d__1 = t111;
		t289 = d__1 * d__1;
		t293 = t250 + t259 + t261 - t265 - t271 - t273;
		t297 = t84 * 4.381079514373337 * t23 * t293 * t26;
/* Computing 2nd power */
		d__1 = t85;
		t298 = d__1 * d__1;
		t299 = t23 * t298;
		t302 = t84 * 140.9129032462046 * t299 * t26;
/* Computing 2nd power */
		d__1 = t33;
		t311 = d__1 * d__1;
		t314 = t40 * t298;
/* Computing 2nd power */
		d__1 = t26;
		t315 = d__1 * d__1;
		t325 = t84 * 20.44503773374224 * t80 * t85 * t26;
		t326 = t100 * sigma;
		t329 = t326 * 281.8258064924092 * t299 * t315;
		t330 = t170 * 1.059417246663683;
		s1 = t55 * .01099109305675375 + t76 * 4. + rho * (t153 * 
			.0337738 - t155 * .0207271272464366 + (t160 * 
			.06053812838078188 + t24 * .1362107888567592 * (-t169 
			+ t171) * t44 - t24 * .1362107888567592 * t98 * (
			-t169 + t171 - t101 * .03710675800196129 * t40 * t165 
			* t26 + t183 * .01649189244531613)) * 
			.0310906908696549 * t116);
		s2 = s1 + rho * (t153 * -.0337738 + t155 * .0207271272464366 
			+ (t160 * -.06053812838078188 + t24 * 
			.1362107888567592 * (-t201 - t202) * t44 - t24 * 
			.1362107888567592 * t98 * (-t201 - t202 - t101 * 
			.03710675800196129 * t40 * t197 * t26 - t183 * 
			.01649189244531613)) * .0310906908696549 * t116);
		v2rhoa2[i__] = s2 + rho * 2. * (t250 + t259 + t261 - t265 - 
			t271 - t273 - t274 * .0310906908696549 * t276 + (t160 
			* 1.059417246663683 - t81 * .6356503479982097 * t93 - 
			t24 * .2724215777135184 * t282 * t111 + t24 * 
			.2724215777135184 * t288 * t289 - t24 * 
			.1362107888567592 * t98 * (t297 - t302 - t101 * 
			11.13933887167696 * t108 * t85 * t26 + t101 * 
			1.19350059339396 * t40 * t293 * t26 + 
			115.1631462675703 / t311 * t35 * t314 * t315 - t101 * 
			38.38771542252344 * t314 * t26 - t325 + t329 + t330 + 
			t183 * .4906338002481548) + t81 * .6356503479982097 * 
			t112 + t24 * .1362107888567592 * (t329 - t325 + t297 
			- t302 + t330) * t44) * .0310906908696549 * t116) + 
			t117 * .1243627634786196;
		t365 = t34 * 4.381079514373337 * t23 * t85 * t26;
		t367 = t28 * .3178251739991049 * t80;
		t378 = t28 * t97;
		v2rhoasigmaaa[i__] = t136 * .1243627634786196 * t116 + rho * 
			4. * ((t80 * -.3178251739991049 * t32 * t44 - sigma * 
			.04329121766895483 * t108 * t125 + t81 * 
			.3178251739991049 * t133 + t23 * .1362107888567592 * 
			t92 * t44 + t24 * .1362107888567592 * (t365 - t367) * 
			t44 - t24 * .1362107888567592 * t282 * t132 - t121 * 
			.1362107888567592 * t97 * t111 - t124 * 
			.01855337900098064 * t378 * t111 + t24 * 
			.2724215777135184 * t32 * t287 * t111 * t132 - t24 * 
			.1362107888567592 * t98 * (t365 - t367 + t326 * 
			2.38700118678792 * t103 - t84 * .1731648706758193 * 
			t108)) * .0310906908696549 * t116 - t115 * 
			.0310906908696549 * t276 * t136);
/* Computing 2nd power */
		d__1 = t132;
		t413 = d__1 * d__1;
/* Computing 2nd power */
		d__1 = t136;
		t427 = d__1 * d__1;
		v2sigmaaa2[i__] = rho * .4974510539144783 * (t40 * 
			.03710675800196129 * t28 * t44 - t121 * 
			.2724215777135184 * t97 * t132 - t124 * 
			.03710675800196129 * t378 * t132 + t24 * 
			.2724215777135184 * t288 * t413 - sigma * 
			.005054340779364009 / t37 / t78 * t98 * t34) * t116 - 
			rho * .4974510539144783 * t427 * t276;
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
} /* rks_c_pbe__ */


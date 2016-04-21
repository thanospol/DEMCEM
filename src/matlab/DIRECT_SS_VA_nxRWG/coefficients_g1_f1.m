function [coef] = coefficients_g1_f1(r1,r2,r3,r4,r5, AreaT)
%% coefficients_g1_f1

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  20 September 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,” Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,” IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,” IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% r1,r2,r3,r4,r5 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r1,r4,r5)
% AreaT = area of outer integral

% OUTPUT DATA
% coef   
%

coef  = zeros(6);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = sqrt(0.3e1);
t2 = 0.1e1 / AreaT;
t3 = t1 * t2;
t4 = bb(3) ^ 2;
t5 = t3 * t4;
t6 = cc(2) * bb(1);
t7 = t6 * dd(3);
t8 = t5 * t7;
t10 = bb(2) ^ 2;
t11 = t3 * t10;
t12 = t11 * t7;
t14 = bb(1) ^ 2;
t15 = t3 * t14;
t16 = cc(2) * bb(3);
t17 = t16 * dd(1);
t18 = t15 * t17;
t20 = cc(1) * bb(3);
t21 = t20 * dd(2);
t22 = t15 * t21;
t24 = cc(3) * bb(2);
t25 = t24 * dd(1);
t26 = t15 * t25;
t28 = cc(1) * bb(2);
t29 = t28 * dd(3);
t30 = t5 * t29;
t32 = cc(3) * bb(1);
t33 = t32 * dd(2);
t34 = t5 * t33;
t36 = t15 * t29;
t38 = t11 * t17;
t40 = t5 * t25;
t42 = t11 * t33;
t44 = t11 * t21;
t46 = cc(3) * cc(1);
t47 = t46 * dd(2);
t48 = t5 * t47;
t50 = t3 * bb(2);
t51 = cc(2) ^ 2;
t52 = t51 * bb(3);
t53 = t52 * dd(1);
t54 = t50 * t53;
t56 = t3 * bb(1);
t57 = cc(1) ^ 2;
t58 = t57 * bb(2);
t59 = t58 * dd(3);
t60 = t56 * t59;
t62 = cc(1) * cc(2);
t63 = t62 * dd(3);
t64 = t11 * t63;
t66 = t3 * bb(3);
t67 = cc(3) ^ 2;
t68 = t67 * bb(1);
t69 = t68 * dd(2);
t70 = t66 * t69;
t72 = t15 * t63;
t74 = t8 / 0.24e2 + t12 / 0.24e2 - t18 / 0.24e2 + t22 / 0.24e2 + t26 / 0.24e2 - t30 / 0.24e2 - t34 / 0.24e2 - t36 / 0.24e2 - t38 / 0.24e2 + t40 / 0.24e2 - t42 / 0.24e2 + t44 / 0.24e2 - t48 / 0.12e2 + t54 / 0.12e2 + t60 / 0.12e2 + t64 / 0.12e2 + t70 / 0.12e2 - t72 / 0.12e2;
t75 = t51 * bb(1);
t76 = t75 * dd(3);
t77 = t50 * t76;
t79 = cc(2) * cc(3);
t80 = t79 * dd(1);
t81 = t5 * t80;
t83 = t67 * bb(2);
t84 = t83 * dd(1);
t85 = t66 * t84;
t87 = t11 * t80;
t89 = t15 * t47;
t91 = t57 * bb(3);
t92 = t91 * dd(2);
t93 = t56 * t92;
t95 = t10 * bb(2);
t96 = t95 * cc(3);
t98 = t3 * t96 * dd(1);
t100 = t14 * bb(1);
t101 = t100 * cc(2);
t103 = t3 * t101 * dd(3);
t105 = t4 * bb(3);
t106 = t105 * cc(2);
t108 = t3 * t106 * dd(1);
t110 = t95 * cc(1);
t112 = t3 * t110 * dd(3);
t114 = t105 * cc(1);
t116 = t3 * t114 * dd(2);
t118 = t100 * cc(3);
t120 = t3 * t118 * dd(2);
t122 = t3 * cc(2);
t123 = bb(3) * dd(2);
t125 = t122 * t28 * t123;
t127 = cc(2) * dd(3);
t129 = t66 * t32 * t127;
t131 = cc(3) * dd(2);
t133 = t50 * t6 * t131;
t135 = t3 * cc(3);
t136 = bb(2) * dd(3);
t138 = t135 * t20 * t136;
t140 = t3 * cc(1);
t141 = bb(3) * dd(1);
t143 = t140 * t6 * t141;
t145 = cc(3) * dd(1);
t147 = t56 * t28 * t145;
t149 = -t77 / 0.12e2 + t81 / 0.12e2 - t85 / 0.12e2 - t87 / 0.12e2 + t89 / 0.12e2 - t93 / 0.12e2 + t98 / 0.24e2 + t103 / 0.24e2 - t108 / 0.24e2 - t112 / 0.24e2 + t116 / 0.24e2 - t120 / 0.24e2 - t125 / 0.12e2 - t129 / 0.12e2 + t133 / 0.12e2 + t138 / 0.12e2 + t143 / 0.12e2 - t147 / 0.12e2;
t151 = t2 * bb(3);
t152 = t151 * cc(3);
t154 = t152 * t7 / 0.12e2;
t155 = t6 * ee(3);
t158 = t2 * bb(2);
t159 = t158 * cc(2);
t161 = t159 * t33 / 0.12e2;
t162 = t32 * ee(2);
t166 = t2 * cc(3) * cc(1);
t167 = bb(3) * bb(2);
t170 = t166 * t167 * dd(3) / 0.12e2;
t175 = t2 * cc(2) * cc(1);
t178 = t175 * t167 * dd(2) / 0.12e2;
t182 = t2 * bb(1);
t183 = t182 * cc(1);
t185 = t183 * t25 / 0.12e2;
t186 = t24 * ee(1);
t189 = bb(1) * bb(3);
t192 = t175 * t189 * dd(1) / 0.12e2;
t196 = t2 * t4;
t197 = t28 * ee(3);
t200 = t2 * t10;
t201 = t200 * t21;
t203 = t20 * ee(2);
t206 = t200 * t7;
t210 = t2 * t14;
t211 = t210 * t29;
t213 = t154 - t152 * t155 / 0.6e1 - t161 + t159 * t162 / 0.6e1 - t170 + t166 * t167 * ee(3) / 0.6e1 + t178 - t175 * t167 * ee(2) / 0.6e1 + t185 - t183 * t186 / 0.6e1 - t192 + t175 * t189 * ee(1) / 0.6e1 - t196 * t197 / 0.12e2 - t201 / 0.24e2 + t200 * t203 / 0.12e2 - t206 / 0.24e2 + t200 * t155 / 0.12e2 + t211 / 0.24e2;
t216 = t210 * t17;
t218 = t16 * ee(1);
t221 = t196 * t7;
t225 = t196 * t25;
t229 = t200 * t17;
t233 = t200 * t33;
t237 = t210 * t25;
t241 = t210 * t21;
t245 = t196 * t29;
t249 = t196 * t33;
t251 = -t210 * t197 / 0.12e2 + t216 / 0.24e2 - t210 * t218 / 0.12e2 - t221 / 0.24e2 + t196 * t155 / 0.12e2 - t225 / 0.24e2 + t196 * t186 / 0.12e2 + t229 / 0.24e2 - t200 * t218 / 0.12e2 + t233 / 0.24e2 - t200 * t162 / 0.12e2 - t237 / 0.24e2 + t210 * t186 / 0.12e2 - t241 / 0.24e2 + t210 * t203 / 0.12e2 + t245 / 0.24e2 - t196 * t162 / 0.12e2 + t249 / 0.24e2;
t254 = t196 * t47 / 0.12e2;
t255 = t62 * ee(3);
t258 = t46 * ee(2);
t261 = t68 * ee(2);
t265 = t200 * t63 / 0.12e2;
t267 = t151 * t69 / 0.12e2;
t269 = t158 * t76 / 0.12e2;
t270 = t75 * ee(3);
t274 = t182 * t59 / 0.12e2;
t275 = t58 * ee(3);
t279 = t210 * t63 / 0.12e2;
t283 = t196 * t80 / 0.12e2;
t284 = t79 * ee(1);
t288 = t151 * t84 / 0.12e2;
t289 = t83 * ee(1);
t293 = t158 * t53 / 0.12e2;
t294 = t52 * ee(1);
t297 = t254 + t200 * t255 / 0.6e1 - t196 * t258 / 0.6e1 + t151 * t261 / 0.6e1 - t265 - t267 + t269 - t158 * t270 / 0.6e1 - t274 + t182 * t275 / 0.6e1 + t279 - t210 * t255 / 0.6e1 - t283 + t196 * t284 / 0.6e1 + t288 - t151 * t289 / 0.6e1 - t293 + t158 * t294 / 0.6e1;
t298 = t200 * t80;
t302 = t210 * t47;
t306 = t182 * t92;
t308 = t91 * ee(2);
t311 = t2 * t95;
t315 = t2 * t100;
t316 = t315 * t127;
t318 = cc(2) * ee(3);
t321 = t2 * t105;
t323 = t321 * cc(2) * dd(1);
t328 = t311 * t145;
t330 = cc(3) * ee(1);
t333 = t315 * t131;
t335 = cc(3) * ee(2);
t342 = t311 * cc(1) * dd(3);
t345 = t321 * cc(1) * dd(2);
t347 = t298 / 0.12e2 - t200 * t284 / 0.6e1 - t302 / 0.12e2 + t210 * t258 / 0.6e1 + t306 / 0.12e2 - t182 * t308 / 0.6e1 - t311 * cc(1) * ee(3) / 0.12e2 - t316 / 0.24e2 + t315 * t318 / 0.12e2 + t323 / 0.24e2 - t321 * cc(2) * ee(1) / 0.12e2 - t328 / 0.24e2 + t311 * t330 / 0.12e2 + t333 / 0.24e2 - t315 * t335 / 0.12e2 + t321 * cc(1) * ee(2) / 0.12e2 + t342 / 0.24e2 - t345 / 0.24e2;
t350 = -t221 + t211 - t206 + t249 + t216 + t333 + t233 + t342 - t237 - t328 + t245 - t225 - t316 + t323 - t345 + t229 - t241 - t201;
t351 = t5 * t155;
t353 = t5 * t186;
t355 = t15 * t186;
t357 = t15 * t203;
t359 = t5 * t162;
t361 = t15 * t218;
t363 = t5 * t197;
t365 = t15 * t197;
t367 = t11 * t203;
t369 = t11 * t155;
t371 = t11 * t162;
t373 = t11 * t218;
t381 = -t351 / 0.24e2 - t353 / 0.24e2 - t355 / 0.24e2 - t357 / 0.24e2 + t359 / 0.24e2 + t361 / 0.24e2 + t363 / 0.24e2 + t365 / 0.24e2 - t367 / 0.24e2 - t369 / 0.24e2 + t371 / 0.24e2 + t373 / 0.24e2 + t8 / 0.48e2 + t12 / 0.48e2 - t18 / 0.48e2 + t22 / 0.48e2 + t26 / 0.48e2 - t30 / 0.48e2;
t395 = t3 * t110 * ee(3);
t398 = t3 * t96 * ee(1);
t401 = t3 * t118 * ee(2);
t404 = t3 * t101 * ee(3);
t407 = t3 * t106 * ee(1);
t410 = t3 * t114 * ee(2);
t412 = -t34 / 0.48e2 - t36 / 0.48e2 - t38 / 0.48e2 + t40 / 0.48e2 - t42 / 0.48e2 + t44 / 0.48e2 + t98 / 0.48e2 + t103 / 0.48e2 - t108 / 0.48e2 - t112 / 0.48e2 + t116 / 0.48e2 - t120 / 0.48e2 + t395 / 0.24e2 - t398 / 0.24e2 + t401 / 0.24e2 - t404 / 0.24e2 + t407 / 0.24e2 - t410 / 0.24e2;
t421 = t154 - t161 - t170 + t178 + t185 - t192 - t201 / 0.48e2 - t206 / 0.48e2 + t211 / 0.48e2 + t216 / 0.48e2 - t221 / 0.48e2 - t225 / 0.48e2 + t229 / 0.48e2;
t427 = t233 / 0.48e2 - t237 / 0.48e2 - t241 / 0.48e2 + t245 / 0.48e2 + t249 / 0.48e2 + t254 - t265 - t267 + t269 - t274 + t279 - t283 + t288 - t293;
t429 = t2 * t67;
t432 = t2 * t51;
t435 = t2 * t57;
t442 = t298 - t302 + t306 - t429 * t7 - t429 * t21 + t432 * t29 + t432 * t33 - t435 * t25 - t435 * t7 + t429 * t17 + t429 * t29 - t432 * t21 - t432 * t25;
t453 = t67 * cc(3);
t454 = t2 * t453;
t458 = t51 * cc(2);
t459 = t2 * t458;
t463 = t57 * cc(1);
t464 = t2 * t463;
t474 = t435 * t33 / 0.12e2 + t435 * t17 / 0.12e2 - t316 / 0.48e2 + t323 / 0.48e2 - t328 / 0.48e2 + t333 / 0.48e2 + t342 / 0.48e2 - t345 / 0.48e2 + t454 * bb(1) * dd(2) / 0.12e2 - t459 * bb(1) * dd(3) / 0.12e2 + t464 * t136 / 0.12e2 - t454 * bb(2) * dd(1) / 0.12e2 + t459 * t141 / 0.12e2 - t464 * t123 / 0.12e2;
t492 = t3 * t57;
t495 = t3 * t51;
t498 = t3 * t67;
t502 = t66 * t32 * t318 - t50 * t6 * t335 + t56 * t28 * t330 - t140 * t6 * bb(3) * ee(1) + t122 * t28 * bb(3) * ee(2) - t135 * t20 * bb(2) * ee(3) - t492 * t155 + t492 * t162 - t495 * t203 + t492 * t218 - t498 * t155 - t495 * t186 + t495 * t197;
t519 = t495 * t162 / 0.18e2 - t492 * t186 / 0.18e2 - t351 / 0.72e2 - t353 / 0.72e2 - t355 / 0.72e2 - t357 / 0.72e2 + t359 / 0.72e2 + t361 / 0.72e2 + t363 / 0.72e2 + t365 / 0.72e2 - t367 / 0.72e2 - t369 / 0.72e2 + t371 / 0.72e2 + t373 / 0.72e2;
t547 = -t66 * t261 / 0.18e2 + t50 * t270 / 0.18e2 + t5 * t258 / 0.18e2 - t11 * t255 / 0.18e2 - t5 * t284 / 0.18e2 + t66 * t289 / 0.18e2 - t50 * t294 / 0.18e2 + t11 * t284 / 0.18e2 - t56 * t275 / 0.18e2 - t15 * t258 / 0.18e2 + t15 * t255 / 0.18e2 + t56 * t308 / 0.18e2 + t495 * t25 / 0.36e2;
t576 = t492 * t25 / 0.36e2 + t495 * t21 / 0.36e2 - t498 * t17 / 0.36e2 + t492 * t7 / 0.36e2 - t492 * t33 / 0.36e2 - t498 * t29 / 0.36e2 + t498 * t21 / 0.36e2 - t492 * t17 / 0.36e2 - t495 * t29 / 0.36e2 - t495 * t33 / 0.36e2 + t498 * t7 / 0.36e2 + t498 * t197 / 0.18e2 - t498 * t203 / 0.18e2 + t498 * t218 / 0.18e2;
t592 = t8 / 0.144e3 + t12 / 0.144e3 - t18 / 0.144e3 + t22 / 0.144e3 + t26 / 0.144e3 - t30 / 0.144e3 - t34 / 0.144e3 - t36 / 0.144e3 - t38 / 0.144e3 + t40 / 0.144e3 - t42 / 0.144e3 + t44 / 0.144e3 - t48 / 0.36e2;
t607 = t54 / 0.36e2 + t60 / 0.36e2 + t64 / 0.36e2 + t70 / 0.36e2 - t72 / 0.36e2 - t77 / 0.36e2 + t81 / 0.36e2 - t85 / 0.36e2 - t87 / 0.36e2 + t89 / 0.36e2 - t93 / 0.36e2 + t98 / 0.144e3 + t103 / 0.144e3 - t108 / 0.144e3;
t618 = t463 * bb(3);
t622 = t463 * bb(2);
t626 = t458 * bb(3);
t630 = t453 * bb(1);
t634 = -t112 / 0.144e3 + t116 / 0.144e3 - t120 / 0.144e3 + t395 / 0.72e2 - t398 / 0.72e2 + t401 / 0.72e2 - t404 / 0.72e2 + t407 / 0.72e2 - t410 / 0.72e2 + t3 * t618 * dd(2) / 0.36e2 - t3 * t622 * dd(3) / 0.36e2 - t3 * t626 * dd(1) / 0.36e2 - t3 * t630 * dd(2) / 0.36e2;
t635 = t453 * bb(2);
t639 = t458 * bb(1);
t667 = t3 * t635 * dd(1) / 0.36e2 + t3 * t639 * dd(3) / 0.36e2 + t3 * t622 * ee(3) / 0.18e2 + t3 * t626 * ee(1) / 0.18e2 + t3 * t630 * ee(2) / 0.18e2 - t3 * t618 * ee(2) / 0.18e2 - t3 * t639 * ee(3) / 0.18e2 - t3 * t635 * ee(1) / 0.18e2 - t125 / 0.36e2 - t129 / 0.36e2 + t133 / 0.36e2 + t138 / 0.36e2 + t143 / 0.36e2 - t147 / 0.36e2;
%
coef(1) = t74 + t149;
coef(2) = t213 + t251 + t297 + t347;
coef(3) = t350 / 0.16e2;
coef(4) = t381 + t412;
coef(5) = t421 + t427 + t442 / 0.12e2 + t474;
coef(6) = t502 / 0.18e2 + t519 + t547 + t576 + t592 + t607 + t634 + t667;
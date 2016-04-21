function [coef] = coefficients_g2_f2(r1,r2,r3,r4,r5, AreaT)
%% coefficients_g2_f2

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  20 September 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, �On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,� IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, �Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,� Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, �Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,� IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, �On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,� IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% r1,r2,r3,r4,r5 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r1,r4,r5)
% AreaT = area of outer integral

% OUTPUT DATA
% coef   
%

coef  = zeros(16);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = 0.1e1 / AreaT;
t2 = bb(3) ^ 2;
t3 = t1 * t2;
t4 = cc(3) * bb(1);
t5 = t4 * dd(2);
t6 = t3 * t5;
t7 = bb(2) ^ 2;
t8 = t1 * t7;
t9 = t8 * t5;
t10 = t7 * bb(2);
t11 = t1 * t10;
t12 = cc(1) * dd(3);
t13 = t11 * t12;
t14 = cc(3) * dd(1);
t15 = t11 * t14;
t16 = cc(2) * bb(1);
t17 = t16 * dd(3);
t18 = t8 * t17;
t19 = t3 * t17;
t20 = bb(1) ^ 2;
t21 = t1 * t20;
t22 = cc(3) * bb(2);
t23 = t22 * dd(1);
t24 = t21 * t23;
t25 = cc(1) * bb(2);
t26 = t25 * dd(3);
t27 = t21 * t26;
t28 = t3 * t26;
t29 = cc(1) * bb(3);
t30 = t29 * dd(2);
t31 = t21 * t30;
t32 = t20 * bb(1);
t33 = t1 * t32;
t34 = cc(2) * dd(3);
t35 = t33 * t34;
t36 = t2 * bb(3);
t37 = t1 * t36;
t38 = cc(2) * dd(1);
t39 = t37 * t38;
t40 = cc(2) * bb(3);
t41 = t40 * dd(1);
t42 = t21 * t41;
t43 = cc(1) * dd(2);
t44 = t37 * t43;
t45 = cc(3) * dd(2);
t46 = t33 * t45;
t47 = t3 * t23;
t48 = t8 * t41;
t49 = t8 * t30;
t50 = t6 + t9 + t13 - t15 - t18 - t19 - t24 + t27 + t28 - t31 - t35 + t39 + t42 - t44 + t46 - t47 + t48 - t49;
t51 = sqrt(0.3e1);
t52 = t43 * t51;
t53 = t37 * t52;
t55 = t12 * t51;
t56 = t11 * t55;
t58 = t45 * t51;
t59 = t33 * t58;
t61 = t14 * t51;
t62 = t11 * t61;
t64 = t34 * t51;
t65 = t33 * t64;
t67 = t38 * t51;
t68 = t37 * t67;
t70 = t51 * t1;
t73 = t70 * t36 * cc(2) * ee(1);
t77 = t70 * t32 * cc(3) * ee(2);
t81 = t70 * t10 * cc(1) * ee(3);
t85 = t70 * t32 * cc(2) * ee(3);
t89 = t70 * t10 * cc(3) * ee(1);
t93 = t70 * t36 * cc(1) * ee(2);
t95 = t21 * cc(2);
t96 = bb(3) * dd(1);
t97 = t96 * t51;
t98 = t95 * t97;
t100 = t21 * cc(1);
t101 = bb(3) * dd(2);
t102 = t101 * t51;
t103 = t100 * t102;
t105 = t3 * cc(3);
t106 = bb(1) * dd(2);
t107 = t106 * t51;
t108 = t105 * t107;
t110 = t8 * cc(2);
t111 = bb(1) * dd(3);
t112 = t111 * t51;
t113 = t110 * t112;
t115 = t3 * cc(2);
t116 = t115 * t112;
t118 = t21 * cc(3);
t119 = bb(2) * dd(1);
t120 = t119 * t51;
t121 = t118 * t120;
t123 = t53 / 0.48e2 - t56 / 0.48e2 - t59 / 0.48e2 + t62 / 0.48e2 + t65 / 0.48e2 - t68 / 0.48e2 + t73 / 0.24e2 + t77 / 0.24e2 + t81 / 0.24e2 - t85 / 0.24e2 - t89 / 0.24e2 - t93 / 0.24e2 - t98 / 0.48e2 + t103 / 0.48e2 - t108 / 0.48e2 + t113 / 0.48e2 + t116 / 0.48e2 + t121 / 0.48e2;
t124 = t3 * cc(1);
t125 = bb(2) * dd(3);
t126 = t125 * t51;
t127 = t124 * t126;
t129 = t105 * t120;
t131 = t100 * t126;
t133 = t8 * cc(1);
t134 = t133 * t102;
t136 = t110 * t97;
t138 = t8 * cc(3);
t139 = t138 * t107;
t141 = t70 * t7;
t142 = t16 * ee(3);
t143 = t141 * t142;
t145 = t70 * t20;
t146 = t22 * ee(1);
t147 = t145 * t146;
t149 = t4 * ee(2);
t150 = t141 * t149;
t152 = t70 * t2;
t153 = t152 * t142;
t155 = t152 * t146;
t157 = t29 * ee(2);
t158 = t141 * t157;
t160 = t25 * ee(3);
t161 = t145 * t160;
t163 = t152 * t149;
t165 = t145 * t157;
t167 = t40 * ee(1);
t168 = t145 * t167;
t170 = t141 * t167;
t172 = t152 * t160;
t174 = -t127 / 0.48e2 + t129 / 0.48e2 - t131 / 0.48e2 + t134 / 0.48e2 - t136 / 0.48e2 - t139 / 0.48e2 - t143 / 0.24e2 - t147 / 0.24e2 + t150 / 0.24e2 - t153 / 0.24e2 - t155 / 0.24e2 - t158 / 0.24e2 + t161 / 0.24e2 + t163 / 0.24e2 - t165 / 0.24e2 + t168 / 0.24e2 + t170 / 0.24e2 + t172 / 0.24e2;
t176 = cc(2) ^ 2;
t177 = t1 * t176;
t178 = t177 * t26;
t180 = cc(3) ^ 2;
t181 = t1 * t180;
t182 = t181 * t30;
t195 = t178 / 0.12e2 - t182 / 0.12e2 + t6 / 0.48e2 + t9 / 0.48e2 - t18 / 0.48e2 - t19 / 0.48e2 - t24 / 0.48e2 + t27 / 0.48e2 + t28 / 0.48e2 - t31 / 0.48e2 + t42 / 0.48e2 - t47 / 0.48e2 + t48 / 0.48e2;
t197 = t177 * t30;
t199 = t177 * t23;
t201 = cc(1) ^ 2;
t202 = t1 * t201;
t203 = t202 * t5;
t205 = t202 * t41;
t207 = t1 * bb(2);
t208 = t176 * bb(3);
t210 = t207 * t208 * dd(1);
t211 = t210 / 0.12e2;
t212 = cc(3) * cc(2);
t213 = t212 * dd(1);
t214 = t8 * t213;
t215 = t214 / 0.12e2;
t216 = cc(3) * cc(1);
t217 = t216 * dd(2);
t218 = t21 * t217;
t219 = t218 / 0.12e2;
t220 = t1 * bb(1);
t221 = t201 * bb(3);
t223 = t220 * t221 * dd(2);
t224 = t223 / 0.12e2;
t225 = t1 * bb(3);
t226 = t180 * bb(1);
t227 = t226 * dd(2);
t228 = t225 * t227;
t229 = t228 / 0.12e2;
t230 = t3 * t217;
t231 = t230 / 0.12e2;
t232 = cc(1) * cc(2);
t233 = t232 * dd(3);
t234 = t8 * t233;
t235 = t234 / 0.12e2;
t236 = t176 * bb(1);
t238 = t207 * t236 * dd(3);
t239 = t238 / 0.12e2;
t240 = t181 * t17;
t242 = -t49 / 0.48e2 - t197 / 0.12e2 - t199 / 0.12e2 + t203 / 0.12e2 + t205 / 0.12e2 - t211 + t215 - t219 + t224 - t229 + t231 - t235 + t239 - t240 / 0.12e2;
t244 = t225 * cc(3);
t245 = t244 * t17;
t246 = t207 * cc(2);
t247 = t246 * t5;
t248 = t220 * cc(1);
t249 = t248 * t23;
t250 = t1 * cc(3);
t251 = t250 * cc(1);
t252 = bb(3) * bb(2);
t254 = t251 * t252 * dd(3);
t256 = t1 * cc(2) * cc(1);
t258 = t256 * t252 * dd(2);
t259 = bb(1) * bb(3);
t261 = t256 * t259 * dd(1);
t263 = t201 * bb(2);
t265 = t220 * t263 * dd(3);
t266 = t21 * t233;
t267 = t3 * t213;
t268 = t180 * bb(2);
t270 = t225 * t268 * dd(1);
t273 = t245 - t247 + t249 - t254 + t258 - t261 + t177 * t5 - t265 + t266 - t267 + t270 - t202 * t23 - t202 * t17;
t274 = t181 * t41;
t276 = t181 * t26;
t284 = t180 * cc(3);
t285 = t1 * t284;
t286 = t285 * t106;
t288 = t176 * cc(2);
t289 = t1 * t288;
t290 = t289 * t111;
t292 = t201 * cc(1);
t293 = t1 * t292;
t294 = t293 * t125;
t296 = t285 * t119;
t298 = t289 * t96;
t300 = t293 * t101;
t302 = t274 / 0.12e2 + t276 / 0.12e2 + t13 / 0.48e2 - t15 / 0.48e2 - t35 / 0.48e2 + t39 / 0.48e2 - t44 / 0.48e2 + t46 / 0.48e2 + t286 / 0.12e2 - t290 / 0.12e2 + t294 / 0.12e2 - t296 / 0.12e2 + t298 / 0.12e2 - t300 / 0.12e2;
t317 = dd(1) * t51;
t319 = t256 * t259 * t317;
t321 = t53 / 0.144e3 - t56 / 0.144e3 - t59 / 0.144e3 + t62 / 0.144e3 + t65 / 0.144e3 - t68 / 0.144e3 + t73 / 0.72e2 + t77 / 0.72e2 + t81 / 0.72e2 - t85 / 0.72e2 - t89 / 0.72e2 - t93 / 0.72e2 + t319 / 0.36e2;
t322 = dd(2) * t51;
t324 = t256 * t252 * t322;
t327 = t248 * t22 * t317;
t330 = t246 * t4 * t322;
t332 = dd(3) * t51;
t334 = t251 * t252 * t332;
t337 = t244 * t16 * t332;
t339 = t70 * bb(1);
t340 = cc(3) * ee(1);
t345 = bb(2) * ee(3);
t350 = bb(3) * ee(2);
t355 = bb(3) * ee(1);
t359 = t70 * bb(3);
t360 = cc(2) * ee(3);
t364 = t70 * bb(2);
t365 = cc(3) * ee(2);
t372 = -t324 / 0.36e2 - t327 / 0.36e2 + t330 / 0.36e2 + t334 / 0.36e2 - t337 / 0.36e2 + t339 * t25 * t340 / 0.18e2 - t70 * cc(3) * t29 * t345 / 0.18e2 + t70 * cc(2) * t25 * t350 / 0.18e2 - t70 * cc(1) * t16 * t355 / 0.18e2 + t359 * t4 * t360 / 0.18e2 - t364 * t16 * t365 / 0.18e2 - t98 / 0.144e3 + t103 / 0.144e3 - t108 / 0.144e3;
t387 = t113 / 0.144e3 + t116 / 0.144e3 + t121 / 0.144e3 - t127 / 0.144e3 + t129 / 0.144e3 - t131 / 0.144e3 + t134 / 0.144e3 - t136 / 0.144e3 - t139 / 0.144e3 - t143 / 0.72e2 - t147 / 0.72e2 + t150 / 0.72e2 - t153 / 0.72e2;
t396 = t359 * t227;
t398 = t118 * t52;
t400 = t133 * t64;
t402 = t220 * t201;
t403 = t402 * t126;
t405 = t95 * t55;
t407 = t115 * t61;
t409 = -t155 / 0.72e2 - t158 / 0.72e2 + t161 / 0.72e2 + t163 / 0.72e2 - t165 / 0.72e2 + t168 / 0.72e2 + t170 / 0.72e2 + t172 / 0.72e2 + t396 / 0.36e2 + t398 / 0.36e2 + t400 / 0.36e2 + t403 / 0.36e2 - t405 / 0.36e2 + t407 / 0.36e2;
t412 = t138 * t67;
t415 = t225 * t180 * t120;
t417 = t402 * t102;
t419 = t207 * t176;
t420 = t419 * t97;
t422 = t419 * t112;
t424 = t124 * t58;
t426 = t263 * ee(3);
t429 = t212 * ee(1);
t432 = t216 * ee(2);
t435 = t226 * ee(2);
t438 = t221 * ee(2);
t443 = t232 * ee(3);
t446 = -t412 / 0.36e2 - t415 / 0.36e2 - t417 / 0.36e2 + t420 / 0.36e2 - t422 / 0.36e2 - t424 / 0.36e2 - t339 * t426 / 0.18e2 - t152 * t429 / 0.18e2 + t152 * t432 / 0.18e2 - t359 * t435 / 0.18e2 + t339 * t438 / 0.18e2 - t145 * t432 / 0.18e2 - t141 * t443 / 0.18e2;
t447 = t208 * ee(1);
t452 = t236 * ee(3);
t457 = t268 * ee(1);
t460 = t70 * t180;
t463 = t70 * t201;
t466 = t70 * t176;
t481 = -t364 * t447 / 0.18e2 + t145 * t443 / 0.18e2 + t364 * t452 / 0.18e2 + t141 * t429 / 0.18e2 + t359 * t457 / 0.18e2 + t460 * t30 / 0.36e2 + t463 * t23 / 0.36e2 + t466 * t23 / 0.36e2 - t460 * t26 / 0.36e2 - t466 * t5 / 0.36e2 + t466 * t30 / 0.36e2 + t463 * t17 / 0.36e2 - t466 * t26 / 0.36e2 - t460 * t41 / 0.36e2;
t509 = t460 * t17 / 0.36e2 - t463 * t5 / 0.36e2 - t463 * t41 / 0.36e2 - t463 * t142 / 0.18e2 + t460 * t160 / 0.18e2 - t460 * t157 / 0.18e2 + t466 * t149 / 0.18e2 - t466 * t146 / 0.18e2 - t466 * t157 / 0.18e2 - t460 * t142 / 0.18e2 + t463 * t149 / 0.18e2 + t460 * t167 / 0.18e2 + t466 * t160 / 0.18e2;
t514 = t284 * bb(2);
t518 = t292 * bb(3);
t522 = t288 * bb(3);
t526 = t292 * bb(2);
t530 = t284 * bb(1);
t534 = t288 * bb(1);
t556 = -t463 * t146 / 0.18e2 + t463 * t167 / 0.18e2 + t70 * t514 * dd(1) / 0.36e2 + t70 * t518 * dd(2) / 0.36e2 - t70 * t522 * dd(1) / 0.36e2 - t70 * t526 * dd(3) / 0.36e2 - t70 * t530 * dd(2) / 0.36e2 + t70 * t534 * dd(3) / 0.36e2 + t70 * t522 * ee(1) / 0.18e2 - t70 * t518 * ee(2) / 0.18e2 - t70 * t514 * ee(1) / 0.18e2 + t70 * t530 * ee(2) / 0.18e2 + t70 * t526 * ee(3) / 0.18e2 - t70 * t534 * ee(3) / 0.18e2;
t560 = t53 / 0.24e2;
t561 = t56 / 0.24e2;
t562 = t59 / 0.24e2;
t563 = t62 / 0.24e2;
t564 = t65 / 0.24e2;
t565 = t68 / 0.24e2;
t572 = t98 / 0.24e2;
t573 = t103 / 0.24e2;
t574 = t108 / 0.24e2;
t575 = t113 / 0.24e2;
t576 = t116 / 0.24e2;
t577 = t121 / 0.24e2;
t578 = -t560 + t561 + t562 - t563 - t564 + t565 - t319 / 0.12e2 + t324 / 0.12e2 + t327 / 0.12e2 - t330 / 0.12e2 - t334 / 0.12e2 + t337 / 0.12e2 + t572 - t573 + t574 - t575 - t576 - t577;
t579 = t127 / 0.24e2;
t580 = t129 / 0.24e2;
t581 = t131 / 0.24e2;
t582 = t134 / 0.24e2;
t583 = t136 / 0.24e2;
t584 = t139 / 0.24e2;
t597 = t579 - t580 + t581 - t582 + t583 + t584 - t396 / 0.12e2 - t398 / 0.12e2 - t400 / 0.12e2 - t403 / 0.12e2 + t405 / 0.12e2 - t407 / 0.12e2 + t412 / 0.12e2 + t415 / 0.12e2 + t417 / 0.12e2 - t420 / 0.12e2 + t422 / 0.12e2 + t424 / 0.12e2;
t598 = t578 + t597;
t605 = bb(2) * ee(1);
t606 = t605 * t332;
t607 = t248 * t606;
t609 = t350 * t317;
t610 = t246 * t609;
t612 = t345 * t317;
t613 = t244 * t612;
t615 = -t560 + t561 + t562 - t563 - t564 + t565 - t73 / 0.12e2 - t77 / 0.12e2 - t81 / 0.12e2 + t85 / 0.12e2 + t89 / 0.12e2 + t93 / 0.12e2 + t607 / 0.12e2 + t610 / 0.12e2 - t613 / 0.12e2;
t617 = t248 * t355 * t322;
t619 = t244 * t606;
t621 = bb(1) * ee(2);
t622 = t621 * t332;
t623 = t246 * t622;
t625 = bb(1) * ee(3);
t626 = t625 * t322;
t627 = t244 * t626;
t629 = ee(1) * dd(2);
t631 = t364 * t40 * t629;
t633 = t248 * t612;
t635 = t248 * t609;
t637 = t244 * t622;
t639 = t246 * t626;
t641 = -t617 / 0.12e2 + t619 / 0.12e2 - t623 / 0.12e2 + t627 / 0.12e2 - t631 / 0.12e2 - t633 / 0.12e2 + t635 / 0.12e2 - t637 / 0.12e2 + t639 / 0.12e2 + t572 - t573 + t574 - t575 - t576 - t577;
t652 = t579 - t580 + t581 - t582 + t583 + t584 + t143 / 0.12e2 + t147 / 0.12e2 - t150 / 0.12e2 + t153 / 0.12e2 + t155 / 0.12e2 + t158 / 0.12e2 - t161 / 0.12e2 - t163 / 0.12e2 + t165 / 0.12e2;
t654 = ee(3) * dd(2) * t51;
t655 = t133 * t654;
t657 = ee(2) * dd(1) * t51;
t658 = t118 * t657;
t660 = ee(1) * dd(3) * t51;
t661 = t95 * t660;
t663 = ee(3) * dd(1) * t51;
t664 = t115 * t663;
t666 = ee(2) * dd(3) * t51;
t667 = t133 * t666;
t668 = t629 * t51;
t669 = t118 * t668;
t670 = t115 * t660;
t671 = t138 * t668;
t672 = t124 * t666;
t673 = t138 * t657;
t674 = t95 * t663;
t675 = t124 * t654;
t676 = -t168 - t170 - t172 - t655 - t658 - t661 + t664 + t667 + t669 - t670 + t671 + t672 - t673 + t674 - t675;
t679 = t6 / 0.24e2;
t680 = t9 / 0.24e2;
t681 = t18 / 0.24e2;
t682 = t19 / 0.24e2;
t683 = t24 / 0.24e2;
t684 = t27 / 0.24e2;
t685 = t28 / 0.24e2;
t686 = t31 / 0.24e2;
t687 = t42 / 0.24e2;
t688 = t47 / 0.24e2;
t689 = t48 / 0.24e2;
t690 = t49 / 0.24e2;
t691 = -t679 - t680 + t681 + t682 + t683 - t684 - t685 + t686 - t687 + t688 - t689 + t690 + t211 - t215 + t219;
t692 = t621 * dd(3);
t694 = t625 * dd(2);
t698 = -t223 + t228 - t230 + t234 - t238 - t245 + t247 - t249 + t254 - t258 + t261 + t244 * t692 - t244 * t694 + t246 * t692 - t246 * t694;
t700 = t345 * dd(1);
t703 = t605 * dd(3);
t707 = t220 * t426 / 0.6e1;
t709 = t21 * t443 / 0.6e1;
t711 = t8 * t167 / 0.12e2;
t713 = t8 * t149 / 0.12e2;
t715 = t21 * t146 / 0.12e2;
t717 = t21 * t157 / 0.12e2;
t718 = t355 * dd(2);
t721 = t350 * dd(1);
t729 = t207 * t447 / 0.6e1;
t731 = t8 * t429 / 0.6e1;
t733 = t21 * t432 / 0.6e1;
t734 = -t181 * t700 / 0.6e1 + t181 * t703 / 0.6e1 - t707 + t709 + t711 + t713 - t715 - t717 - t177 * t718 / 0.6e1 + t177 * t721 / 0.6e1 - t202 * t718 / 0.6e1 + t202 * t721 / 0.6e1 - t729 + t731 - t733;
t736 = t220 * t438 / 0.6e1;
t738 = t3 * t429 / 0.6e1;
t740 = t225 * t457 / 0.6e1;
t741 = t265 / 0.12e2;
t742 = t266 / 0.12e2;
t743 = t267 / 0.12e2;
t744 = t270 / 0.12e2;
t745 = t340 * dd(2);
t748 = t365 * dd(1);
t756 = t3 * t149 / 0.12e2;
t758 = t3 * t160 / 0.12e2;
t759 = cc(1) * ee(2);
t760 = t759 * dd(3);
t763 = cc(1) * ee(3);
t764 = t763 * dd(2);
t767 = t736 - t738 + t740 + t741 - t742 + t743 - t744 - t8 * t745 / 0.12e2 + t8 * t748 / 0.12e2 - t21 * t745 / 0.12e2 + t21 * t748 / 0.12e2 + t756 + t758 - t3 * t760 / 0.12e2 + t3 * t764 / 0.12e2;
t771 = t8 * t157 / 0.12e2;
t777 = t8 * t142 / 0.12e2;
t787 = t225 * t435 / 0.6e1;
t789 = t3 * t432 / 0.6e1;
t791 = t8 * t443 / 0.6e1;
t793 = t207 * t452 / 0.6e1;
t795 = t21 * t160 / 0.12e2;
t797 = t21 * t167 / 0.12e2;
t798 = t360 * dd(1);
t801 = -t771 - t8 * t760 / 0.12e2 + t8 * t764 / 0.12e2 - t777 - t181 * t692 / 0.6e1 + t181 * t694 / 0.6e1 - t177 * t692 / 0.6e1 + t177 * t694 / 0.6e1 - t787 + t789 - t791 + t793 + t795 + t797 - t21 * t798 / 0.12e2;
t802 = cc(2) * ee(1);
t803 = t802 * dd(3);
t807 = t3 * t142 / 0.12e2;
t813 = t3 * t146 / 0.12e2;
t831 = t244 * t142 / 0.6e1;
t833 = t246 * t149 / 0.6e1;
t834 = t21 * t803 / 0.12e2 - t807 - t3 * t798 / 0.12e2 + t3 * t803 / 0.12e2 - t813 - t202 * t700 / 0.6e1 + t202 * t703 / 0.6e1 + t248 * t700 / 0.12e2 - t248 * t703 / 0.12e2 + t246 * t718 / 0.12e2 - t246 * t721 / 0.12e2 + t248 * t718 / 0.12e2 - t248 * t721 / 0.12e2 + t831 - t833;
t838 = t251 * t252 * ee(3) / 0.6e1;
t842 = bb(3) * ee(3);
t848 = t256 * t252 * ee(2) / 0.6e1;
t849 = bb(2) * ee(2);
t857 = t248 * t146 / 0.6e1;
t864 = t256 * t259 * ee(1) / 0.6e1;
t868 = bb(1) * ee(1);
t872 = t250 * cc(2);
t882 = -t838 + t251 * t350 * dd(3) / 0.6e1 - t251 * t842 * dd(2) / 0.6e1 + t848 + t256 * t849 * dd(3) / 0.6e1 - t256 * t345 * dd(2) / 0.6e1 + t857 + t244 * t700 / 0.12e2 - t244 * t703 / 0.12e2 - t864 + t256 * t625 * dd(1) / 0.6e1 - t256 * t868 * dd(3) / 0.6e1 + t872 * t842 * dd(1) / 0.6e1 - t872 * t355 * dd(3) / 0.6e1 + t872 * t605 * dd(2) / 0.6e1;
t892 = t13 / 0.24e2;
t893 = t15 / 0.24e2;
t894 = t35 / 0.24e2;
t895 = t39 / 0.24e2;
t896 = t44 / 0.24e2;
t897 = t46 / 0.24e2;
t899 = t37 * t759 / 0.12e2;
t901 = t11 * t763 / 0.12e2;
t903 = t33 * t360 / 0.12e2;
t905 = t37 * t802 / 0.12e2;
t907 = t11 * t340 / 0.12e2;
t909 = t33 * t365 / 0.12e2;
t910 = -t872 * t849 * dd(1) / 0.6e1 + t251 * t868 * dd(2) / 0.6e1 - t251 * t621 * dd(1) / 0.6e1 - t892 + t893 + t894 - t895 + t896 - t897 - t899 + t901 - t903 + t905 - t907 + t909;
t914 = t679 + t680 - t681 - t682 - t683 + t684 + t685 - t686 + t687 - t688 + t689 - t690 - t211 + t215 - t219 + t224 - t229 + t231;
t921 = -t235 + t239 + t245 / 0.12e2 - t247 / 0.12e2 + t249 / 0.12e2 - t254 / 0.12e2 + t258 / 0.12e2 - t261 / 0.12e2 + t707 - t709 - t711 - t713 + t715 + t717 + t729 - t731 + t733 - t736;
t923 = t738 - t740 - t741 + t742 - t743 + t744 - t756 - t758 + t771 + t777 + t787 - t789 + t791 - t793 - t795 - t797 + t807 + t813;
t924 = -t831 + t833 + t838 - t848 - t857 + t864 + t892 - t893 - t894 + t895 - t896 + t897 + t899 - t901 + t903 - t905 + t907 - t909;
t927 = -t607 + t655 + t658 + t661 - t610 - t664 + t613 - t667 + t617 - t669 - t619 + t623;
t928 = -t627 + t670 - t671 - t672 + t673 - t674 + t631 + t633 + t675 - t635 + t637 - t639;
t948 = t53 / 0.12e2 - t56 / 0.12e2 - t59 / 0.12e2 + t62 / 0.12e2 + t65 / 0.12e2 - t68 / 0.12e2 + t319 / 0.6e1 - t324 / 0.6e1 - t327 / 0.6e1 + t330 / 0.6e1 + t334 / 0.6e1 - t337 / 0.6e1 - t98 / 0.12e2 + t103 / 0.12e2 - t108 / 0.12e2 + t113 / 0.12e2 + t116 / 0.12e2 + t121 / 0.12e2;
t967 = -t127 / 0.12e2 + t129 / 0.12e2 - t131 / 0.12e2 + t134 / 0.12e2 - t136 / 0.12e2 - t139 / 0.12e2 + t396 / 0.6e1 + t398 / 0.6e1 + t400 / 0.6e1 + t403 / 0.6e1 - t405 / 0.6e1 + t407 / 0.6e1 - t412 / 0.6e1 - t415 / 0.6e1 - t417 / 0.6e1 + t420 / 0.6e1 - t422 / 0.6e1 - t424 / 0.6e1;
t968 = t948 + t967;
t971 = -t178 / 0.6e1 + t182 / 0.6e1 - t679 - t680 + t681 + t682 + t683 - t684 - t685 + t686 - t687 + t688 - t689;
t985 = t690 + t197 / 0.6e1 + t199 / 0.6e1 - t203 / 0.6e1 - t205 / 0.6e1 + t210 / 0.6e1 - t214 / 0.6e1 + t218 / 0.6e1 - t223 / 0.6e1 + t228 / 0.6e1 - t230 / 0.6e1 + t234 / 0.6e1 - t238 / 0.6e1 + t240 / 0.6e1;
t995 = -t274 / 0.6e1 - t276 / 0.6e1 - t892 + t893 + t894 - t895 + t896 - t897 - t286 / 0.6e1 + t290 / 0.6e1 - t294 / 0.6e1 + t296 / 0.6e1 - t298 / 0.6e1 + t300 / 0.6e1;
%
coef(1) = t50 / 0.16e2;
coef(2) = t123 + t174;
coef(3) = t195 + t242 + t273 / 0.12e2 + t302;
coef(4) = t321 + t372 + t387 + t409 + t446 + t481 + t509 + t556;
coef(5) = t598;
coef(6) = -t50 / 0.8e1;
coef(7) = t615 + t641 + t652 + t676 / 0.12e2;
coef(8) = t691 + t698 / 0.12e2 + t734 + t767 + t801 + t834 + t882 + t910;
coef(9) = -t598;
coef(10) = t914 + t921 + t923 + t924;
coef(11) = t927 / 0.6e1 + t928 / 0.6e1;
coef(12) = t50 / 0.4e1;
coef(13) = t968;
coef(14) = -t50 / 0.8e1;
coef(15) = t971 + t985 - t273 / 0.6e1 + t995;
coef(16) = -t968;
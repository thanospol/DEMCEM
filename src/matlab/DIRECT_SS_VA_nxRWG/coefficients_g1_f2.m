function [coef] = coefficients_g1_f2(r1,r2,r3,r4,r5, AreaT)
%% coefficients_g1_f2

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

coef  = zeros(11);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = 0.1e1 / AreaT;
t3 = t1 * cc(3) * bb(3);
t4 = cc(2) * ee(3);
t5 = t4 * dd(1);
t8 = cc(2) * ee(1);
t9 = t8 * dd(3);
t12 = t1 * bb(2);
t13 = t12 * bb(3);
t17 = cc(2) * ee(2);
t21 = t1 * bb(1);
t22 = t21 * bb(3);
t23 = cc(1) * ee(1);
t27 = cc(1) * ee(2);
t32 = t1 * cc(2) * bb(2);
t33 = cc(3) * ee(1);
t34 = t33 * dd(2);
t37 = cc(3) * ee(2);
t38 = t37 * dd(1);
t42 = t1 * cc(1) * bb(1);
t50 = cc(3) * ee(3);
t54 = t3 * t5 / 0.6e1 - t3 * t9 / 0.6e1 + t13 * t8 * dd(2) / 0.12e2 - t13 * t17 * dd(1) / 0.12e2 + t22 * t23 * dd(2) / 0.12e2 - t22 * t27 * dd(1) / 0.12e2 + t32 * t34 / 0.6e1 - t32 * t38 / 0.6e1 + t42 * t34 / 0.6e1 - t42 * t38 / 0.6e1 + t22 * t37 * dd(3) / 0.12e2 - t22 * t50 * dd(2) / 0.12e2;
t55 = t12 * bb(1);
t62 = t27 * dd(3);
t65 = cc(1) * ee(3);
t66 = t65 * dd(2);
t89 = t55 * t17 * dd(3) / 0.12e2 - t55 * t4 * dd(2) / 0.12e2 + t3 * t62 / 0.6e1 - t3 * t66 / 0.6e1 + t32 * t62 / 0.6e1 - t32 * t66 / 0.6e1 + t55 * t65 * dd(1) / 0.12e2 - t55 * t23 * dd(3) / 0.12e2 + t13 * t50 * dd(1) / 0.12e2 - t13 * t33 * dd(3) / 0.12e2 + t42 * t5 / 0.6e1 - t42 * t9 / 0.6e1;
t91 = cc(2) ^ 2;
t92 = t1 * t91;
t94 = bb(1) * ee(3) * dd(2);
t98 = bb(1) * ee(2) * dd(3);
t101 = cc(3) ^ 2;
t102 = t1 * t101;
t107 = bb(3) ^ 2;
t108 = t1 * t107;
t113 = bb(2) ^ 2;
t114 = t1 * t113;
t119 = bb(1) ^ 2;
t120 = t1 * t119;
t129 = t92 * t94 / 0.6e1 - t92 * t98 / 0.6e1 + t102 * t94 / 0.6e1 - t102 * t98 / 0.6e1 - t108 * t62 / 0.12e2 + t108 * t66 / 0.12e2 - t114 * t62 / 0.12e2 + t114 * t66 / 0.12e2 - t120 * t5 / 0.12e2 + t120 * t9 / 0.12e2 - t108 * t5 / 0.12e2 + t108 * t9 / 0.12e2;
t130 = cc(1) ^ 2;
t131 = t1 * t130;
t132 = bb(2) * ee(3);
t133 = t132 * dd(1);
t137 = bb(2) * ee(1) * dd(3);
t152 = bb(3) * ee(1);
t153 = t152 * dd(2);
t156 = bb(3) * ee(2);
t157 = t156 * dd(1);
t164 = -t131 * t133 / 0.6e1 + t131 * t137 / 0.6e1 - t102 * t133 / 0.6e1 + t102 * t137 / 0.6e1 - t114 * t34 / 0.12e2 + t114 * t38 / 0.12e2 - t120 * t34 / 0.12e2 + t120 * t38 / 0.12e2 - t92 * t153 / 0.6e1 + t92 * t157 / 0.6e1 - t131 * t153 / 0.6e1 + t131 * t157 / 0.6e1;
t167 = t113 * bb(2);
t168 = t1 * t167;
t169 = cc(1) * dd(3);
t170 = t168 * t169;
t171 = t170 / 0.24e2;
t172 = t101 * cc(3);
t173 = t1 * t172;
t174 = bb(1) * dd(2);
t175 = t173 * t174;
t177 = t91 * cc(2);
t178 = t1 * t177;
t179 = bb(1) * dd(3);
t180 = t178 * t179;
t182 = t119 * bb(1);
t183 = t1 * t182;
t184 = cc(2) * dd(3);
t185 = t183 * t184;
t186 = t185 / 0.24e2;
t187 = t107 * bb(3);
t188 = t1 * t187;
t189 = cc(2) * dd(1);
t190 = t188 * t189;
t191 = t190 / 0.24e2;
t192 = t130 * cc(1);
t193 = t1 * t192;
t194 = bb(2) * dd(3);
t195 = t193 * t194;
t197 = bb(2) * dd(1);
t198 = t173 * t197;
t200 = cc(3) * dd(1);
t201 = t168 * t200;
t202 = t201 / 0.24e2;
t203 = cc(3) * dd(2);
t204 = t183 * t203;
t205 = t204 / 0.24e2;
t206 = bb(3) * dd(1);
t207 = t178 * t206;
t209 = bb(3) * dd(2);
t210 = t193 * t209;
t212 = cc(1) * dd(2);
t213 = t188 * t212;
t214 = t213 / 0.24e2;
t215 = cc(3) * cc(2);
t217 = t22 * t215 * dd(3);
t219 = -t171 - t175 / 0.6e1 + t180 / 0.6e1 + t186 - t191 - t195 / 0.6e1 + t198 / 0.6e1 + t202 - t205 - t207 / 0.6e1 + t210 / 0.6e1 + t214 - t217 / 0.6e1;
t221 = t55 * t215 * dd(2);
t223 = bb(2) * cc(1);
t224 = t223 * dd(3);
t225 = t3 * t224;
t227 = bb(3) * cc(1);
t228 = t227 * dd(2);
t229 = t32 * t228;
t231 = cc(1) * cc(3);
t233 = t55 * t231 * dd(1);
t235 = t108 * t224;
t236 = t235 / 0.24e2;
t237 = bb(1) * cc(3);
t238 = t237 * dd(2);
t239 = t108 * t238;
t240 = t239 / 0.24e2;
t241 = bb(1) * cc(2);
t242 = t241 * dd(3);
t243 = t114 * t242;
t244 = t243 / 0.24e2;
t245 = t114 * t228;
t246 = t245 / 0.24e2;
t247 = t102 * t228;
t249 = t102 * t242;
t251 = t92 * t238;
t253 = t92 * t224;
t255 = bb(3) * cc(2);
t256 = t255 * dd(1);
t257 = t120 * t256;
t258 = t257 / 0.24e2;
t259 = t120 * t224;
t260 = t259 / 0.24e2;
t261 = t221 / 0.6e1 + t225 / 0.6e1 - t229 / 0.6e1 - t233 / 0.6e1 - t236 - t240 + t244 + t246 + t247 / 0.6e1 + t249 / 0.6e1 - t251 / 0.6e1 - t253 / 0.6e1 - t258 - t260;
t263 = bb(2) * cc(3);
t264 = t263 * dd(1);
t265 = t108 * t264;
t266 = t265 / 0.24e2;
t267 = t108 * t242;
t268 = t267 / 0.24e2;
t269 = t131 * t242;
t271 = t131 * t264;
t273 = t102 * t224;
t275 = t102 * t256;
t277 = t114 * t238;
t278 = t277 / 0.24e2;
t279 = t114 * t256;
t280 = t279 / 0.24e2;
t281 = t120 * t228;
t282 = t281 / 0.24e2;
t283 = t120 * t264;
t284 = t283 / 0.24e2;
t285 = t92 * t264;
t287 = t92 * t228;
t289 = t131 * t256;
t291 = t266 + t268 + t269 / 0.6e1 + t271 / 0.6e1 - t273 / 0.6e1 - t275 / 0.6e1 - t278 - t280 + t282 + t284 + t285 / 0.6e1 + t287 / 0.6e1 - t289 / 0.6e1;
t293 = t231 * dd(2);
t294 = t108 * t293;
t295 = t1 * bb(3);
t296 = bb(1) * t101;
t297 = t296 * dd(2);
t298 = t295 * t297;
t299 = bb(1) * t91;
t300 = t299 * dd(3);
t301 = t12 * t300;
t302 = cc(1) * cc(2);
t303 = t302 * dd(3);
t304 = t114 * t303;
t305 = t120 * t303;
t306 = bb(2) * t130;
t307 = t306 * dd(3);
t308 = t21 * t307;
t309 = bb(2) * t101;
t310 = t309 * dd(1);
t311 = t295 * t310;
t312 = t215 * dd(1);
t313 = t108 * t312;
t314 = t114 * t312;
t315 = bb(3) * t91;
t316 = t315 * dd(1);
t317 = t12 * t316;
t318 = bb(3) * t130;
t319 = t318 * dd(2);
t320 = t21 * t319;
t321 = t120 * t293;
t322 = t42 * t256;
t323 = -t131 * t238 - t294 + t298 - t301 + t304 - t305 + t308 - t311 + t313 - t314 + t317 - t320 + t321 + t322;
t326 = -t239 + t267 + t245 - t170 + t265 + t213 - t277 - t235 + t201 - t204 + t283 - t190 - t259 - t279 + t281 - t257 + t243 + t185;
t327 = sqrt(0.3e1);
t328 = t184 * t327;
t329 = t183 * t328;
t331 = t169 * t327;
t332 = t168 * t331;
t334 = t212 * t327;
t335 = t188 * t334;
t337 = t200 * t327;
t338 = t168 * t337;
t340 = t203 * t327;
t341 = t183 * t340;
t343 = t189 * t327;
t344 = t188 * t343;
t346 = t327 * t1;
t349 = t346 * t182 * cc(3) * ee(2);
t353 = t346 * t182 * cc(2) * ee(3);
t357 = t346 * t187 * cc(2) * ee(1);
t361 = t346 * t167 * cc(1) * ee(3);
t365 = t346 * t187 * cc(1) * ee(2);
t369 = t346 * t167 * cc(3) * ee(1);
t372 = t120 * bb(2) * t331;
t375 = t108 * bb(2) * t337;
t377 = t114 * cc(1);
t379 = t377 * t209 * t327;
t381 = t108 * cc(2);
t383 = t381 * t179 * t327;
t385 = t108 * cc(1);
t387 = t385 * t194 * t327;
t389 = t120 * cc(3);
t391 = t389 * t197 * t327;
t393 = t329 / 0.48e2 - t332 / 0.48e2 + t335 / 0.48e2 + t338 / 0.48e2 - t341 / 0.48e2 - t344 / 0.48e2 + t349 / 0.24e2 - t353 / 0.24e2 + t357 / 0.24e2 + t361 / 0.24e2 - t365 / 0.24e2 - t369 / 0.24e2 - t372 / 0.48e2 + t375 / 0.48e2 + t379 / 0.48e2 + t383 / 0.48e2 - t387 / 0.48e2 + t391 / 0.48e2;
t394 = t114 * cc(3);
t396 = t394 * t174 * t327;
t399 = t114 * bb(3) * t343;
t402 = t108 * bb(1) * t340;
t405 = t120 * bb(3) * t334;
t408 = t114 * bb(1) * t328;
t410 = t120 * cc(2);
t412 = t410 * t206 * t327;
t414 = t346 * t107;
t415 = t237 * ee(2);
t416 = t414 * t415;
t418 = t346 * t113;
t419 = t241 * ee(3);
t420 = t418 * t419;
t422 = t223 * ee(3);
t423 = t414 * t422;
t425 = t346 * t119;
t426 = t425 * t422;
t428 = t418 * t415;
t430 = t227 * ee(2);
t431 = t418 * t430;
t433 = t425 * t430;
t435 = t255 * ee(1);
t436 = t418 * t435;
t438 = t414 * t419;
t440 = t263 * ee(1);
t441 = t425 * t440;
t443 = t425 * t435;
t445 = t414 * t440;
t447 = -t396 / 0.48e2 - t399 / 0.48e2 - t402 / 0.48e2 + t405 / 0.48e2 + t408 / 0.48e2 - t412 / 0.48e2 + t416 / 0.24e2 - t420 / 0.24e2 + t423 / 0.24e2 + t426 / 0.24e2 + t428 / 0.24e2 - t431 / 0.24e2 - t433 / 0.24e2 + t436 / 0.24e2 - t438 / 0.24e2 - t441 / 0.24e2 + t443 / 0.24e2 - t445 / 0.24e2;
t461 = t217 / 0.12e2;
t462 = t170 / 0.48e2 + t175 / 0.12e2 - t180 / 0.12e2 - t185 / 0.48e2 + t190 / 0.48e2 + t195 / 0.12e2 - t198 / 0.12e2 - t201 / 0.48e2 + t204 / 0.48e2 + t207 / 0.12e2 - t210 / 0.12e2 - t213 / 0.48e2 + t461;
t463 = t221 / 0.12e2;
t464 = t225 / 0.12e2;
t465 = t229 / 0.12e2;
t466 = t233 / 0.12e2;
t477 = -t463 - t464 + t465 + t466 + t235 / 0.48e2 + t239 / 0.48e2 - t243 / 0.48e2 - t245 / 0.48e2 - t247 / 0.12e2 - t249 / 0.12e2 + t251 / 0.12e2 + t253 / 0.12e2 + t257 / 0.48e2 + t259 / 0.48e2;
t492 = -t265 / 0.48e2 - t267 / 0.48e2 - t269 / 0.12e2 - t271 / 0.12e2 + t273 / 0.12e2 + t275 / 0.12e2 + t277 / 0.48e2 + t279 / 0.48e2 - t281 / 0.48e2 - t283 / 0.48e2 - t285 / 0.12e2 - t287 / 0.12e2 + t289 / 0.12e2;
t495 = t346 * cc(1);
t497 = t495 * t241 * t206;
t499 = t346 * cc(2);
t501 = t499 * t223 * t209;
t503 = t346 * cc(3);
t505 = t503 * t227 * t194;
t507 = t346 * bb(2);
t509 = t507 * t241 * t203;
t511 = t346 * bb(3);
t513 = t511 * t237 * t184;
t515 = t346 * bb(1);
t517 = t515 * t223 * t200;
t538 = t497 / 0.36e2 - t501 / 0.36e2 + t505 / 0.36e2 + t509 / 0.36e2 - t513 / 0.36e2 - t517 / 0.36e2 + t499 * t223 * t156 / 0.18e2 - t507 * t241 * t37 / 0.18e2 + t511 * t237 * t4 / 0.18e2 - t495 * t241 * t152 / 0.18e2 + t515 * t223 * t33 / 0.18e2 - t503 * t227 * t132 / 0.18e2 + t329 / 0.144e3;
t550 = t172 * bb(2);
t554 = t192 * bb(3);
t558 = t177 * bb(3);
t562 = -t332 / 0.144e3 + t335 / 0.144e3 + t338 / 0.144e3 - t341 / 0.144e3 - t344 / 0.144e3 + t349 / 0.72e2 - t353 / 0.72e2 + t357 / 0.72e2 + t361 / 0.72e2 - t365 / 0.72e2 - t369 / 0.72e2 + t346 * t550 * dd(1) / 0.36e2 + t346 * t554 * dd(2) / 0.36e2 - t346 * t558 * dd(1) / 0.36e2;
t564 = t192 * bb(2);
t568 = t172 * bb(1);
t572 = t177 * bb(1);
t598 = -t346 * t564 * dd(3) / 0.36e2 - t346 * t568 * dd(2) / 0.36e2 + t346 * t572 * dd(3) / 0.36e2 + t346 * t568 * ee(2) / 0.18e2 - t346 * t550 * ee(1) / 0.18e2 - t346 * t572 * ee(3) / 0.18e2 + t346 * t564 * ee(3) / 0.18e2 - t346 * t554 * ee(2) / 0.18e2 + t346 * t558 * ee(1) / 0.18e2 - t372 / 0.144e3 + t375 / 0.144e3 + t379 / 0.144e3 + t383 / 0.144e3;
t613 = -t387 / 0.144e3 + t391 / 0.144e3 - t396 / 0.144e3 - t399 / 0.144e3 - t402 / 0.144e3 + t405 / 0.144e3 + t408 / 0.144e3 - t412 / 0.144e3 + t416 / 0.72e2 - t420 / 0.72e2 + t423 / 0.72e2 + t426 / 0.72e2 + t428 / 0.72e2 - t431 / 0.72e2;
t622 = t511 * t310;
t624 = t425 * t293;
t626 = t414 * t312;
t628 = t414 * t293;
t630 = t515 * t307;
t632 = t507 * t316;
t634 = t418 * t312;
t636 = -t433 / 0.72e2 + t436 / 0.72e2 - t438 / 0.72e2 - t441 / 0.72e2 + t443 / 0.72e2 - t445 / 0.72e2 - t622 / 0.36e2 + t624 / 0.36e2 + t626 / 0.36e2 - t628 / 0.36e2 + t630 / 0.36e2 + t632 / 0.36e2 - t634 / 0.36e2;
t637 = t511 * t297;
t639 = t507 * t300;
t641 = t425 * t303;
t643 = t418 * t303;
t645 = t515 * t319;
t647 = t215 * ee(1);
t650 = t296 * ee(2);
t653 = t231 * ee(2);
t658 = t315 * ee(1);
t661 = t302 * ee(3);
t664 = t306 * ee(3);
t671 = t637 / 0.36e2 - t639 / 0.36e2 - t641 / 0.36e2 + t643 / 0.36e2 - t645 / 0.36e2 - t414 * t647 / 0.18e2 - t511 * t650 / 0.18e2 - t425 * t653 / 0.18e2 + t414 * t653 / 0.18e2 - t507 * t658 / 0.18e2 - t418 * t661 / 0.18e2 - t515 * t664 / 0.18e2 + t425 * t661 / 0.18e2 + t418 * t647 / 0.18e2;
t673 = t309 * ee(1);
t676 = t299 * ee(3);
t679 = t318 * ee(2);
t682 = t346 * t91;
t685 = t346 * t130;
t690 = t346 * t101;
t705 = t511 * t673 / 0.18e2 + t507 * t676 / 0.18e2 + t515 * t679 / 0.18e2 + t682 * t264 / 0.36e2 - t685 * t238 / 0.36e2 + t682 * t228 / 0.36e2 + t690 * t242 / 0.36e2 - t690 * t224 / 0.36e2 + t685 * t264 / 0.36e2 + t685 * t242 / 0.36e2 - t690 * t256 / 0.36e2 - t682 * t224 / 0.36e2 - t682 * t238 / 0.36e2;
t734 = -t685 * t256 / 0.36e2 + t690 * t228 / 0.36e2 + t685 * t415 / 0.18e2 - t690 * t419 / 0.18e2 - t685 * t440 / 0.18e2 - t685 * t419 / 0.18e2 + t685 * t435 / 0.18e2 - t682 * t430 / 0.18e2 + t682 * t415 / 0.18e2 + t682 * t422 / 0.18e2 + t690 * t422 / 0.18e2 - t682 * t440 / 0.18e2 - t690 * t430 / 0.18e2 + t690 * t435 / 0.18e2;
t738 = dd(1) * t327;
t741 = dd(3) * t327;
t745 = ee(2) * dd(1) * t327;
t747 = dd(2) * t327;
t751 = ee(3) * dd(1) * t327;
t754 = ee(3) * dd(2) * t327;
t757 = ee(2) * dd(3) * t327;
t762 = ee(1) * dd(2) * t327;
t765 = ee(1) * dd(3) * t327;
t770 = -t55 * t65 * t738 - t22 * t37 * t741 - t389 * t745 + t55 * t4 * t747 + t381 * t751 - t385 * t754 + t377 * t757 + t22 * t27 * t738 + t394 * t762 - t410 * t765 + t13 * t33 * t741 - t394 * t745;
t790 = t13 * t17 * t738 + t389 * t762 + t410 * t751 - t55 * t17 * t741 - t381 * t765 + t55 * t23 * t741 - t22 * t23 * t747 + t385 * t757 + t22 * t50 * t747 - t377 * t754 - t13 * t50 * t738 - t13 * t8 * t747;
t810 = -t497 / 0.6e1 + t501 / 0.6e1 - t505 / 0.6e1 - t509 / 0.6e1 + t513 / 0.6e1 + t517 / 0.6e1 - t329 / 0.12e2 + t332 / 0.12e2 - t335 / 0.12e2 - t338 / 0.12e2 + t341 / 0.12e2 + t344 / 0.12e2 + t372 / 0.12e2 - t375 / 0.12e2 - t379 / 0.12e2 - t383 / 0.12e2 + t387 / 0.12e2 - t391 / 0.12e2;
t829 = t396 / 0.12e2 + t399 / 0.12e2 + t402 / 0.12e2 - t405 / 0.12e2 - t408 / 0.12e2 + t412 / 0.12e2 + t622 / 0.6e1 - t624 / 0.6e1 - t626 / 0.6e1 + t628 / 0.6e1 - t630 / 0.6e1 - t632 / 0.6e1 + t634 / 0.6e1 - t637 / 0.6e1 + t639 / 0.6e1 + t641 / 0.6e1 - t643 / 0.6e1 + t645 / 0.6e1;
t849 = t497 / 0.12e2 - t501 / 0.12e2 + t505 / 0.12e2 + t509 / 0.12e2 - t513 / 0.12e2 - t517 / 0.12e2 + t329 / 0.24e2 - t332 / 0.24e2 + t335 / 0.24e2 + t338 / 0.24e2 - t341 / 0.24e2 - t344 / 0.24e2 - t372 / 0.24e2 + t375 / 0.24e2 + t379 / 0.24e2 + t383 / 0.24e2 - t387 / 0.24e2 + t391 / 0.24e2;
t868 = -t396 / 0.24e2 - t399 / 0.24e2 - t402 / 0.24e2 + t405 / 0.24e2 + t408 / 0.24e2 - t412 / 0.24e2 - t622 / 0.12e2 + t624 / 0.12e2 + t626 / 0.12e2 - t628 / 0.12e2 + t630 / 0.12e2 + t632 / 0.12e2 - t634 / 0.12e2 + t637 / 0.12e2 - t639 / 0.12e2 - t641 / 0.12e2 + t643 / 0.12e2 - t645 / 0.12e2;
t882 = t171 - t186 + t191 - t202 + t205 - t214 - t168 * t65 / 0.12e2 + t183 * t4 / 0.12e2 - t188 * t8 / 0.12e2 + t168 * t33 / 0.12e2 - t183 * t37 / 0.12e2 + t188 * t27 / 0.12e2 + t461 - t463 - t464 + t465 + t466 + t236;
t890 = t240 - t244 - t246 + t258 + t260 - t266 - t268 + t278 + t280 - t282 - t284 + t294 / 0.12e2 - t298 / 0.12e2 + t301 / 0.12e2 - t304 / 0.12e2 + t305 / 0.12e2 - t308 / 0.12e2 + t311 / 0.12e2;
t923 = -t313 / 0.12e2 + t314 / 0.12e2 - t317 / 0.12e2 + t320 / 0.12e2 - t321 / 0.12e2 - t108 * t422 / 0.12e2 - t108 * t415 / 0.12e2 + t114 * t419 / 0.12e2 + t114 * t430 / 0.12e2 - t120 * t435 / 0.12e2 - t120 * t422 / 0.12e2 + t108 * t440 / 0.12e2 + t108 * t419 / 0.12e2 - t114 * t415 / 0.12e2 - t114 * t435 / 0.12e2 + t120 * t430 / 0.12e2 + t120 * t440 / 0.12e2 - t108 * t653 / 0.6e1;
t962 = t295 * t650 / 0.6e1 - t12 * t676 / 0.6e1 + t114 * t661 / 0.6e1 - t120 * t661 / 0.6e1 + t21 * t664 / 0.6e1 - t295 * t673 / 0.6e1 + t108 * t647 / 0.6e1 - t114 * t647 / 0.6e1 + t12 * t658 / 0.6e1 - t21 * t679 / 0.6e1 + t120 * t653 / 0.6e1 - t322 / 0.12e2 - t22 * t215 * ee(3) / 0.6e1 + t55 * t215 * ee(2) / 0.6e1 + t3 * t422 / 0.6e1 - t32 * t430 / 0.6e1 - t55 * t231 * ee(1) / 0.6e1 + t42 * t435 / 0.6e1;
%
coef(1) = t54 + t89 + t129 + t164;
coef(2) = t219 + t261 + t291 + t323 / 0.6e1;
coef(3) = t326 / 0.8e1;
coef(4) = -t326 / 0.16e2;
coef(5) = t393 + t447;
coef(6) = t462 + t477 + t492 - t323 / 0.12e2;
coef(7) = t538 + t562 + t598 + t613 + t636 + t671 + t705 + t734;
coef(8) = t770 / 0.12e2 + t790 / 0.12e2;
coef(9) = t810 + t829;
coef(10) = t849 + t868;
coef(11) = t882 + t890 + t923 + t962;
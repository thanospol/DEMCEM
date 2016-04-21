function [coef,coefm] = coefficients_g3_f2(r1,r2,r3,r4,ko,AreaT)
%% coefficients_g3_f2

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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% ko = wavenumber
% AreaT = Area of outer triangle

% OUTPUT DATA
% coef   
% coefm
%

c  = zeros(6);
cm = zeros(6);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
%
t1 = sqrt(0.3e1);
t2 = 0.1e1 / AreaT;
t3 = t1 * t2;
t4 = bb(3) ^ 2;
t5 = t3 * t4;
t6 = cc(2) * bb(1);
t7 = t6 * dd(3);
t8 = t5 * t7;
t10 = cc(1) * bb(2) * dd(3);
t11 = t5 * t10;
t12 = cc(3) * bb(2);
t13 = t12 * dd(1);
t14 = t5 * t13;
t15 = bb(2) ^ 2;
t16 = t15 * bb(2);
t19 = t3 * t16 * cc(3) * dd(1);
t20 = t3 * t15;
t21 = t20 * t7;
t22 = t4 * bb(3);
t25 = t3 * t22 * cc(1) * dd(2);
t27 = cc(2) * bb(3) * dd(1);
t28 = t20 * t27;
t29 = bb(1) ^ 2;
t30 = t3 * t29;
t31 = t30 * t13;
t33 = cc(1) * bb(3) * dd(2);
t34 = t30 * t33;
t37 = t3 * t22 * cc(2) * dd(1);
t38 = t20 * t33;
t39 = cc(3) * bb(1);
t40 = t39 * dd(2);
t41 = t20 * t40;
t42 = t5 * t40;
t43 = t30 * t10;
t44 = t30 * t27;
t47 = t3 * t16 * cc(1) * dd(3);
t48 = t29 * bb(1);
t51 = t3 * t48 * cc(2) * dd(3);
t54 = t3 * t48 * cc(3) * dd(2);
t55 = -t8 + t11 - t14 - t19 - t21 - t25 + t28 - t31 - t34 + t37 - t38 + t41 + t42 + t43 + t44 + t47 - t51 + t54;
t68 = cc(3) ^ 2;
t69 = t2 * t68;
t71 = bb(3) * dd(2);
t73 = t69 * bb(1) * t71 * t1;
t75 = -t8 / 0.72e2 + t11 / 0.72e2 - t14 / 0.72e2 - t21 / 0.72e2 + t28 / 0.72e2 - t31 / 0.72e2 - t34 / 0.72e2 - t38 / 0.72e2 + t41 / 0.72e2 + t42 / 0.72e2 + t43 / 0.72e2 + t44 / 0.72e2 - t73 / 0.18e2;
t76 = t2 * cc(2);
t77 = t76 * cc(1);
t80 = t77 * t15 * dd(3) * t1;
t81 = t2 * cc(3);
t82 = t81 * cc(1);
t85 = t82 * t4 * dd(2) * t1;
t86 = cc(2) ^ 2;
t87 = t2 * t86;
t89 = bb(2) * dd(3);
t91 = t87 * bb(1) * t89 * t1;
t93 = bb(3) * dd(1);
t95 = t69 * bb(2) * t93 * t1;
t96 = cc(1) ^ 2;
t97 = t2 * t96;
t99 = bb(1) * dd(3);
t101 = t97 * bb(2) * t99 * t1;
t104 = t77 * t29 * dd(3) * t1;
t105 = t81 * cc(2);
t108 = t105 * t4 * dd(1) * t1;
t110 = bb(2) * dd(1);
t112 = t87 * bb(3) * t110 * t1;
t115 = t105 * t15 * dd(1) * t1;
t118 = t82 * t29 * dd(2) * t1;
t120 = bb(1) * dd(2);
t122 = t97 * bb(3) * t120 * t1;
t123 = t3 * t96;
t125 = t3 * t86;
t127 = t3 * t68;
t129 = -t80 + t85 + t91 + t95 - t101 + t104 - t108 - t112 + t115 - t118 + t122 + t123 * t27 + t125 * t40 - t127 * t33;
t141 = t2 * bb(3) * cc(3);
t142 = dd(3) * t1;
t144 = t141 * t6 * t142;
t145 = bb(3) * bb(1);
t146 = dd(1) * t1;
t148 = t77 * t145 * t146;
t149 = bb(2) * bb(3);
t150 = dd(2) * t1;
t152 = t77 * t149 * t150;
t154 = t2 * bb(1) * cc(1);
t156 = t154 * t12 * t146;
t157 = t125 * t10 - t127 * t7 + t123 * t40 - t123 * t13 - t123 * t7 + t127 * t27 - t125 * t33 - t125 * t13 + t127 * t10 + t144 - t148 + t152 + t156;
t159 = t2 * bb(2) * cc(2);
t161 = t159 * t39 * t150;
t164 = t82 * t149 * t142;
t172 = t96 * cc(1);
t177 = t68 * cc(3);
t186 = t86 * cc(2);
t199 = -t161 / 0.18e2 - t164 / 0.18e2 - t19 / 0.72e2 - t25 / 0.72e2 + t37 / 0.72e2 + t47 / 0.72e2 - t51 / 0.72e2 + t54 / 0.72e2 + t3 * t172 * bb(2) * dd(3) / 0.18e2 - t3 * t177 * bb(2) * dd(1) / 0.18e2 + t3 * t177 * bb(1) * dd(2) / 0.18e2 - t3 * t186 * bb(1) * dd(3) / 0.18e2 + t3 * t186 * bb(3) * dd(1) / 0.18e2 - t3 * t172 * bb(3) * dd(2) / 0.18e2;
t220 = -t8 / 0.24e2 + t11 / 0.24e2 - t14 / 0.24e2 - t21 / 0.24e2 + t28 / 0.24e2 - t31 / 0.24e2 - t34 / 0.24e2 - t38 / 0.24e2 + t41 / 0.24e2 + t42 / 0.24e2 + t43 / 0.24e2 + t44 / 0.24e2 - t73 / 0.12e2 - t80 / 0.12e2 + t85 / 0.12e2 + t91 / 0.12e2 + t95 / 0.12e2 - t101 / 0.12e2;
t239 = t104 / 0.12e2 - t108 / 0.12e2 - t112 / 0.12e2 + t115 / 0.12e2 - t118 / 0.12e2 + t122 / 0.12e2 + t144 / 0.12e2 - t148 / 0.12e2 + t152 / 0.12e2 + t156 / 0.12e2 - t161 / 0.12e2 - t164 / 0.12e2 - t19 / 0.24e2 - t25 / 0.24e2 + t37 / 0.24e2 + t47 / 0.24e2 - t51 / 0.24e2 + t54 / 0.24e2;
t241 = t2 * t22;
t244 = t241 * cc(1) * dd(2) / 0.12e2;
t245 = t2 * t16;
t248 = t245 * cc(1) * dd(3) / 0.12e2;
t249 = t2 * t177;
t252 = t2 * t186;
t255 = t2 * t48;
t258 = t255 * cc(2) * dd(3) / 0.12e2;
t261 = t241 * cc(2) * dd(1) / 0.12e2;
t262 = t2 * t172;
t269 = t245 * cc(3) * dd(1) / 0.12e2;
t272 = t255 * cc(3) * dd(2) / 0.12e2;
t277 = t141 * t7;
t279 = t244 - t248 - t249 * t120 / 0.6e1 + t252 * t99 / 0.6e1 + t258 - t261 - t262 * t89 / 0.6e1 + t249 * t110 / 0.6e1 + t269 - t272 - t252 * t93 / 0.6e1 + t262 * t71 / 0.6e1 - t277 / 0.4e1;
t280 = t159 * t40;
t283 = t82 * t149 * dd(3);
t286 = t77 * t149 * dd(2);
t288 = t154 * t13;
t291 = t77 * t145 * dd(1);
t295 = t81 * cc(2) * t4 * dd(1);
t297 = t149 * dd(1);
t298 = t69 * t297;
t300 = t2 * t15;
t302 = t300 * t27 / 0.12e2;
t304 = t300 * t40 / 0.12e2;
t305 = t2 * t29;
t307 = t305 * t13 / 0.12e2;
t309 = t305 * t33 / 0.12e2;
t316 = t280 / 0.4e1 + t283 / 0.4e1 - t286 / 0.4e1 - t288 / 0.4e1 + t291 / 0.4e1 + t295 / 0.4e1 - t298 / 0.4e1 - t302 - t304 + t307 + t309 + t87 * t33 / 0.6e1 + t87 * t13 / 0.6e1 - t97 * t40 / 0.6e1;
t320 = t87 * t297;
t324 = t76 * cc(3) * t15 * dd(1);
t326 = t2 * cc(1);
t329 = t326 * cc(3) * t29 * dd(2);
t331 = t145 * dd(2);
t332 = t97 * t331;
t334 = t2 * t4;
t336 = t334 * t40 / 0.12e2;
t338 = t334 * t10 / 0.12e2;
t340 = t300 * t33 / 0.12e2;
t342 = t300 * t7 / 0.12e2;
t351 = -t97 * t27 / 0.6e1 + t320 / 0.4e1 - t324 / 0.4e1 + t329 / 0.4e1 - t332 / 0.4e1 - t336 - t338 + t340 + t342 + t69 * t7 / 0.6e1 + t69 * t33 / 0.6e1 - t87 * t10 / 0.6e1 - t87 * t40 / 0.6e1;
t352 = t69 * t331;
t356 = t81 * cc(1) * t4 * dd(2);
t360 = t76 * cc(1) * t15 * dd(3);
t363 = bb(1) * bb(2) * dd(3);
t364 = t87 * t363;
t367 = t305 * t10 / 0.12e2;
t369 = t305 * t27 / 0.12e2;
t371 = t334 * t7 / 0.12e2;
t373 = t334 * t13 / 0.12e2;
t382 = t97 * t363;
t386 = t326 * cc(2) * t29 * dd(3);
t388 = t352 / 0.4e1 - t356 / 0.4e1 + t360 / 0.4e1 - t364 / 0.4e1 - t367 - t369 + t371 + t373 + t97 * t13 / 0.6e1 + t97 * t7 / 0.6e1 - t69 * t27 / 0.6e1 - t69 * t10 / 0.6e1 + t382 / 0.4e1 - t386 / 0.4e1;
t399 = t244 - t248 + t258 - t261 + t269 - t272 - t277 / 0.6e1 + t280 / 0.6e1 + t283 / 0.6e1 - t286 / 0.6e1 - t288 / 0.6e1 + t291 / 0.6e1 + t295 / 0.6e1 - t298 / 0.6e1 - t302 - t304 + t307 + t309;
t410 = t320 / 0.6e1 - t324 / 0.6e1 + t329 / 0.6e1 - t332 / 0.6e1 - t336 - t338 + t340 + t342 + t352 / 0.6e1 - t356 / 0.6e1 + t360 / 0.6e1 - t364 / 0.6e1 - t367 - t369 + t371 + t373 + t382 / 0.6e1 - t386 / 0.6e1;
t412 = -t8 + t11 - t14 - t21 + t28 - t31 - t34 - t38 + t41 + t42 + t43 + t44 - t73 - t80 + t85 + t91 + t95 - t101;
t413 = t104 - t108 - t112 + t115 - t118 + t122 + t144 - t148 + t152 + t156 - t161 - t164 - t19 - t25 + t37 + t47 - t51 + t54;
%
c(1) = t55 / 0.24e2;
c(2) = t75 + t129 / 0.18e2 + t157 / 0.18e2 + t199;
c(3) = t220 + t239;
c(4) = t279 + t316 + t351 + t388;
c(5) = t399 + t410;
c(6) = t412 / 0.12e2 + t413 / 0.12e2;
%
t1 = (ko ^ 2);
t3 = (j ^ 2);
t5 = (0.1e1 / t1 / t3);
t7 = 3 * t5 * c(5);
t10 = (0.1e1 / ko / j);
t14 = 8 * t5 * c(2);
t16 = 2 * t10 * c(6);
t18 = 4 * t10 * c(2);
t20 = 2 * t10 * c(1);
t22 = 3 * t5 * c(4);
t29 = 8 / t1 / ko / t3 / j * c(2);
t31 = 2 * t10 * c(3);
%
coef(1) = c(6);
coef(2) = t7;
coef(3) = 3 * t10 * c(4);
coef(4) = c(4);
coef(5) = c(1);
coef(6) = t14;
coef(7) = t16;
coef(8) = -t7;
coef(9) = -t14;
coef(10) = t18;
coef(11) = -t18;
coef(12) = t20;
coef(13) = t22;
coef(14) = -c(2);
coef(15) = c(2);
coef(16) = c(3);
coef(17) = -t29;
coef(18) = t29;
coef(19) = -t31;
coef(20) = c(5);
coef(21) = 3 * t10 * c(5);
coef(22) = -t29;
coef(23) = t31;
coef(24) = -t20;
coef(25) = -t16;
coef(26) = t29;
coef(27) = -t22;
%
t1 = 0.1e1 / AreaT;
t3 = t1 * bb(2) * cc(2);
t4 = cc(3) * bb(1);
t5 = sqrt(0.3e1);
t6 = dd(2) * t5;
t8 = t3 * t4 * t6;
t9 = t1 * cc(3);
t10 = t9 * cc(1);
t11 = bb(3) * bb(2);
t12 = dd(3) * t5;
t14 = t10 * t11 * t12;
t16 = t1 * bb(3) * cc(3);
t17 = cc(2) * bb(1);
t19 = t16 * t17 * t12;
t21 = t1 * bb(1) * cc(1);
t22 = cc(3) * bb(2);
t23 = dd(1) * t5;
t25 = t21 * t22 * t23;
t26 = t1 * cc(1);
t27 = t26 * cc(2);
t28 = bb(1) * bb(3);
t30 = t27 * t28 * t23;
t32 = t27 * t11 * t6;
t33 = cc(1) ^ 2;
t34 = t1 * t33;
t36 = bb(1) * dd(3);
t38 = t34 * bb(2) * t36 * t5;
t39 = bb(1) ^ 2;
t42 = t10 * t39 * dd(2) * t5;
t43 = bb(2) ^ 2;
t46 = t27 * t43 * dd(3) * t5;
t49 = t27 * t39 * dd(3) * t5;
t50 = cc(3) ^ 2;
t51 = t1 * t50;
t53 = bb(3) * dd(1);
t55 = t51 * bb(2) * t53 * t5;
t57 = bb(3) * dd(2);
t59 = t51 * bb(1) * t57 * t5;
t60 = cc(2) ^ 2;
t61 = t1 * t60;
t63 = bb(2) * dd(3);
t65 = t61 * bb(1) * t63 * t5;
t66 = t1 * cc(2);
t67 = t66 * cc(3);
t70 = t67 * t43 * dd(1) * t5;
t71 = bb(3) ^ 2;
t74 = t10 * t71 * dd(2) * t5;
t76 = bb(1) * dd(2);
t78 = t34 * bb(3) * t76 * t5;
t81 = t67 * t71 * dd(1) * t5;
t83 = bb(2) * dd(1);
t85 = t61 * bb(3) * t83 * t5;
t86 = -t8 - t14 + t19 + t25 - t30 + t32 - t38 - t42 - t46 + t49 + t55 - t59 + t65 + t70 + t74 + t78 - t81 - t85;
t87 = t5 * t1;
t88 = t87 * t39;
t90 = cc(2) * bb(3) * dd(1);
t91 = t88 * t90;
t93 = cc(1) * bb(2) * dd(3);
t94 = t88 * t93;
t95 = t43 * bb(2);
t98 = t87 * t95 * cc(1) * dd(3);
t99 = t71 * bb(3);
t102 = t87 * t99 * cc(1) * dd(2);
t103 = t87 * t71;
t104 = t103 * t93;
t105 = t39 * bb(1);
t108 = t87 * t105 * cc(2) * dd(3);
t110 = cc(1) * bb(3) * dd(2);
t111 = t88 * t110;
t112 = t22 * dd(1);
t113 = t103 * t112;
t114 = t87 * t43;
t115 = t114 * t90;
t116 = t17 * dd(3);
t117 = t114 * t116;
t118 = t4 * dd(2);
t119 = t114 * t118;
t120 = t103 * t116;
t121 = t88 * t112;
t124 = t87 * t95 * cc(3) * dd(1);
t125 = t103 * t118;
t128 = t87 * t99 * cc(2) * dd(1);
t131 = t87 * t105 * cc(3) * dd(2);
t132 = t114 * t110;
t133 = t91 + t94 + t98 - t102 + t104 - t108 - t111 - t113 + t115 - t117 + t119 - t120 - t121 - t124 + t125 + t128 + t131 - t132;
t135 = t1 * t99;
t138 = t135 * cc(1) * dd(2) / 0.12e2;
t139 = t1 * t95;
t142 = t139 * cc(1) * dd(3) / 0.12e2;
t143 = t50 * cc(3);
t144 = t1 * t143;
t147 = t60 * cc(2);
t148 = t1 * t147;
t151 = t1 * t105;
t154 = t151 * cc(2) * dd(3) / 0.12e2;
t157 = t135 * cc(2) * dd(1) / 0.12e2;
t158 = t33 * cc(1);
t159 = t1 * t158;
t166 = t139 * cc(3) * dd(1) / 0.12e2;
t169 = t151 * cc(3) * dd(2) / 0.12e2;
t174 = t3 * t118;
t176 = t138 - t142 - t144 * t76 / 0.6e1 + t148 * t36 / 0.6e1 + t154 - t157 - t159 * t63 / 0.6e1 + t144 * t83 / 0.6e1 + t166 - t169 - t148 * t53 / 0.6e1 + t159 * t57 / 0.6e1 + t174 / 0.4e1;
t178 = t10 * t11 * dd(3);
t180 = t16 * t116;
t183 = t27 * t11 * dd(2);
t185 = t21 * t112;
t188 = t27 * t28 * dd(1);
t192 = t1 * t71;
t194 = t192 * t118 / 0.12e2;
t196 = t192 * t93 / 0.12e2;
t197 = t1 * t43;
t199 = t197 * t110 / 0.12e2;
t201 = t197 * t116 / 0.12e2;
t208 = t28 * dd(2);
t209 = t51 * t208;
t211 = t178 / 0.4e1 - t180 / 0.4e1 - t183 / 0.4e1 - t185 / 0.4e1 + t188 / 0.4e1 + t51 * t110 / 0.6e1 - t194 - t196 + t199 + t201 + t51 * t116 / 0.6e1 - t61 * t93 / 0.6e1 - t61 * t118 / 0.6e1 + t209 / 0.4e1;
t215 = t9 * cc(1) * t71 * dd(2);
t219 = t66 * cc(1) * t43 * dd(3);
t222 = bb(1) * bb(2) * dd(3);
t223 = t61 * t222;
t225 = t1 * t39;
t227 = t225 * t93 / 0.12e2;
t229 = t225 * t90 / 0.12e2;
t231 = t192 * t116 / 0.12e2;
t233 = t192 * t112 / 0.12e2;
t242 = t34 * t222;
t246 = t26 * cc(2) * t39 * dd(3);
t248 = -t215 / 0.4e1 + t219 / 0.4e1 - t223 / 0.4e1 - t227 - t229 + t231 + t233 + t34 * t112 / 0.6e1 + t34 * t116 / 0.6e1 - t51 * t90 / 0.6e1 - t51 * t93 / 0.6e1 + t242 / 0.4e1 - t246 / 0.4e1;
t251 = t9 * cc(2) * t71 * dd(1);
t253 = t11 * dd(1);
t254 = t51 * t253;
t257 = t197 * t90 / 0.12e2;
t259 = t197 * t118 / 0.12e2;
t261 = t225 * t112 / 0.12e2;
t263 = t225 * t110 / 0.12e2;
t272 = t61 * t253;
t276 = t66 * cc(3) * t43 * dd(1);
t280 = t26 * cc(3) * t39 * dd(2);
t282 = t34 * t208;
t284 = t251 / 0.4e1 - t254 / 0.4e1 - t257 - t259 + t261 + t263 + t61 * t110 / 0.6e1 + t61 * t112 / 0.6e1 - t34 * t118 / 0.6e1 - t34 * t90 / 0.6e1 + t272 / 0.4e1 - t276 / 0.4e1 + t280 / 0.4e1 - t282 / 0.4e1;
t295 = -t138 + t142 - t154 + t157 - t166 + t169 - t174 / 0.6e1 - t178 / 0.6e1 + t180 / 0.6e1 + t183 / 0.6e1 + t185 / 0.6e1 - t188 / 0.6e1 + t194 + t196 - t199 - t201 - t209 / 0.6e1 + t215 / 0.6e1;
t306 = -t219 / 0.6e1 + t223 / 0.6e1 + t227 + t229 - t231 - t233 - t242 / 0.6e1 + t246 / 0.6e1 - t251 / 0.6e1 + t254 / 0.6e1 + t257 + t259 - t261 - t263 - t272 / 0.6e1 + t276 / 0.6e1 - t280 / 0.6e1 + t282 / 0.6e1;
t309 = -t8 - t14 + t19 + t25 - t30 + t32 - t38 - t42 - t46 + t49 + t55 - t59 + t65;
t324 = t70 / 0.18e2 + t74 / 0.18e2 + t78 / 0.18e2 - t81 / 0.18e2 - t85 / 0.18e2 + t104 / 0.72e2 - t132 / 0.72e2 - t117 / 0.72e2 + t91 / 0.72e2 - t120 / 0.72e2 - t113 / 0.72e2 + t94 / 0.72e2 - t111 / 0.72e2 + t115 / 0.72e2;
t351 = t119 / 0.72e2 - t121 / 0.72e2 + t125 / 0.72e2 - t102 / 0.72e2 + t98 / 0.72e2 - t108 / 0.72e2 + t128 / 0.72e2 - t124 / 0.72e2 + t131 / 0.72e2 - t87 * t147 * bb(1) * dd(3) / 0.18e2 + t87 * t143 * bb(1) * dd(2) / 0.18e2 - t87 * t143 * bb(2) * dd(1) / 0.18e2 - t87 * t158 * bb(3) * dd(2) / 0.18e2;
t358 = t87 * t33;
t361 = t87 * t50;
t364 = t87 * t60;
t373 = t87 * t158 * bb(2) * dd(3) + t87 * t147 * bb(3) * dd(1) + t358 * t90 + t358 * t118 - t361 * t116 + t361 * t93 - t364 * t112 + t364 * t118 - t361 * t110 + t364 * t93 + t361 * t90 - t364 * t110 - t358 * t112 - t358 * t116;
%
cm(1) = t86 / 0.12e2 + t133 / 0.24e2;
cm(2) = t176 + t211 + t248 + t284;
cm(3) = t295 + t306;
cm(4) = -t86 / 0.12e2 - t133 / 0.12e2;
cm(5) = t133 / 0.24e2;
cm(6) = t309 / 0.18e2 + t324 + t351 + t373 / 0.18e2;
%
t1 = (ko ^ 2);
t4 = (j ^ 2);
t9 = 8 / t1 / ko / t4 / j * cm(6);
t12 = 1 / t1 / t4;
t14 = 3 * t12 * cm(2);
t17 = 1 / ko / j;
t19 = 2 * t17 * cm(5);
t21 = 2 * t17 * cm(4);
t23 = 2 * t17 * cm(1);
t29 = 3 * t12 * cm(3);
t31 = 4 * t17 * cm(6);
t33 = 8 * t12 * cm(6);
%
coefm(1) = t9;
coefm(2) = -t9;
coefm(3) = -t14;
coefm(4) = -t19;
coefm(5) = -t21;
coefm(6) = t23;
coefm(7) = t9;
coefm(8) = -t23;
coefm(9) = -t9;
coefm(10) = cm(3);
coefm(11) = 3 * t17 * cm(3);
coefm(12) = 3 * t17 * cm(2);
coefm(13) = t29;
coefm(14) = cm(5);
coefm(15) = cm(4);
coefm(16) = cm(2);
coefm(17) = -t31;
coefm(18) = t19;
coefm(19) = -t33;
coefm(20) = t21;
coefm(21) = t31;
coefm(22) = -t29;
coefm(23) = t33;
coefm(24) = t14;
coefm(25) = cm(1);
coefm(26) = cm(6);
coefm(27) = -cm(6);
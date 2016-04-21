function [coef,coefm] = coefficients_g3_f1(r1,r2,r3,r4,ko,AreaT)
%% coefficients_g3_f1

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
t1 = 0.1e1 / AreaT;
t2 = bb(2) ^ 2;
t3 = t1 * t2;
t5 = bb(3) * dd(2);
t6 = sqrt(0.3e1);
t7 = t5 * t6;
t8 = t3 * cc(1) * t7;
t9 = t3 * cc(2);
t10 = bb(1) * dd(3);
t11 = t10 * t6;
t12 = t9 * t11;
t13 = bb(1) ^ 2;
t14 = t1 * t13;
t15 = t14 * cc(1);
t16 = t15 * t7;
t17 = bb(3) ^ 2;
t18 = t1 * t17;
t19 = t18 * cc(3);
t20 = bb(1) * dd(2);
t21 = t20 * t6;
t22 = t19 * t21;
t24 = bb(2) * dd(3);
t25 = t24 * t6;
t26 = t18 * cc(1) * t25;
t27 = t15 * t25;
t29 = t18 * cc(2) * t11;
t30 = bb(2) * dd(1);
t31 = t30 * t6;
t32 = t19 * t31;
t34 = t1 * t17 * bb(3);
t35 = cc(1) * dd(2);
t37 = t34 * t35 * t6;
t39 = t14 * cc(3) * t31;
t40 = cc(2) * dd(1);
t42 = t34 * t40 * t6;
t44 = t1 * t2 * bb(2);
t45 = cc(1) * dd(3);
t47 = t44 * t45 * t6;
t49 = t3 * cc(3) * t21;
t50 = bb(3) * dd(1);
t51 = t50 * t6;
t52 = t9 * t51;
t54 = t14 * cc(2) * t51;
t56 = t1 * t13 * bb(1);
t57 = cc(3) * dd(2);
t59 = t56 * t57 * t6;
t60 = cc(3) * dd(1);
t62 = t44 * t60 * t6;
t63 = cc(2) * dd(3);
t65 = t56 * t63 * t6;
t66 = -t8 - t12 - t16 + t22 + t26 + t27 - t29 - t32 - t37 - t39 + t42 + t47 + t49 + t52 + t54 + t59 - t62 - t65;
t67 = t1 * cc(3);
t68 = t67 * cc(1);
t69 = bb(3) * dd(3);
t72 = t68 * t69 * bb(2) * t6;
t75 = t1 * bb(2) * cc(2);
t76 = cc(3) * bb(1);
t79 = t75 * t76 * dd(2) * t6;
t81 = t1 * cc(2);
t82 = t81 * cc(1);
t83 = bb(2) * dd(2);
t84 = bb(3) * t6;
t86 = t82 * t83 * t84;
t88 = bb(1) * dd(1);
t90 = t82 * t88 * t84;
t92 = t67 * cc(2);
t95 = t92 * t69 * bb(1) * t6;
t98 = t1 * bb(1) * cc(1);
t99 = cc(3) * bb(2);
t102 = t98 * t99 * dd(1) * t6;
t111 = -t72 / 0.18e2 - t79 / 0.18e2 + t86 / 0.18e2 - t90 / 0.18e2 + t95 / 0.18e2 + t102 / 0.18e2 - t8 / 0.72e2 - t12 / 0.72e2 - t16 / 0.72e2 + t22 / 0.72e2 + t26 / 0.72e2 + t27 / 0.72e2 - t29 / 0.72e2;
t117 = cc(3) ^ 2;
t118 = t1 * t117;
t120 = t118 * bb(1) * t7;
t122 = cc(2) ^ 2;
t123 = t1 * t122;
t125 = t123 * bb(3) * t31;
t128 = t123 * bb(1) * t25;
t132 = t68 * t17 * dd(2) * t6;
t136 = t82 * t2 * dd(3) * t6;
t140 = t92 * t17 * dd(1) * t6;
t143 = t118 * bb(2) * t51;
t147 = t92 * t2 * dd(1) * t6;
t151 = t68 * t13 * dd(2) * t6;
t153 = -t32 / 0.72e2 - t39 / 0.72e2 + t49 / 0.72e2 + t52 / 0.72e2 + t54 / 0.72e2 - t120 / 0.18e2 - t125 / 0.18e2 + t128 / 0.18e2 + t132 / 0.18e2 - t136 / 0.18e2 - t140 / 0.18e2 + t143 / 0.18e2 + t147 / 0.18e2 - t151 / 0.18e2;
t155 = cc(1) ^ 2;
t156 = t1 * t155;
t158 = t156 * bb(3) * t21;
t160 = t156 * bb(2) * t11;
t163 = t82 * t13 * dd(3) * t6;
t164 = t6 * t1;
t165 = t164 * t122;
t166 = t76 * dd(2);
t169 = cc(1) * bb(2) * dd(3);
t171 = t164 * t155;
t173 = bb(3) * cc(2) * dd(1);
t176 = bb(3) * cc(1) * dd(2);
t178 = t164 * t117;
t181 = bb(1) * cc(2) * dd(3);
t183 = t99 * dd(1);
t188 = t158 - t160 + t163 + t165 * t166 + t165 * t169 + t171 * t173 - t165 * t176 + t178 * t173 - t178 * t181 - t171 * t183 + t171 * t166 - t171 * t181 + t178 * t169;
t199 = t122 * cc(2);
t208 = t155 * cc(1);
t213 = t117 * cc(3);
t226 = -t178 * t176 / 0.18e2 - t165 * t183 / 0.18e2 - t37 / 0.72e2 + t42 / 0.72e2 + t47 / 0.72e2 + t59 / 0.72e2 - t62 / 0.72e2 - t65 / 0.72e2 + t164 * t199 * bb(3) * dd(1) / 0.18e2 - t164 * t199 * bb(1) * dd(3) / 0.18e2 + t164 * t208 * bb(2) * dd(3) / 0.18e2 - t164 * t213 * bb(2) * dd(1) / 0.18e2 + t164 * t213 * bb(1) * dd(2) / 0.18e2 - t164 * t208 * bb(3) * dd(2) / 0.18e2;
t229 = -t72 - t140 + t128 + t132 - t125 + t102 + t143 - t136 + t147 - t151 + t158 - t79 - t120 - t160 + t163 + t86 - t90 + t95;
t248 = t72 / 0.12e2 + t79 / 0.12e2 - t86 / 0.12e2 + t90 / 0.12e2 - t95 / 0.12e2 - t102 / 0.12e2 + t8 / 0.24e2 + t12 / 0.24e2 + t16 / 0.24e2 - t22 / 0.24e2 - t26 / 0.24e2 - t27 / 0.24e2 + t29 / 0.24e2 + t32 / 0.24e2 + t39 / 0.24e2 - t49 / 0.24e2 - t52 / 0.24e2 - t54 / 0.24e2;
t267 = t120 / 0.12e2 + t125 / 0.12e2 - t128 / 0.12e2 - t132 / 0.12e2 + t136 / 0.12e2 + t140 / 0.12e2 - t143 / 0.12e2 - t147 / 0.12e2 + t151 / 0.12e2 - t158 / 0.12e2 + t160 / 0.12e2 - t163 / 0.12e2 + t37 / 0.24e2 - t42 / 0.24e2 - t47 / 0.24e2 - t59 / 0.24e2 + t62 / 0.24e2 + t65 / 0.24e2;
t269 = t1 * t213;
t272 = t1 * t199;
t275 = t1 * t208;
t285 = t82 * t88 * bb(3);
t288 = t92 * t69 * bb(1);
t291 = t68 * t69 * bb(2);
t294 = t82 * t83 * bb(3);
t296 = t75 * t166;
t298 = t98 * t183;
t300 = t50 * bb(2);
t301 = t123 * t300;
t305 = t81 * cc(3) * t2 * dd(1);
t307 = t1 * cc(1);
t310 = t307 * cc(3) * t13 * dd(2);
t312 = t5 * bb(1);
t313 = t156 * t312;
t315 = t118 * t300;
t319 = -t269 * t20 / 0.6e1 + t272 * t10 / 0.6e1 - t275 * t24 / 0.6e1 + t269 * t30 / 0.6e1 - t272 * t50 / 0.6e1 + t275 * t5 / 0.6e1 + t285 / 0.12e2 - t288 / 0.12e2 + t291 / 0.12e2 - t294 / 0.12e2 + t296 / 0.12e2 - t298 / 0.12e2 + t301 / 0.12e2 - t305 / 0.12e2 + t310 / 0.12e2 - t313 / 0.12e2 - t315 / 0.12e2 - t118 * t169 / 0.6e1;
t322 = t24 * bb(1);
t323 = t156 * t322;
t329 = t307 * cc(2) * t13 * dd(3);
t335 = t67 * cc(2) * t17 * dd(1);
t339 = t118 * t312;
t345 = t67 * cc(1) * t17 * dd(2);
t351 = t81 * cc(1) * t2 * dd(3);
t355 = t123 * t322;
t365 = t156 * t183 / 0.6e1 + t323 / 0.12e2 + t156 * t181 / 0.6e1 - t329 / 0.12e2 - t118 * t173 / 0.6e1 + t335 / 0.12e2 + t118 * t181 / 0.6e1 + t339 / 0.12e2 + t118 * t176 / 0.6e1 - t345 / 0.12e2 - t123 * t169 / 0.6e1 + t351 / 0.12e2 - t123 * t166 / 0.6e1 - t355 / 0.12e2 + t123 * t176 / 0.6e1 + t123 * t183 / 0.6e1 - t156 * t166 / 0.6e1 - t156 * t173 / 0.6e1;
t391 = t34 * t35 / 0.12e2 - t44 * t45 / 0.12e2 + t56 * t63 / 0.12e2 - t34 * t40 / 0.12e2 + t44 * t60 / 0.12e2 - t56 * t57 / 0.12e2 + t285 / 0.6e1 - t288 / 0.6e1 + t291 / 0.6e1 - t294 / 0.6e1 + t296 / 0.6e1 - t298 / 0.6e1 + t301 / 0.6e1 - t305 / 0.6e1 + t310 / 0.6e1 - t313 / 0.6e1 - t315 / 0.6e1 + t323 / 0.6e1;
t422 = -t329 / 0.6e1 + t335 / 0.6e1 + t339 / 0.6e1 - t345 / 0.6e1 + t351 / 0.6e1 - t355 / 0.6e1 - t3 * t173 / 0.12e2 + t18 * t183 / 0.12e2 - t18 * t166 / 0.12e2 - t18 * t169 / 0.12e2 + t3 * t176 / 0.12e2 + t3 * t181 / 0.12e2 - t14 * t169 / 0.12e2 - t14 * t173 / 0.12e2 + t18 * t181 / 0.12e2 - t3 * t166 / 0.12e2 + t14 * t183 / 0.12e2 + t14 * t176 / 0.12e2;
%
c(1) = t66 / 0.24e2;
c(2) = t111 + t153 + t188 / 0.18e2 + t226;
c(3) = t229 / 0.12e2;
c(4) = t248 + t267;
c(5) = t319 + t365;
c(6) = t391 + t422;
%
t3 = (0.1e1 / ko / j);
t6 = (ko ^ 2);
t8 = (j ^ 2);
t10 = (0.1e1 / t6 / t8);
t12 = 3 * t10 * c(6);
t21 = 8 / t6 / ko / t8 / j * c(2);
t23 = 2 * t3 * c(1);
t25 = 2 * t3 * c(4);
t27 = 3 * t10 * c(5);
t29 = 2 * t3 * c(3);
t31 = 8 * t10 * c(2);
t33 = 4 * t3 * c(2);
%
coef(1) = c(5);
coef(2) = c(3);
coef(3) = c(1);
coef(4) = 3 * t3 * c(5);
coef(5) = t12;
coef(6) = 3 * t3 * c(6);
coef(7) = c(6);
coef(8) = t21;
coef(9) = -t21;
coef(10) = -t23;
coef(11) = t25;
coef(12) = -t21;
coef(13) = t21;
coef(14) = -t27;
coef(15) = -t29;
coef(16) = c(4);
coef(17) = -c(2);
coef(18) = c(2);
coef(19) = -t31;
coef(20) = t33;
coef(21) = -t33;
coef(22) = t23;
coef(23) = -t12;
coef(24) = t29;
coef(25) = t27;
coef(26) = t31;
coef(27) = -t25;
%
t1 = 0.1e1 / AreaT;
t2 = t1 * bb(1);
t3 = cc(1) ^ 2;
t4 = t2 * t3;
t5 = bb(2) * dd(3);
t6 = sqrt(0.3e1);
t7 = t5 * t6;
t8 = t4 * t7;
t10 = bb(3) ^ 2;
t11 = t1 * t10;
t12 = t11 * cc(2);
t13 = cc(3) * dd(1);
t14 = t13 * t6;
t15 = t12 * t14;
t17 = bb(1) ^ 2;
t18 = t1 * t17;
t19 = t18 * cc(2);
t20 = cc(1) * dd(3);
t21 = t20 * t6;
t22 = t19 * t21;
t24 = t1 * bb(2);
t25 = cc(2) ^ 2;
t26 = t24 * t25;
t27 = bb(3) * dd(1);
t28 = t27 * t6;
t29 = t26 * t28;
t31 = t1 * bb(3);
t32 = cc(3) ^ 2;
t33 = t31 * t32;
t34 = bb(2) * dd(1);
t35 = t34 * t6;
t36 = t33 * t35;
t38 = t11 * cc(1);
t39 = cc(3) * dd(2);
t40 = t39 * t6;
t41 = t38 * t40;
t43 = bb(1) * dd(2);
t44 = t43 * t6;
t45 = t33 * t44;
t47 = t18 * cc(3);
t48 = cc(1) * dd(2);
t49 = t48 * t6;
t50 = t47 * t49;
t52 = bb(1) * dd(3);
t53 = t52 * t6;
t54 = t26 * t53;
t56 = t11 * cc(3);
t57 = t56 * t44;
t59 = t38 * t7;
t61 = bb(2) ^ 2;
t62 = t1 * t61;
t63 = t62 * cc(1);
t64 = bb(3) * dd(2);
t65 = t64 * t6;
t66 = t63 * t65;
t68 = t62 * cc(2);
t69 = t68 * t53;
t71 = t18 * cc(1);
t72 = t71 * t7;
t74 = t19 * t28;
t76 = t12 * t53;
t78 = t56 * t35;
t80 = t68 * t28;
t82 = t8 / 0.12e2 + t15 / 0.12e2 - t22 / 0.12e2 + t29 / 0.12e2 - t36 / 0.12e2 - t41 / 0.12e2 + t45 / 0.12e2 + t50 / 0.12e2 - t54 / 0.12e2 - t57 / 0.24e2 - t59 / 0.24e2 + t66 / 0.24e2 + t69 / 0.24e2 - t72 / 0.24e2 - t74 / 0.24e2 + t76 / 0.24e2 + t78 / 0.24e2 - t80 / 0.24e2;
t83 = t62 * cc(3);
t84 = t83 * t44;
t86 = t47 * t35;
t88 = cc(2) * dd(1);
t89 = t88 * t6;
t90 = t83 * t89;
t92 = cc(2) * dd(3);
t93 = t92 * t6;
t94 = t63 * t93;
t96 = t4 * t65;
t98 = t71 * t65;
t101 = t1 * t10 * bb(3);
t102 = t101 * t49;
t105 = t1 * t61 * bb(2);
t106 = t105 * t21;
t109 = t1 * t17 * bb(1);
t110 = t109 * t93;
t112 = t101 * t89;
t114 = t105 * t14;
t116 = t109 * t40;
t119 = t1 * cc(2) * cc(1);
t120 = bb(2) * bb(3);
t121 = dd(2) * t6;
t123 = t119 * t120 * t121;
t125 = bb(1) * bb(3);
t126 = dd(1) * t6;
t128 = t119 * t125 * t126;
t131 = t1 * cc(3) * cc(1);
t132 = dd(3) * t6;
t134 = t131 * t120 * t132;
t136 = t24 * cc(2);
t137 = cc(3) * bb(1);
t139 = t136 * t137 * t121;
t141 = t2 * cc(1);
t142 = cc(3) * bb(2);
t144 = t141 * t142 * t126;
t146 = t31 * cc(3);
t147 = cc(2) * bb(1);
t149 = t146 * t147 * t132;
t151 = -t84 / 0.24e2 + t86 / 0.24e2 - t90 / 0.12e2 + t94 / 0.12e2 - t96 / 0.12e2 + t98 / 0.24e2 + t102 / 0.24e2 - t106 / 0.24e2 + t110 / 0.24e2 - t112 / 0.24e2 + t114 / 0.24e2 - t116 / 0.24e2 - t123 / 0.12e2 + t128 / 0.12e2 + t134 / 0.12e2 + t139 / 0.12e2 - t144 / 0.12e2 - t149 / 0.12e2;
t154 = t131 * t120 * dd(3);
t157 = t119 * t120 * dd(2);
t159 = t142 * dd(1);
t160 = t141 * t159;
t163 = t119 * t125 * dd(1);
t165 = t137 * dd(2);
t166 = t136 * t165;
t168 = t147 * dd(3);
t169 = t146 * t168;
t171 = t32 * cc(3);
t172 = t1 * t171;
t175 = t25 * cc(2);
t176 = t1 * t175;
t179 = t3 * cc(1);
t180 = t1 * t179;
t191 = t2 * t3 * bb(2) * dd(3);
t194 = cc(2) * cc(1) * dd(3);
t195 = t18 * t194;
t198 = cc(2) * cc(3) * dd(1);
t199 = t11 * t198;
t203 = t31 * t32 * bb(2) * dd(1);
t207 = t24 * t25 * bb(3) * dd(1);
t209 = t62 * t198;
t211 = t154 / 0.12e2 - t157 / 0.12e2 - t160 / 0.12e2 + t163 / 0.12e2 + t166 / 0.12e2 - t169 / 0.12e2 - t172 * t43 / 0.6e1 + t176 * t52 / 0.6e1 - t180 * t5 / 0.6e1 + t172 * t34 / 0.6e1 - t176 * t27 / 0.6e1 + t180 * t64 / 0.6e1 + t191 / 0.12e2 - t195 / 0.12e2 + t199 / 0.12e2 - t203 / 0.12e2 + t207 / 0.12e2 - t209 / 0.12e2;
t213 = cc(3) * cc(1) * dd(2);
t214 = t18 * t213;
t218 = t2 * t3 * bb(3) * dd(2);
t220 = t1 * t32;
t224 = cc(1) * bb(3) * dd(2);
t227 = t1 * t25;
t229 = cc(1) * bb(2) * dd(3);
t234 = t1 * t3;
t240 = cc(2) * bb(3) * dd(1);
t255 = t24 * t25 * bb(1) * dd(3);
t257 = t62 * t194;
t261 = t31 * t32 * bb(1) * dd(2);
t263 = t11 * t213;
t265 = t214 / 0.12e2 - t218 / 0.12e2 + t220 * t168 / 0.6e1 + t220 * t224 / 0.6e1 - t227 * t229 / 0.6e1 - t227 * t165 / 0.6e1 + t234 * t159 / 0.6e1 + t234 * t168 / 0.6e1 - t220 * t240 / 0.6e1 - t220 * t229 / 0.6e1 + t227 * t224 / 0.6e1 + t227 * t159 / 0.6e1 - t234 * t165 / 0.6e1 - t234 * t240 / 0.6e1 - t255 / 0.12e2 + t257 / 0.12e2 + t261 / 0.12e2 - t263 / 0.12e2;
t291 = -t154 / 0.6e1 + t157 / 0.6e1 + t160 / 0.6e1 - t163 / 0.6e1 - t166 / 0.6e1 + t169 / 0.6e1 - t101 * t48 / 0.12e2 + t105 * t20 / 0.12e2 - t109 * t92 / 0.12e2 + t101 * t88 / 0.12e2 - t105 * t13 / 0.12e2 + t109 * t39 / 0.12e2 - t191 / 0.6e1 + t195 / 0.6e1 - t199 / 0.6e1 + t203 / 0.6e1 - t207 / 0.6e1 + t209 / 0.6e1;
t322 = -t214 / 0.6e1 + t218 / 0.6e1 + t255 / 0.6e1 - t257 / 0.6e1 - t261 / 0.6e1 + t263 / 0.6e1 + t11 * t165 / 0.12e2 + t11 * t229 / 0.12e2 - t62 * t224 / 0.12e2 - t62 * t168 / 0.12e2 + t18 * t229 / 0.12e2 + t18 * t240 / 0.12e2 - t11 * t168 / 0.12e2 - t11 * t159 / 0.12e2 + t62 * t240 / 0.12e2 + t62 * t165 / 0.12e2 - t18 * t159 / 0.12e2 - t18 * t224 / 0.12e2;
t324 = -t86 + t84 - t78 - t114 + t59 + t80 - t98 + t106 + t74 - t102 - t69 + t72 + t57 - t110 - t66 + t116 - t76 + t112;
t338 = -t8 / 0.18e2 - t15 / 0.18e2 + t22 / 0.18e2 - t29 / 0.18e2 + t36 / 0.18e2 + t41 / 0.18e2 - t45 / 0.18e2 - t50 / 0.18e2 + t54 / 0.18e2 + t57 / 0.72e2 + t59 / 0.72e2 - t66 / 0.72e2 - t69 / 0.72e2;
t350 = t6 * t1;
t351 = t350 * t3;
t356 = t350 * t32;
t359 = t72 / 0.72e2 + t74 / 0.72e2 - t76 / 0.72e2 - t78 / 0.72e2 + t80 / 0.72e2 + t84 / 0.72e2 - t86 / 0.72e2 + t90 / 0.18e2 - t94 / 0.18e2 + t96 / 0.18e2 - t98 / 0.72e2 + t351 * t240 / 0.18e2 + t351 * t165 / 0.18e2 - t356 * t224 / 0.18e2;
t361 = t350 * t25;
t384 = t361 * t165 / 0.18e2 - t361 * t159 / 0.18e2 - t351 * t168 / 0.18e2 - t361 * t224 / 0.18e2 + t356 * t240 / 0.18e2 - t351 * t159 / 0.18e2 + t361 * t229 / 0.18e2 + t356 * t229 / 0.18e2 - t356 * t168 / 0.18e2 - t102 / 0.72e2 + t106 / 0.72e2 - t110 / 0.72e2 + t112 / 0.72e2;
t417 = -t114 / 0.72e2 + t116 / 0.72e2 + t350 * t179 * bb(2) * dd(3) / 0.18e2 + t350 * t171 * bb(1) * dd(2) / 0.18e2 - t350 * t171 * bb(2) * dd(1) / 0.18e2 - t350 * t175 * bb(1) * dd(3) / 0.18e2 - t350 * t179 * bb(3) * dd(2) / 0.18e2 + t350 * t175 * bb(3) * dd(1) / 0.18e2 + t123 / 0.18e2 - t128 / 0.18e2 - t134 / 0.18e2 - t139 / 0.18e2 + t144 / 0.18e2 + t149 / 0.18e2;
t420 = t8 + t15 + t29 - t36 + t94 - t96 - t54 - t22 + t134 + t139 - t144 + t128 - t123 - t41 + t45 - t149 + t50 - t90;
%
cm(1) = t82 + t151;
cm(2) = t211 + t265;
cm(3) = t291 + t322;
cm(4) = t324 / 0.24e2;
cm(5) = t338 + t359 + t384 + t417;
cm(6) = t420 / 0.12e2;
%
t1 = (ko ^ 2);
t3 = (j ^ 2);
t5 = (0.1e1 / t1 / t3);
t7 = 3 * t5 * cm(3);
t10 = (0.1e1 / ko / j);
t19 = 8 / t1 / ko / t3 / j * cm(5);
t21 = 2 * t10 * cm(1);
t25 = 3 * t5 * cm(2);
t27 = 2 * t10 * cm(6);
t29 = 2 * t10 * cm(4);
t31 = 8 * t5 * cm(5);
t33 = 4 * t10 * cm(5);
%
coefm(1) = cm(2);
coefm(2) = cm(6);
coefm(3) = cm(4);
coefm(4) = t7;
coefm(5) = 3 * t10 * cm(2);
coefm(6) = t19;
coefm(7) = -t19;
coefm(8) = -t21;
coefm(9) = cm(3);
coefm(10) = 3 * t10 * cm(3);
coefm(11) = t19;
coefm(12) = -t19;
coefm(13) = -t25;
coefm(14) = -t27;
coefm(15) = -t29;
coefm(16) = t21;
coefm(17) = cm(5);
coefm(18) = cm(1);
coefm(19) = -cm(5);
coefm(20) = t31;
coefm(21) = t25;
coefm(22) = -t7;
coefm(23) = -t31;
coefm(24) = t33;
coefm(25) = t29;
coefm(26) = t27;
coefm(27) = -t33;
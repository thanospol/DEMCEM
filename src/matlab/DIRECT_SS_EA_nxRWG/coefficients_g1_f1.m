function [coef,coefm] = coefficients_g1_f1(r1,r2,r3,r4,ko,AreaT)
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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% ko = wavenumber
% AreaT = Area of outer triangle

% OUTPUT DATA
% coef   
% coefm
%

c  = zeros(4);
cm = zeros(4);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
% 
t1 = 0.1e1 / AreaT;
t2 = t1 * cc(2);
t3 = t2 * cc(1);
t4 = bb(2) ^ 2;
t6 = sqrt(0.3e1);
t10 = cc(3) ^ 2;
t11 = t1 * t10;
t13 = bb(3) * dd(2);
t14 = t13 * t6;
t17 = cc(1) ^ 2;
t18 = t1 * t17;
t20 = bb(1) * dd(3);
t21 = t20 * t6;
t24 = bb(1) ^ 2;
t29 = t1 * cc(3);
t30 = t29 * cc(2);
t31 = bb(3) ^ 2;
t36 = t29 * cc(1);
t41 = cc(2) ^ 2;
t42 = t1 * t41;
t45 = dd(3) * bb(2) * t6;
t49 = bb(3) * dd(1);
t50 = t49 * t6;
t63 = dd(2) * bb(1) * t6;
t68 = dd(1) * bb(2) * t6;
t71 = t1 * t4;
t73 = t71 * cc(3) * t63;
t75 = -t3 * t4 * dd(3) * t6 / 0.18e2 - t11 * bb(1) * t14 / 0.18e2 - t18 * bb(2) * t21 / 0.18e2 + t3 * t24 * dd(3) * t6 / 0.18e2 - t30 * t31 * dd(1) * t6 / 0.18e2 + t36 * t31 * dd(2) * t6 / 0.18e2 + t42 * bb(1) * t45 / 0.18e2 + t11 * bb(2) * t50 / 0.18e2 + t30 * t4 * dd(1) * t6 / 0.18e2 - t36 * t24 * dd(2) * t6 / 0.18e2 + t18 * bb(3) * t63 / 0.18e2 - t42 * bb(3) * t68 / 0.18e2 + t73 / 0.72e2;
t76 = t1 * t24;
t78 = t76 * cc(3) * t68;
t80 = t76 * cc(1);
t81 = t80 * t14;
t83 = t71 * cc(2);
t84 = t83 * t50;
t86 = t83 * t21;
t88 = t1 * t31;
t90 = t88 * cc(1) * t45;
t92 = t80 * t45;
t94 = t88 * cc(3);
t95 = t94 * t63;
t98 = t76 * cc(2) * t50;
t101 = t88 * cc(2) * t21;
t104 = t71 * cc(1) * t14;
t106 = t94 * t68;
t108 = t6 * t1;
t109 = t108 * t17;
t111 = bb(3) * cc(2) * dd(1);
t114 = cc(2) * bb(1);
t115 = t114 * dd(3);
t118 = t108 * t41;
t119 = cc(3) * bb(1);
t120 = t119 * dd(2);
t123 = -t78 / 0.72e2 - t81 / 0.72e2 + t84 / 0.72e2 - t86 / 0.72e2 + t90 / 0.72e2 + t92 / 0.72e2 + t95 / 0.72e2 + t98 / 0.72e2 - t101 / 0.72e2 - t104 / 0.72e2 - t106 / 0.72e2 + t109 * t111 / 0.18e2 - t109 * t115 / 0.18e2 + t118 * t120 / 0.18e2;
t126 = bb(3) * cc(1) * dd(2);
t129 = cc(1) * bb(2) * dd(3);
t131 = cc(3) * bb(2);
t132 = t131 * dd(1);
t134 = t108 * t10;
t142 = t1 * bb(3) * cc(3);
t143 = dd(3) * t6;
t146 = bb(2) * bb(3);
t147 = dd(2) * t6;
t151 = t1 * bb(2) * cc(2);
t156 = -t118 * t126 + t118 * t129 - t118 * t132 + t134 * t111 - t109 * t132 + t109 * t120 + t134 * t129 - t134 * t126 - t134 * t115 + t142 * t114 * t143 + t3 * t146 * t147 - t151 * t119 * t147 - t36 * t146 * t143;
t157 = bb(3) * bb(1);
t158 = dd(1) * t6;
t163 = t1 * bb(1) * cc(1);
t168 = t1 * t24 * bb(1);
t169 = cc(2) * dd(3);
t171 = t168 * t169 * t6;
t174 = t1 * t31 * bb(3);
t175 = cc(2) * dd(1);
t177 = t174 * t175 * t6;
t180 = t1 * t4 * bb(2);
t181 = cc(3) * dd(1);
t183 = t180 * t181 * t6;
t185 = cc(1) * dd(2);
t187 = t174 * t185 * t6;
t189 = cc(3) * dd(2);
t191 = t168 * t189 * t6;
t193 = cc(1) * dd(3);
t195 = t180 * t193 * t6;
t197 = t10 * cc(3);
t206 = t41 * cc(2);
t211 = t17 * cc(1);
t224 = -t3 * t157 * t158 / 0.18e2 + t163 * t131 * t158 / 0.18e2 - t171 / 0.72e2 + t177 / 0.72e2 - t183 / 0.72e2 - t187 / 0.72e2 + t191 / 0.72e2 + t195 / 0.72e2 + t108 * t197 * bb(1) * dd(2) / 0.18e2 - t108 * t197 * bb(2) * dd(1) / 0.18e2 + t108 * t206 * bb(3) * dd(1) / 0.18e2 - t108 * t211 * bb(3) * dd(2) / 0.18e2 - t108 * t206 * bb(1) * dd(3) / 0.18e2 + t108 * t211 * bb(2) * dd(3) / 0.18e2;
t227 = t171 - t98 + t86 + t187 + t106 + t101 + t183 - t177 - t90 - t195 + t104 + t81 - t92 - t191 - t73 - t84 - t95 + t78;
t261 = t49 * bb(2);
t268 = t1 * cc(1);
t273 = t174 * t185 / 0.12e2 - t180 * t193 / 0.12e2 + t168 * t169 / 0.12e2 - t174 * t175 / 0.12e2 + t180 * t181 / 0.12e2 - t168 * t189 / 0.12e2 - t142 * t115 / 0.6e1 + t151 * t120 / 0.6e1 + t36 * t146 * dd(3) / 0.6e1 - t3 * t146 * dd(2) / 0.6e1 - t163 * t132 / 0.6e1 + t3 * t157 * dd(1) / 0.6e1 - t71 * t120 / 0.12e2 + t76 * t132 / 0.12e2 + t76 * t126 / 0.12e2 + t42 * t261 / 0.6e1 - t2 * cc(3) * t4 * dd(1) / 0.6e1 + t268 * cc(3) * t24 * dd(2) / 0.6e1;
t274 = t13 * bb(1);
t295 = t20 * bb(2);
t320 = -t18 * t274 / 0.6e1 - t88 * t120 / 0.12e2 - t88 * t129 / 0.12e2 + t71 * t126 / 0.12e2 + t71 * t115 / 0.12e2 + t11 * t274 / 0.6e1 - t29 * cc(1) * t31 * dd(2) / 0.6e1 + t2 * cc(1) * t4 * dd(3) / 0.6e1 - t42 * t295 / 0.6e1 - t76 * t129 / 0.12e2 - t76 * t111 / 0.12e2 + t88 * t115 / 0.12e2 + t88 * t132 / 0.12e2 + t18 * t295 / 0.6e1 - t268 * cc(2) * t24 * dd(3) / 0.6e1 + t29 * cc(2) * t31 * dd(1) / 0.6e1 - t11 * t261 / 0.6e1 - t71 * t111 / 0.12e2;
%
c(1) = t75 + t123 + t156 / 0.18e2 + t224;
c(2) = t227 / 0.24e2;
c(3) = t273 + t320;
c(4) = -t227 / 0.24e2;
%
t3 = (0.1e1 / ko / j);
t6 = (ko ^ 2);
t9 = (j ^ 2);
t14 = 8 / t6 / ko / t9 / j * c(1);
t16 = 2 * t3 * c(2);
t18 = 2 * t3 * c(4);
t21 = 1 / t6 / t9;
t23 = 3 * t21 * c(3);
t25 = 8 * t21 * c(1);
t27 = 4 * t3 * c(1);
%
coef(1) = c(3);
coef(2) = 3 * t3 * c(3);
coef(3) = -t14;
coef(4) = t14;
coef(5) = t16;
coef(6) = -t18;
coef(7) = -t14;
coef(8) = -t16;
coef(9) = t14;
coef(10) = c(4);
coef(11) = t23;
coef(12) = -c(1);
coef(13) = c(1);
coef(14) = c(2);
coef(15) = -t23;
coef(16) = t18;
coef(17) = t25;
coef(18) = -t25;
coef(19) = t27;
coef(20) = -t27;
%
t1 = 0.1e1 / AreaT;
t2 = bb(3) ^ 2;
t4 = t1 * t2 * bb(3);
t5 = cc(2) * dd(1);
t6 = sqrt(0.3e1);
t7 = t5 * t6;
t8 = t4 * t7;
t9 = t1 * t2;
t10 = t9 * cc(1);
t12 = bb(2) * dd(3) * t6;
t13 = t10 * t12;
t14 = bb(1) ^ 2;
t15 = t1 * t14;
t16 = t15 * cc(1);
t17 = t16 * t12;
t18 = bb(2) ^ 2;
t20 = t1 * t18 * bb(2);
t21 = cc(3) * dd(1);
t22 = t21 * t6;
t23 = t20 * t22;
t24 = t1 * t18;
t25 = t24 * cc(1);
t27 = bb(3) * dd(2) * t6;
t28 = t25 * t27;
t30 = t1 * t14 * bb(1);
t31 = cc(2) * dd(3);
t32 = t31 * t6;
t33 = t30 * t32;
t34 = t9 * cc(3);
t36 = bb(2) * dd(1) * t6;
t37 = t34 * t36;
t38 = t16 * t27;
t39 = t9 * cc(2);
t41 = bb(1) * dd(3) * t6;
t42 = t39 * t41;
t43 = t24 * cc(2);
t45 = bb(3) * dd(1) * t6;
t46 = t43 * t45;
t47 = t24 * cc(3);
t49 = bb(1) * dd(2) * t6;
t50 = t47 * t49;
t51 = t43 * t41;
t52 = t15 * cc(2);
t53 = t52 * t45;
t54 = t34 * t49;
t55 = cc(1) * dd(3);
t56 = t55 * t6;
t57 = t20 * t56;
t58 = t15 * cc(3);
t59 = t58 * t36;
t60 = cc(3) * dd(2);
t61 = t60 * t6;
t62 = t30 * t61;
t63 = cc(1) * dd(2);
t64 = t63 * t6;
t65 = t4 * t64;
t66 = t8 + t13 + t17 - t23 - t28 - t33 - t37 - t38 - t42 + t46 + t50 - t51 + t53 + t54 + t57 - t59 + t62 - t65;
t79 = t1 * bb(2);
t80 = t79 * cc(2);
t81 = cc(3) * bb(1);
t82 = t81 * dd(2);
t85 = t1 * bb(3);
t86 = t85 * cc(3);
t87 = cc(2) * bb(1);
t88 = t87 * dd(3);
t92 = t1 * cc(3) * cc(1);
t93 = bb(3) * bb(2);
t98 = t1 * cc(2) * cc(1);
t102 = t1 * bb(1);
t103 = t102 * cc(1);
t104 = cc(3) * bb(2);
t105 = t104 * dd(1);
t108 = bb(1) * bb(3);
t112 = cc(2) ^ 2;
t120 = cc(1) * cc(3) * dd(2);
t123 = cc(3) ^ 2;
t129 = cc(1) * cc(2) * dd(3);
t132 = cc(1) ^ 2;
t137 = -t4 * t63 / 0.12e2 + t20 * t55 / 0.12e2 - t30 * t31 / 0.12e2 + t4 * t5 / 0.12e2 - t20 * t21 / 0.12e2 + t30 * t60 / 0.12e2 - t80 * t82 / 0.6e1 + t86 * t88 / 0.6e1 - t92 * t93 * dd(3) / 0.6e1 + t98 * t93 * dd(2) / 0.6e1 + t103 * t105 / 0.6e1 - t98 * t108 * dd(1) / 0.6e1 + t79 * t112 * bb(1) * dd(3) / 0.6e1 - t24 * t88 / 0.12e2 + t9 * t120 / 0.6e1 - t85 * t123 * bb(1) * dd(2) / 0.6e1 - t24 * t129 / 0.6e1 - t102 * t132 * bb(2) * dd(3) / 0.6e1;
t141 = cc(2) * cc(3) * dd(1);
t161 = cc(1) * bb(3) * dd(2);
t165 = cc(1) * bb(2) * dd(3);
t173 = cc(2) * bb(3) * dd(1);
t188 = t15 * t129 / 0.6e1 - t9 * t141 / 0.6e1 + t85 * t123 * bb(2) * dd(1) / 0.6e1 - t79 * t112 * bb(3) * dd(1) / 0.6e1 + t24 * t141 / 0.6e1 - t15 * t120 / 0.6e1 + t102 * t132 * bb(3) * dd(2) / 0.6e1 - t24 * t161 / 0.12e2 + t9 * t165 / 0.12e2 + t9 * t82 / 0.12e2 + t15 * t165 / 0.12e2 + t15 * t173 / 0.12e2 - t9 * t88 / 0.12e2 - t9 * t105 / 0.12e2 + t24 * t173 / 0.12e2 + t24 * t82 / 0.12e2 - t15 * t161 / 0.12e2 - t15 * t105 / 0.12e2;
t204 = t13 / 0.72e2 + t17 / 0.72e2 - t28 / 0.72e2 - t37 / 0.72e2 - t38 / 0.72e2 - t42 / 0.72e2 + t46 / 0.72e2 + t50 / 0.72e2 - t51 / 0.72e2 + t53 / 0.72e2 + t54 / 0.72e2 - t59 / 0.72e2 + t10 * t61 / 0.18e2;
t205 = t85 * t123;
t207 = t102 * t132;
t211 = t79 * t112;
t219 = t6 * t1;
t220 = t219 * t123;
t222 = t219 * t132;
t224 = t219 * t112;
t226 = -t205 * t49 + t207 * t27 + t205 * t36 + t47 * t7 + t211 * t41 - t207 * t12 - t39 * t22 - t25 * t32 - t58 * t64 + t52 * t56 - t211 * t45 + t220 * t173 - t222 * t88 - t224 * t105;
t250 = -t220 * t88 / 0.18e2 + t224 * t165 / 0.18e2 + t222 * t82 / 0.18e2 - t220 * t161 / 0.18e2 + t222 * t173 / 0.18e2 + t220 * t165 / 0.18e2 - t222 * t105 / 0.18e2 - t224 * t161 / 0.18e2 + t224 * t82 / 0.18e2 + t8 / 0.72e2 - t23 / 0.72e2 - t33 / 0.72e2 + t57 / 0.72e2;
t253 = t123 * cc(3);
t258 = t132 * cc(1);
t263 = t112 * cc(2);
t280 = dd(1) * t6;
t284 = dd(3) * t6;
t288 = dd(2) * t6;
t301 = t62 / 0.72e2 - t65 / 0.72e2 + t219 * t253 * bb(1) * dd(2) / 0.18e2 + t219 * t258 * bb(2) * dd(3) / 0.18e2 + t219 * t263 * bb(3) * dd(1) / 0.18e2 - t219 * t253 * bb(2) * dd(1) / 0.18e2 - t219 * t263 * bb(1) * dd(3) / 0.18e2 - t219 * t258 * bb(3) * dd(2) / 0.18e2 + t103 * t104 * t280 / 0.18e2 - t92 * t93 * t284 / 0.18e2 + t98 * t93 * t288 / 0.18e2 - t80 * t81 * t288 / 0.18e2 - t98 * t108 * t280 / 0.18e2 + t86 * t87 * t284 / 0.18e2;
%
cm(1) = t66 / 0.24e2;
cm(2) = t137 + t188;
cm(3) = t204 + t226 / 0.18e2 + t250 + t301;
cm(4) = -t66 / 0.24e2;
%
t1 = (ko ^ 2);
t3 = (j ^ 2);
t5 = (0.1e1 / t1 / t3);
t7 = 3 * t5 * cm(2);
t10 = (0.1e1 / ko / j);
t19 = 8 / t1 / ko / t3 / j * cm(3);
t21 = 2 * t10 * cm(4);
t23 = 2 * t10 * cm(1);
t25 = 8 * t5 * cm(3);
t27 = 4 * t10 * cm(3);
%
coefm(1) = t7;
coefm(2) = cm(1);
coefm(3) = 3 * t10 * cm(2);
coefm(4) = cm(2);
coefm(5) = t19;
coefm(6) = t21;
coefm(7) = t19;
coefm(8) = -t19;
coefm(9) = -t23;
coefm(10) = cm(3);
coefm(11) = -cm(3);
coefm(12) = -t7;
coefm(13) = t23;
coefm(14) = t25;
coefm(15) = -t25;
coefm(16) = t27;
coefm(17) = -t27;
coefm(18) = cm(4);
coefm(19) = -t19;
coefm(20) = -t21;
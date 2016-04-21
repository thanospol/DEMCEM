function [coef,coefm] = coefficients_g2_f2(r1,r2,r3,r4,ko,AreaT)
%% coefficients_g2_f2

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
t2 = t1 * bb(3);
t3 = cc(3) ^ 2;
t4 = t2 * t3;
t6 = sqrt(0.3e1);
t7 = bb(1) * dd(2) * t6;
t10 = bb(3) ^ 2;
t11 = t1 * t10;
t12 = t11 * cc(1);
t13 = cc(3) * dd(2);
t14 = t13 * t6;
t17 = t1 * bb(1);
t18 = cc(1) ^ 2;
t19 = t17 * t18;
t21 = bb(2) * dd(3) * t6;
t24 = bb(1) ^ 2;
t25 = t1 * t24;
t26 = t25 * cc(2);
t27 = cc(1) * dd(3);
t28 = t27 * t6;
t31 = t11 * cc(2);
t32 = cc(3) * dd(1);
t33 = t32 * t6;
t36 = bb(2) ^ 2;
t37 = t1 * t36;
t38 = t37 * cc(1);
t39 = cc(2) * dd(3);
t40 = t39 * t6;
t43 = t1 * bb(2);
t44 = cc(2) ^ 2;
t45 = t43 * t44;
t47 = bb(1) * dd(3) * t6;
t51 = bb(2) * dd(1) * t6;
t54 = t37 * cc(3);
t55 = cc(2) * dd(1);
t56 = t55 * t6;
t59 = t25 * cc(3);
t60 = cc(1) * dd(2);
t61 = t60 * t6;
t65 = bb(3) * dd(2) * t6;
t69 = bb(3) * dd(1) * t6;
t72 = t31 * t47;
t74 = -t4 * t7 / 0.18e2 + t12 * t14 / 0.18e2 - t19 * t21 / 0.18e2 + t26 * t28 / 0.18e2 - t31 * t33 / 0.18e2 - t38 * t40 / 0.18e2 + t45 * t47 / 0.18e2 + t4 * t51 / 0.18e2 + t54 * t56 / 0.18e2 - t59 * t61 / 0.18e2 + t19 * t65 / 0.18e2 - t45 * t69 / 0.18e2 - t72 / 0.72e2;
t75 = t11 * cc(3);
t76 = t75 * t51;
t78 = t75 * t7;
t80 = t37 * cc(2);
t81 = t80 * t69;
t83 = t12 * t21;
t85 = t54 * t7;
t87 = t38 * t65;
t89 = t59 * t51;
t91 = t25 * cc(1);
t92 = t91 * t65;
t94 = t80 * t47;
t96 = t91 * t21;
t98 = t26 * t69;
t100 = t6 * t1;
t101 = t100 * t44;
t102 = cc(3) * bb(2);
t103 = t102 * dd(1);
t106 = t100 * t3;
t107 = cc(2) * bb(1);
t108 = t107 * dd(3);
t112 = cc(1) * bb(3) * dd(2);
t115 = -t76 / 0.72e2 + t78 / 0.72e2 + t81 / 0.72e2 + t83 / 0.72e2 + t85 / 0.72e2 - t87 / 0.72e2 - t89 / 0.72e2 - t92 / 0.72e2 - t94 / 0.72e2 + t96 / 0.72e2 + t98 / 0.72e2 - t101 * t103 / 0.18e2 - t106 * t108 / 0.18e2 - t106 * t112 / 0.18e2;
t118 = t100 * t18;
t121 = cc(3) * bb(1);
t122 = t121 * dd(2);
t125 = bb(3) * cc(2) * dd(1);
t128 = cc(1) * bb(2) * dd(3);
t134 = t1 * cc(2) * cc(1);
t135 = bb(2) * bb(3);
t136 = dd(2) * t6;
t139 = t2 * cc(3);
t140 = dd(3) * t6;
t143 = t43 * cc(2);
t147 = t1 * cc(3) * cc(1);
t150 = -t101 * t112 - t118 * t103 - t118 * t108 + t101 * t122 + t118 * t125 + t101 * t128 + t118 * t122 + t106 * t125 + t106 * t128 + t134 * t135 * t136 + t139 * t107 * t140 - t143 * t121 * t136 - t147 * t135 * t140;
t151 = bb(3) * bb(1);
t152 = dd(1) * t6;
t156 = t17 * cc(1);
t161 = t1 * t24 * bb(1);
t162 = t161 * t40;
t165 = t1 * t36 * bb(2);
t166 = t165 * t28;
t169 = t1 * t10 * bb(3);
t170 = t169 * t56;
t172 = t169 * t61;
t174 = t161 * t14;
t176 = t165 * t33;
t178 = t3 * cc(3);
t187 = t18 * cc(1);
t192 = t44 * cc(2);
t205 = -t134 * t151 * t152 / 0.18e2 + t156 * t102 * t152 / 0.18e2 - t162 / 0.72e2 + t166 / 0.72e2 + t170 / 0.72e2 - t172 / 0.72e2 + t174 / 0.72e2 - t176 / 0.72e2 - t100 * t178 * bb(2) * dd(1) / 0.18e2 + t100 * t178 * bb(1) * dd(2) / 0.18e2 + t100 * t187 * bb(2) * dd(3) / 0.18e2 + t100 * t192 * bb(3) * dd(1) / 0.18e2 - t100 * t187 * bb(3) * dd(2) / 0.18e2 - t100 * t192 * bb(1) * dd(3) / 0.18e2;
t208 = t162 - t96 + t87 + t172 + t72 - t98 + t176 - t170 - t78 - t166 - t83 + t89 + t94 - t174 - t81 + t76 + t92 - t85;
t245 = cc(3) * cc(2) * dd(1);
t251 = cc(3) * cc(1) * dd(2);
t254 = t169 * t60 / 0.12e2 - t165 * t27 / 0.12e2 + t161 * t39 / 0.12e2 - t169 * t55 / 0.12e2 + t165 * t32 / 0.12e2 - t161 * t13 / 0.12e2 - t139 * t108 / 0.6e1 + t143 * t122 / 0.6e1 + t147 * t135 * dd(3) / 0.6e1 - t134 * t135 * dd(2) / 0.6e1 - t156 * t103 / 0.6e1 + t134 * t151 * dd(1) / 0.6e1 - t37 * t125 / 0.12e2 + t43 * t44 * bb(3) * dd(1) / 0.6e1 - t37 * t122 / 0.12e2 - t37 * t245 / 0.6e1 + t25 * t103 / 0.12e2 + t25 * t251 / 0.6e1;
t276 = cc(1) * cc(2) * dd(3);
t303 = t25 * t112 / 0.12e2 - t17 * t18 * bb(3) * dd(2) / 0.6e1 - t11 * t122 / 0.12e2 - t11 * t128 / 0.12e2 + t37 * t112 / 0.12e2 + t37 * t108 / 0.12e2 + t2 * t3 * bb(1) * dd(2) / 0.6e1 - t11 * t251 / 0.6e1 + t37 * t276 / 0.6e1 - t43 * t44 * bb(1) * dd(3) / 0.6e1 - t25 * t128 / 0.12e2 - t25 * t125 / 0.12e2 + t11 * t108 / 0.12e2 + t11 * t103 / 0.12e2 + t17 * t18 * bb(2) * dd(3) / 0.6e1 - t25 * t276 / 0.6e1 + t11 * t245 / 0.6e1 - t2 * t3 * bb(2) * dd(1) / 0.6e1;
%
c(1) = t74 + t115 + t150 / 0.18e2 + t205;
c(2) = t208 / 0.24e2;
c(3) = t254 + t303;
c(4) = -t208 / 0.24e2;
%
t3 = (0.1e1 / ko / j);
t7 = 2 * t3 * c(2);
t9 = 2 * t3 * c(4);
t10 = (ko ^ 2);
t13 = (j ^ 2);
t18 = 8 / t10 / ko / t13 / j * c(1);
t21 = 1 / t10 / t13;
t23 = 3 * t21 * c(3);
t25 = 8 * t21 * c(1);
t27 = 4 * t3 * c(1);
%
coef(1) = 3 * t3 * c(3);
coef(2) = c(3);
coef(3) = -t7;
coef(4) = t7;
coef(5) = -t9;
coef(6) = -t18;
coef(7) = t18;
coef(8) = t18;
coef(9) = -t18;
coef(10) = c(4);
coef(11) = t23;
coef(12) = -c(1);
coef(13) = c(1);
coef(14) = c(2);
coef(15) = -t23;
coef(16) = t9;
coef(17) = t25;
coef(18) = -t25;
coef(19) = t27;
coef(20) = -t27;
%
t1 = 0.1e1 / AreaT;
t2 = bb(3) ^ 2;
t4 = t1 * t2 * bb(3);
t5 = cc(1) * dd(2);
t8 = bb(2) ^ 2;
t10 = t1 * t8 * bb(2);
t11 = cc(1) * dd(3);
t14 = bb(1) ^ 2;
t16 = t1 * t14 * bb(1);
t17 = cc(2) * dd(3);
t20 = cc(2) * dd(1);
t23 = cc(3) * dd(1);
t26 = cc(3) * dd(2);
t29 = t1 * bb(2);
t30 = t29 * cc(2);
t31 = cc(3) * bb(1);
t32 = t31 * dd(2);
t35 = t1 * bb(3);
t36 = t35 * cc(3);
t37 = cc(2) * bb(1);
t38 = t37 * dd(3);
t42 = t1 * cc(3) * cc(1);
t43 = bb(3) * bb(2);
t48 = t1 * cc(2) * cc(1);
t52 = t1 * bb(1);
t53 = t52 * cc(1);
t54 = cc(3) * bb(2);
t55 = t54 * dd(1);
t58 = bb(1) * bb(3);
t62 = cc(2) ^ 2;
t67 = cc(3) ^ 2;
t72 = t1 * t2;
t74 = cc(1) * cc(3) * dd(2);
t77 = t1 * t8;
t79 = cc(1) * cc(2) * dd(3);
t82 = cc(1) ^ 2;
t87 = t1 * t14;
t90 = -t4 * t5 / 0.12e2 + t10 * t11 / 0.12e2 - t16 * t17 / 0.12e2 + t4 * t20 / 0.12e2 - t10 * t23 / 0.12e2 + t16 * t26 / 0.12e2 - t30 * t32 / 0.6e1 + t36 * t38 / 0.6e1 - t42 * t43 * dd(3) / 0.6e1 + t48 * t43 * dd(2) / 0.6e1 + t53 * t55 / 0.6e1 - t48 * t58 * dd(1) / 0.6e1 + t29 * t62 * bb(1) * dd(3) / 0.6e1 - t35 * t67 * bb(1) * dd(2) / 0.6e1 + t72 * t74 / 0.6e1 - t77 * t79 / 0.6e1 - t52 * t82 * bb(2) * dd(3) / 0.6e1 + t87 * t79 / 0.6e1;
t92 = cc(2) * cc(3) * dd(1);
t112 = cc(1) * bb(2) * dd(3);
t120 = cc(1) * bb(3) * dd(2);
t124 = cc(2) * bb(3) * dd(1);
t141 = -t72 * t92 / 0.6e1 + t35 * t67 * bb(2) * dd(1) / 0.6e1 - t29 * t62 * bb(3) * dd(1) / 0.6e1 + t77 * t92 / 0.6e1 - t87 * t74 / 0.6e1 + t52 * t82 * bb(3) * dd(2) / 0.6e1 + t72 * t112 / 0.12e2 + t72 * t32 / 0.12e2 - t77 * t38 / 0.12e2 - t77 * t120 / 0.12e2 + t87 * t124 / 0.12e2 + t87 * t112 / 0.12e2 - t72 * t55 / 0.12e2 - t72 * t38 / 0.12e2 - t87 * t55 / 0.12e2 + t77 * t32 / 0.12e2 + t77 * t124 / 0.12e2 - t87 * t120 / 0.12e2;
t143 = sqrt(0.3e1);
t144 = t26 * t143;
t145 = t16 * t144;
t146 = t87 * cc(1);
t148 = bb(2) * dd(3) * t143;
t149 = t146 * t148;
t150 = t77 * cc(2);
t152 = bb(3) * dd(1) * t143;
t153 = t150 * t152;
t154 = t77 * cc(1);
t156 = bb(3) * dd(2) * t143;
t157 = t154 * t156;
t158 = t72 * cc(1);
t159 = t158 * t148;
t160 = t17 * t143;
t161 = t16 * t160;
t162 = t72 * cc(3);
t164 = bb(2) * dd(1) * t143;
t165 = t162 * t164;
t166 = t20 * t143;
t167 = t4 * t166;
t168 = t72 * cc(2);
t170 = bb(1) * dd(3) * t143;
t171 = t168 * t170;
t172 = t23 * t143;
t173 = t10 * t172;
t174 = t11 * t143;
t175 = t10 * t174;
t176 = t87 * cc(2);
t177 = t176 * t152;
t179 = bb(1) * dd(2) * t143;
t180 = t162 * t179;
t181 = t150 * t170;
t182 = t5 * t143;
t183 = t4 * t182;
t184 = t77 * cc(3);
t185 = t184 * t179;
t186 = t87 * cc(3);
t187 = t186 * t164;
t188 = t146 * t156;
t189 = t145 + t149 + t153 - t157 + t159 - t161 - t165 + t167 - t171 - t173 + t175 + t177 + t180 - t181 - t183 + t185 - t187 - t188;
t190 = dd(2) * t143;
t194 = dd(3) * t143;
t204 = dd(1) * t143;
t218 = -t30 * t31 * t190 / 0.18e2 + t36 * t37 * t194 / 0.18e2 - t42 * t43 * t194 / 0.18e2 + t48 * t43 * t190 / 0.18e2 + t53 * t54 * t204 / 0.18e2 - t48 * t58 * t204 / 0.18e2 + t149 / 0.72e2 + t153 / 0.72e2 - t157 / 0.72e2 + t159 / 0.72e2 - t165 / 0.72e2 - t171 / 0.72e2 + t177 / 0.72e2;
t224 = t29 * t62;
t227 = t35 * t67;
t234 = t52 * t82;
t245 = t180 / 0.72e2 - t181 / 0.72e2 + t185 / 0.72e2 - t187 / 0.72e2 - t188 / 0.72e2 + t224 * t170 / 0.18e2 - t227 * t179 / 0.18e2 + t158 * t144 / 0.18e2 - t154 * t160 / 0.18e2 + t234 * t156 / 0.18e2 - t168 * t172 / 0.18e2 - t234 * t148 / 0.18e2 + t176 * t174 / 0.18e2 + t227 * t164 / 0.18e2;
t250 = t143 * t1;
t251 = t250 * t82;
t253 = t250 * t67;
t255 = t250 * t62;
t264 = -t186 * t182 - t224 * t152 + t184 * t166 + t251 * t32 - t253 * t38 - t255 * t55 + t255 * t32 + t253 * t112 + t251 * t124 - t255 * t120 - t253 * t120 + t255 * t112 + t253 * t124;
t275 = t82 * cc(1);
t280 = t62 * cc(2);
t289 = t67 * cc(3);
t302 = -t251 * t38 / 0.18e2 - t251 * t55 / 0.18e2 + t145 / 0.72e2 - t161 / 0.72e2 + t167 / 0.72e2 - t173 / 0.72e2 + t175 / 0.72e2 - t183 / 0.72e2 - t250 * t275 * bb(3) * dd(2) / 0.18e2 - t250 * t280 * bb(1) * dd(3) / 0.18e2 + t250 * t280 * bb(3) * dd(1) / 0.18e2 + t250 * t289 * bb(1) * dd(2) / 0.18e2 + t250 * t275 * bb(2) * dd(3) / 0.18e2 - t250 * t289 * bb(2) * dd(1) / 0.18e2;
%
cm(1) = t90 + t141;
cm(2) = t189 / 0.24e2;
cm(3) = t218 + t245 + t264 / 0.18e2 + t302;
cm(4) = -t189 / 0.24e2;
%
t3 = (0.1e1 / ko / j);
t5 = 2 * t3 * cm(4);
t8 = (ko ^ 2);
t11 = (j ^ 2);
t16 = 8 / t8 / ko / t11 / j * cm(3);
t18 = 2 * t3 * cm(2);
t21 = 1 / t8 / t11;
t23 = 3 * t21 * cm(1);
t25 = 8 * t21 * cm(3);
t27 = 4 * t3 * cm(3);
%
coefm(1) = t5;
coefm(2) = 3 * t3 * cm(1);
coefm(3) = cm(1);
coefm(4) = -t16;
coefm(5) = t16;
coefm(6) = -t5;
coefm(7) = -t18;
coefm(8) = cm(4);
coefm(9) = cm(3);
coefm(10) = -cm(3);
coefm(11) = -t23;
coefm(12) = t18;
coefm(13) = t25;
coefm(14) = -t25;
coefm(15) = t27;
coefm(16) = -t27;
coefm(17) = -t16;
coefm(18) = t16;
coefm(19) = cm(2);
coefm(20) = t23;
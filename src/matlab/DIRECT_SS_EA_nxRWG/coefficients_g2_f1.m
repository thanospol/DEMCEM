function [coef,coefm] = coefficients_g2_f1(r1,r2,r3,r4,ko,AreaT)
%% coefficients_g2_f1

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
t2 = bb(3) ^ 2;
t3 = t1 * t2;
t5 = bb(2) * dd(3);
t6 = sqrt(0.3e1);
t7 = t5 * t6;
t8 = t3 * cc(1) * t7;
t9 = t3 * cc(3);
t10 = bb(2) * dd(1);
t11 = t10 * t6;
t12 = t9 * t11;
t13 = bb(2) ^ 2;
t14 = t1 * t13;
t17 = bb(3) * dd(2) * t6;
t18 = t14 * cc(1) * t17;
t19 = t14 * cc(2);
t21 = bb(1) * dd(3) * t6;
t22 = t19 * t21;
t24 = t1 * t13 * bb(2);
t25 = cc(3) * dd(1);
t27 = t24 * t25 * t6;
t28 = bb(1) ^ 2;
t29 = t1 * t28;
t30 = t29 * cc(1);
t31 = t30 * t17;
t32 = bb(1) * dd(2);
t33 = t32 * t6;
t34 = t9 * t33;
t36 = t1 * t28 * bb(1);
t37 = cc(2) * dd(3);
t39 = t36 * t37 * t6;
t42 = bb(3) * dd(1) * t6;
t43 = t29 * cc(2) * t42;
t45 = t1 * t2 * bb(3);
t46 = cc(2) * dd(1);
t48 = t45 * t46 * t6;
t49 = t30 * t7;
t50 = cc(3) * dd(2);
t52 = t36 * t50 * t6;
t53 = t19 * t42;
t54 = cc(1) * dd(3);
t56 = t24 * t54 * t6;
t58 = t3 * cc(2) * t21;
t59 = cc(1) * dd(2);
t61 = t45 * t59 * t6;
t63 = t29 * cc(3) * t11;
t65 = t14 * cc(3) * t33;
t66 = t8 - t12 - t18 - t22 - t27 - t31 + t34 - t39 + t43 + t48 + t49 + t52 + t53 + t56 - t58 - t61 - t63 + t65;
t79 = t1 * cc(1);
t80 = t79 * cc(3);
t85 = t8 / 0.72e2 - t12 / 0.72e2 - t18 / 0.72e2 - t22 / 0.72e2 - t31 / 0.72e2 + t34 / 0.72e2 + t43 / 0.72e2 + t49 / 0.72e2 + t53 / 0.72e2 - t58 / 0.72e2 - t63 / 0.72e2 + t65 / 0.72e2 - t80 * t28 * dd(2) * t6 / 0.18e2;
t86 = cc(1) ^ 2;
t87 = t1 * t86;
t93 = t1 * cc(2);
t94 = t93 * cc(1);
t98 = cc(2) ^ 2;
t99 = t1 * t98;
t102 = cc(3) ^ 2;
t103 = t1 * t102;
t113 = t1 * cc(3);
t114 = t113 * cc(2);
t118 = t1 * bb(2);
t124 = t6 * t1;
t125 = t124 * t86;
t126 = cc(2) * bb(1);
t127 = t126 * dd(3);
t129 = t124 * t102;
t132 = bb(2) * cc(1) * dd(3);
t134 = t87 * bb(3) * t33 + t80 * t2 * dd(2) * t6 - t94 * t13 * dd(3) * t6 + t99 * bb(1) * t7 - t103 * bb(1) * t17 - t87 * bb(2) * t21 + t94 * t28 * dd(3) * t6 + t103 * bb(2) * t42 - t114 * t2 * dd(1) * t6 - t118 * t98 * t42 + t114 * t13 * dd(1) * t6 - t125 * t127 - t129 * t127 + t129 * t132;
t136 = t124 * t98;
t138 = bb(3) * cc(1) * dd(2);
t140 = cc(3) * bb(1);
t141 = t140 * dd(2);
t145 = cc(2) * bb(3) * dd(1);
t150 = cc(3) * bb(2);
t151 = t150 * dd(1);
t154 = t118 * cc(2);
t155 = dd(2) * t6;
t158 = bb(2) * bb(3);
t162 = t1 * bb(1) * cc(1);
t163 = dd(1) * t6;
t167 = t1 * bb(3) * cc(3);
t168 = dd(3) * t6;
t171 = -t136 * t138 + t125 * t141 + t136 * t132 + t129 * t145 - t129 * t138 + t125 * t145 + t136 * t141 - t125 * t151 - t136 * t151 - t154 * t140 * t155 + t94 * t158 * t155 + t162 * t150 * t163 + t167 * t126 * t168;
t172 = bb(3) * bb(1);
t185 = t86 * cc(1);
t194 = t98 * cc(2);
t199 = t102 * cc(3);
t212 = -t94 * t172 * t163 / 0.18e2 - t80 * t158 * t168 / 0.18e2 - t27 / 0.72e2 - t39 / 0.72e2 + t48 / 0.72e2 + t52 / 0.72e2 + t56 / 0.72e2 - t61 / 0.72e2 - t124 * t185 * bb(3) * dd(2) / 0.18e2 + t124 * t185 * bb(2) * dd(3) / 0.18e2 + t124 * t194 * bb(3) * dd(1) / 0.18e2 + t124 * t199 * bb(1) * dd(2) / 0.18e2 - t124 * t199 * bb(2) * dd(1) / 0.18e2 - t124 * t194 * bb(1) * dd(3) / 0.18e2;
t250 = t5 * bb(1);
t257 = -t45 * t59 / 0.12e2 + t24 * t54 / 0.12e2 - t36 * t37 / 0.12e2 + t45 * t46 / 0.12e2 - t24 * t25 / 0.12e2 + t36 * t50 / 0.12e2 + t167 * t127 / 0.6e1 - t154 * t141 / 0.6e1 - t80 * t158 * dd(3) / 0.6e1 + t94 * t158 * dd(2) / 0.6e1 + t162 * t151 / 0.6e1 - t94 * t172 * dd(1) / 0.6e1 - t29 * t151 / 0.12e2 - t29 * t138 / 0.12e2 + t14 * t145 / 0.12e2 + t14 * t141 / 0.12e2 - t87 * t250 / 0.6e1 + t79 * cc(2) * t28 * dd(3) / 0.6e1;
t273 = t32 * bb(3);
t308 = -t113 * cc(2) * t2 * dd(1) / 0.6e1 + t103 * t10 * bb(3) / 0.6e1 + t3 * t141 / 0.12e2 + t3 * t132 / 0.12e2 - t14 * t138 / 0.12e2 - t14 * t127 / 0.12e2 - t103 * t273 / 0.6e1 + t113 * cc(1) * t2 * dd(2) / 0.6e1 - t93 * cc(1) * t13 * dd(3) / 0.6e1 + t99 * t250 / 0.6e1 + t29 * t132 / 0.12e2 + t29 * t145 / 0.12e2 - t3 * t127 / 0.12e2 - t3 * t151 / 0.12e2 - t118 * t98 * bb(3) * dd(1) / 0.6e1 + t93 * cc(3) * t13 * dd(1) / 0.6e1 - t79 * cc(3) * t28 * dd(2) / 0.6e1 + t87 * t273 / 0.6e1;
t309 = t257 + t308;
%
c(1) = t66 / 0.24e2;
c(2) = t85 + t134 / 0.18e2 + t171 / 0.18e2 + t212;
c(3) = t66 / 0.24e2;
c(4) = t309;
c(5) = -t309;
c(6) = -t66 / 0.12e2;
%
t3 = (0.1e1 / ko / j);
t6 = (ko ^ 2);
t9 = (j ^ 2);
t14 = 8 / t6 / ko / t9 / j * c(2);
t16 = 2 * t3 * c(1);
t18 = 2 * t3 * c(6);
t21 = 1 / t6 / t9;
t23 = 3 * t21 * c(4);
t25 = 2 * t3 * c(3);
t29 = 3 * t21 * c(5);
t31 = 4 * t3 * c(2);
t33 = 8 * t21 * c(2);
%
coef(1) = c(5);
coef(2) = 3 * t3 * c(5);
coef(3) = t14;
coef(4) = -t16;
coef(5) = -t18;
coef(6) = -t23;
coef(7) = t25;
coef(8) = -t14;
coef(9) = t14;
coef(10) = -t25;
coef(11) = -t14;
coef(12) = c(1);
coef(13) = c(6);
coef(14) = c(4);
coef(15) = 3 * t3 * c(4);
coef(16) = t29;
coef(17) = -c(2);
coef(18) = c(2);
coef(19) = c(3);
coef(20) = t16;
coef(21) = -t31;
coef(22) = -t33;
coef(23) = t18;
coef(24) = -t29;
coef(25) = t31;
coef(26) = t33;
coef(27) = t23;
%
t1 = 0.1e1 / AreaT;
t2 = t1 * bb(2);
t3 = t2 * cc(2);
t4 = cc(3) * bb(1);
t5 = t4 * dd(2);
t9 = t1 * cc(3) * cc(1);
t10 = bb(3) * bb(2);
t15 = t1 * cc(2) * cc(1);
t19 = bb(3) ^ 2;
t20 = t1 * t19;
t22 = cc(1) * cc(3) * dd(2);
t25 = bb(2) ^ 2;
t26 = t1 * t25;
t28 = cc(1) * cc(2) * dd(3);
t31 = cc(2) ^ 2;
t36 = t1 * bb(1);
t37 = cc(1) ^ 2;
t42 = bb(1) ^ 2;
t43 = t1 * t42;
t47 = cc(2) * cc(3) * dd(1);
t50 = t1 * bb(3);
t51 = cc(3) ^ 2;
t69 = cc(1) * bb(2) * dd(3);
t73 = cc(1) * bb(3) * dd(2);
t76 = cc(2) * bb(1);
t77 = t76 * dd(3);
t82 = -t3 * t5 / 0.6e1 - t9 * t10 * dd(3) / 0.6e1 + t15 * t10 * dd(2) / 0.6e1 + t20 * t22 / 0.6e1 - t26 * t28 / 0.6e1 + t2 * t31 * bb(1) * dd(3) / 0.6e1 - t36 * t37 * bb(2) * dd(3) / 0.6e1 + t43 * t28 / 0.6e1 - t20 * t47 / 0.6e1 + t50 * t51 * bb(2) * dd(1) / 0.6e1 - t2 * t31 * bb(3) * dd(1) / 0.6e1 + t26 * t47 / 0.6e1 - t43 * t22 / 0.6e1 + t36 * t37 * bb(3) * dd(2) / 0.6e1 + t20 * t69 / 0.12e2 - t26 * t73 / 0.12e2 - t26 * t77 / 0.12e2 + t43 * t69 / 0.12e2;
t84 = cc(2) * bb(3) * dd(1);
t89 = cc(3) * bb(2);
t90 = t89 * dd(1);
t108 = t1 * t19 * bb(3);
t109 = cc(1) * dd(2);
t113 = t1 * t25 * bb(2);
t114 = cc(1) * dd(3);
t118 = t1 * t42 * bb(1);
t119 = cc(2) * dd(3);
t122 = cc(2) * dd(1);
t125 = cc(3) * dd(1);
t128 = cc(3) * dd(2);
t131 = t50 * cc(3);
t134 = bb(1) * bb(3);
t138 = t36 * cc(1);
t141 = t43 * t84 / 0.12e2 - t20 * t77 / 0.12e2 - t20 * t90 / 0.12e2 + t26 * t84 / 0.12e2 + t26 * t5 / 0.12e2 - t43 * t90 / 0.12e2 - t43 * t73 / 0.12e2 + t20 * t5 / 0.12e2 - t50 * t51 * bb(1) * dd(2) / 0.6e1 - t108 * t109 / 0.12e2 + t113 * t114 / 0.12e2 - t118 * t119 / 0.12e2 + t108 * t122 / 0.12e2 - t113 * t125 / 0.12e2 + t118 * t128 / 0.12e2 + t131 * t77 / 0.6e1 - t15 * t134 * dd(1) / 0.6e1 + t138 * t90 / 0.6e1;
t142 = t82 + t141;
t143 = t26 * cc(2);
t145 = sqrt(0.3e1);
t146 = bb(1) * dd(3) * t145;
t147 = t143 * t146;
t148 = t20 * cc(3);
t150 = bb(2) * dd(1) * t145;
t151 = t148 * t150;
t153 = bb(3) * dd(1) * t145;
t154 = t143 * t153;
t156 = bb(1) * dd(2) * t145;
t157 = t148 * t156;
t158 = t109 * t145;
t159 = t108 * t158;
t160 = t20 * cc(1);
t162 = bb(2) * dd(3) * t145;
t163 = t160 * t162;
t164 = t114 * t145;
t165 = t113 * t164;
t166 = t26 * cc(1);
t168 = bb(3) * dd(2) * t145;
t169 = t166 * t168;
t170 = t43 * cc(1);
t171 = t170 * t162;
t172 = t26 * cc(3);
t173 = t172 * t156;
t174 = t43 * cc(2);
t175 = t174 * t153;
t176 = t43 * cc(3);
t177 = t176 * t150;
t178 = t20 * cc(2);
t179 = t178 * t146;
t180 = t170 * t168;
t181 = t128 * t145;
t182 = t118 * t181;
t183 = t122 * t145;
t184 = t108 * t183;
t185 = t119 * t145;
t186 = t118 * t185;
t187 = t125 * t145;
t188 = t113 * t187;
t189 = -t147 - t151 + t154 + t157 - t159 + t163 + t165 - t169 + t171 + t173 + t175 - t177 - t179 - t180 + t182 + t184 - t186 - t188;
t190 = dd(2) * t145;
t197 = dd(1) * t145;
t204 = dd(3) * t145;
t217 = t145 * t1;
t218 = t31 * cc(2);
t223 = -t3 * t4 * t190 / 0.18e2 + t15 * t10 * t190 / 0.18e2 - t15 * t134 * t197 / 0.18e2 + t138 * t89 * t197 / 0.18e2 + t131 * t76 * t204 / 0.18e2 - t9 * t10 * t204 / 0.18e2 - t159 / 0.72e2 + t165 / 0.72e2 + t182 / 0.72e2 + t184 / 0.72e2 - t186 / 0.72e2 - t188 / 0.72e2 + t217 * t218 * bb(3) * dd(1) / 0.18e2;
t224 = t37 * cc(1);
t229 = t51 * cc(3);
t255 = -t217 * t224 * bb(3) * dd(2) / 0.18e2 + t217 * t229 * bb(1) * dd(2) / 0.18e2 - t217 * t218 * bb(1) * dd(3) / 0.18e2 + t217 * t224 * bb(2) * dd(3) / 0.18e2 - t217 * t229 * bb(2) * dd(1) / 0.18e2 - t147 / 0.72e2 - t151 / 0.72e2 + t154 / 0.72e2 + t157 / 0.72e2 + t163 / 0.72e2 - t169 / 0.72e2 + t171 / 0.72e2 + t173 / 0.72e2 + t175 / 0.72e2;
t264 = t50 * t51;
t267 = t2 * t31;
t272 = t36 * t37;
t283 = -t177 / 0.72e2 - t179 / 0.72e2 - t180 / 0.72e2 + t174 * t164 / 0.18e2 - t178 * t187 / 0.18e2 - t264 * t156 / 0.18e2 + t267 * t146 / 0.18e2 + t172 * t183 / 0.18e2 + t272 * t168 / 0.18e2 - t166 * t185 / 0.18e2 - t176 * t158 / 0.18e2 - t272 * t162 / 0.18e2 - t267 * t153 / 0.18e2;
t286 = t217 * t37;
t288 = t217 * t31;
t290 = t217 * t51;
t301 = t160 * t181 + t264 * t150 + t286 * t84 - t288 * t90 + t290 * t69 - t290 * t73 - t286 * t77 - t290 * t77 + t288 * t5 - t288 * t73 + t290 * t84 + t286 * t5 - t286 * t90 + t288 * t69;
cm(1) = t142;
cm(2) = t142;
cm(3) = t189 / 0.12e2;
cm(4) = t189 / 0.24e2;
cm(5) = t189 / 0.24e2;
cm(6) = t223 + t255 + t283 + t301 / 0.18e2;
t3 = (0.1e1 / ko / j);
t6 = (ko ^ 2);
t8 = (j ^ 2);
t10 = (0.1e1 / t6 / t8);
t12 = 3 * t10 * cm(2);
t21 = 8 / t6 / ko / t8 / j * cm(6);
t23 = 2 * t3 * cm(4);
t25 = 3 * t10 * cm(1);
t27 = 2 * t3 * cm(5);
t29 = 2 * t3 * cm(3);
t31 = 4 * t3 * cm(6);
t33 = 8 * t10 * cm(6);
%
coefm(1) = 3 * t3 * cm(2);
coefm(2) = cm(2);
coefm(3) = cm(3);
coefm(4) = cm(5);
coefm(5) = cm(1);
coefm(6) = t12;
coefm(7) = 3 * t3 * cm(1);
coefm(8) = t21;
coefm(9) = -t21;
coefm(10) = -t23;
coefm(11) = t23;
coefm(12) = -t21;
coefm(13) = -t25;
coefm(14) = t21;
coefm(15) = -t27;
coefm(16) = -t29;
coefm(17) = cm(4);
coefm(18) = cm(6);
coefm(19) = -cm(6);
coefm(20) = t25;
coefm(21) = t29;
coefm(22) = -t31;
coefm(23) = t31;
coefm(24) = t27;
coefm(25) = -t33;
coefm(26) = -t12;
coefm(27) = t33;
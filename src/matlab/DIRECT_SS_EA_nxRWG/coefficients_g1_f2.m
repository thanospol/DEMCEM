function [coef,coefm] = coefficients_g1_f2(r1,r2,r3,r4,ko,AreaT)
%% coefficients_g1_f2

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
t1 = sqrt(0.3e1);
t2 = 0.1e1 / AreaT;
t3 = t1 * t2;
t4 = bb(3) ^ 2;
t5 = t4 * bb(3);
t8 = t3 * t5 * cc(1) * dd(2);
t9 = bb(2) ^ 2;
t10 = t9 * bb(2);
t13 = t3 * t10 * cc(1) * dd(3);
t14 = bb(1) ^ 2;
t15 = t14 * bb(1);
t18 = t3 * t15 * cc(2) * dd(3);
t21 = t3 * t5 * cc(2) * dd(1);
t24 = t3 * t10 * cc(3) * dd(1);
t27 = t3 * t15 * cc(3) * dd(2);
t28 = t3 * t14;
t29 = cc(3) * bb(2);
t30 = t29 * dd(1);
t31 = t28 * t30;
t33 = cc(1) * bb(3) * dd(2);
t34 = t28 * t33;
t35 = t3 * t9;
t36 = t35 * t33;
t37 = cc(2) * bb(1);
t38 = t37 * dd(3);
t39 = t35 * t38;
t40 = t3 * t4;
t41 = cc(3) * bb(1);
t42 = t41 * dd(2);
t43 = t40 * t42;
t45 = cc(1) * bb(2) * dd(3);
t46 = t40 * t45;
t47 = t40 * t38;
t48 = t40 * t30;
t49 = t28 * t45;
t51 = cc(2) * bb(3) * dd(1);
t52 = t28 * t51;
t53 = t35 * t51;
t54 = t35 * t42;
t55 = -t8 + t13 - t18 + t21 - t24 + t27 - t31 - t34 - t36 - t39 + t43 + t46 - t47 - t48 + t49 + t52 + t53 + t54;
t56 = t2 * bb(3);
t57 = t56 * cc(3);
t58 = dd(3) * t1;
t62 = t2 * bb(2);
t63 = t62 * cc(2);
t64 = dd(2) * t1;
t69 = t2 * cc(3) * cc(1);
t70 = bb(2) * bb(3);
t75 = t2 * cc(2) * cc(1);
t79 = t2 * bb(1);
t80 = t79 * cc(1);
t81 = dd(1) * t1;
t85 = bb(3) * bb(1);
t95 = cc(3) ^ 2;
t96 = t95 * cc(3);
t101 = t57 * t37 * t58 / 0.18e2 - t63 * t41 * t64 / 0.18e2 - t69 * t70 * t58 / 0.18e2 + t75 * t70 * t64 / 0.18e2 + t80 * t29 * t81 / 0.18e2 - t75 * t85 * t81 / 0.18e2 - t8 / 0.72e2 + t13 / 0.72e2 - t18 / 0.72e2 + t21 / 0.72e2 - t24 / 0.72e2 + t27 / 0.72e2 + t3 * t96 * bb(1) * dd(2) / 0.18e2;
t102 = cc(2) ^ 2;
t103 = t102 * cc(2);
t108 = cc(1) ^ 2;
t109 = t108 * cc(1);
t135 = -t3 * t103 * bb(1) * dd(3) / 0.18e2 + t3 * t109 * bb(2) * dd(3) / 0.18e2 - t3 * t96 * bb(2) * dd(1) / 0.18e2 + t3 * t103 * bb(3) * dd(1) / 0.18e2 - t3 * t109 * bb(3) * dd(2) / 0.18e2 - t31 / 0.72e2 - t34 / 0.72e2 - t36 / 0.72e2 - t39 / 0.72e2 + t43 / 0.72e2 + t46 / 0.72e2 - t47 / 0.72e2 - t48 / 0.72e2 + t49 / 0.72e2;
t140 = t2 * t9;
t142 = cc(2) * dd(3);
t146 = t2 * t4;
t148 = cc(3) * dd(2);
t152 = t56 * t95;
t157 = t62 * t102;
t162 = t79 * t108;
t167 = t2 * t14;
t169 = cc(1) * dd(3);
t174 = cc(3) * dd(1);
t187 = cc(2) * dd(1);
t191 = t52 / 0.72e2 + t53 / 0.72e2 + t54 / 0.72e2 - t140 * cc(1) * t142 * t1 / 0.18e2 + t146 * cc(1) * t148 * t1 / 0.18e2 - t152 * bb(1) * dd(2) * t1 / 0.18e2 + t157 * bb(1) * dd(3) * t1 / 0.18e2 - t162 * bb(2) * dd(3) * t1 / 0.18e2 + t167 * cc(2) * t169 * t1 / 0.18e2 - t146 * cc(2) * t174 * t1 / 0.18e2 + t152 * bb(2) * dd(1) * t1 / 0.18e2 - t157 * bb(3) * dd(1) * t1 / 0.18e2 + t140 * cc(3) * t187 * t1 / 0.18e2;
t193 = cc(1) * dd(2);
t199 = t3 * t95;
t201 = t3 * t108;
t207 = t3 * t102;
t214 = -t167 * cc(3) * t193 * t1 + t162 * bb(3) * dd(2) * t1 + t199 * t51 - t201 * t38 - t199 * t33 - t199 * t38 + t201 * t42 + t201 * t51 + t207 * t42 - t207 * t30 + t199 * t45 - t201 * t30 + t207 * t45 - t207 * t33;
t217 = t2 * t5;
t219 = t2 * t10;
t221 = t2 * t15;
t238 = t217 * t193 - t219 * t169 + t221 * t142 - t217 * t187 + t219 * t174 - t221 * t148 - t146 * t42 - t146 * t45 + t140 * t33 + t140 * t38 - t167 * t45 - t167 * t51 + t146 * t38 + t146 * t30 - t140 * t51 - t140 * t42 + t167 * t30 + t167 * t33;
t240 = cc(1) * cc(2) * dd(3);
t243 = cc(1) * cc(3) * dd(2);
t256 = cc(2) * cc(3) * dd(1);
t278 = t140 * t240 - t146 * t243 + t56 * t95 * bb(1) * dd(2) - t62 * t102 * bb(1) * dd(3) + t79 * t108 * bb(2) * dd(3) - t167 * t240 + t146 * t256 - t56 * t95 * bb(2) * dd(1) + t62 * t102 * bb(3) * dd(1) - t140 * t256 + t167 * t243 - t79 * t108 * bb(3) * dd(2) - t57 * t38 + t63 * t42 + t69 * t70 * dd(3) - t75 * t70 * dd(2) - t80 * t30 + t75 * t85 * dd(1);
t279 = t238 / 0.12e2 + t278 / 0.6e1;
%
c(1) = t55 / 0.24e2;
c(2) = t101 + t135 + t191 + t214 / 0.18e2;
c(3) = t279;
c(4) = t279;
c(5) = t55 / 0.24e2;
c(6) = t55 / 0.12e2;
%
t3 = (0.1e1 / ko / j);
t6 = (ko ^ 2);
t8 = (j ^ 2);
t10 = (0.1e1 / t6 / t8);
t12 = 3 * t10 * c(3);
t19 = 8 / t6 / ko / t8 / j * c(2);
t21 = 2 * t3 * c(6);
t23 = 2 * t3 * c(1);
t25 = 2 * t3 * c(5);
t27 = 3 * t10 * c(4);
t29 = 8 * t10 * c(2);
t31 = 4 * t3 * c(2);
%
coef(1) = c(4);
coef(2) = 3 * t3 * c(4);
coef(3) = -t12;
coef(4) = t19;
coef(5) = -t21;
coef(6) = -t23;
coef(7) = -t19;
coef(8) = t25;
coef(9) = -t19;
coef(10) = -t25;
coef(11) = t19;
coef(12) = c(2);
coef(13) = -c(2);
coef(14) = -t27;
coef(15) = t29;
coef(16) = -t29;
coef(17) = t21;
coef(18) = t31;
coef(19) = -t31;
coef(20) = t23;
coef(21) = t12;
coef(22) = c(1);
coef(23) = t27;
coef(24) = 3 * t3 * c(3);
coef(25) = c(3);
coef(26) = c(6);
coef(27) = c(5);
%
t1 = sqrt(0.3e1);
t2 = 0.1e1 / AreaT;
t3 = t1 * t2;
t4 = bb(3) ^ 2;
t5 = t3 * t4;
t6 = cc(3) * bb(1);
t7 = t6 * dd(2);
t8 = t5 * t7;
t10 = cc(1) * bb(2) * dd(3);
t11 = t5 * t10;
t12 = bb(2) ^ 2;
t13 = t3 * t12;
t14 = cc(2) * bb(1);
t15 = t14 * dd(3);
t16 = t13 * t15;
t18 = cc(1) * bb(3) * dd(2);
t19 = t13 * t18;
t20 = t5 * t15;
t21 = bb(1) ^ 2;
t22 = t3 * t21;
t24 = cc(2) * bb(3) * dd(1);
t25 = t22 * t24;
t26 = t22 * t10;
t27 = cc(3) * bb(2);
t28 = t27 * dd(1);
t29 = t5 * t28;
t30 = t22 * t18;
t31 = t22 * t28;
t32 = t13 * t7;
t33 = t13 * t24;
t34 = t4 * bb(3);
t37 = t3 * t34 * cc(2) * dd(1);
t38 = t12 * bb(2);
t41 = t3 * t38 * cc(1) * dd(3);
t44 = t3 * t34 * cc(1) * dd(2);
t47 = t3 * t38 * cc(3) * dd(1);
t48 = t21 * bb(1);
t51 = t3 * t48 * cc(2) * dd(3);
t54 = t3 * t48 * cc(3) * dd(2);
t55 = t8 + t11 - t16 - t19 - t20 + t25 + t26 - t29 - t30 - t31 + t32 + t33 + t37 + t41 - t44 - t47 - t51 + t54;
t56 = t2 * bb(3);
t57 = t56 * cc(3);
t59 = t2 * bb(2);
t60 = t59 * cc(2);
t63 = t2 * cc(3) * cc(1);
t64 = bb(3) * bb(2);
t68 = t2 * cc(2) * cc(1);
t71 = t2 * bb(1);
t72 = t71 * cc(1);
t74 = bb(1) * bb(3);
t77 = t2 * t4;
t79 = cc(1) * cc(3) * dd(2);
t81 = cc(3) ^ 2;
t85 = t2 * t12;
t87 = cc(1) * cc(2) * dd(3);
t89 = cc(2) ^ 2;
t93 = cc(1) ^ 2;
t97 = t2 * t21;
t100 = cc(2) * cc(3) * dd(1);
t113 = t57 * t15 - t60 * t7 - t63 * t64 * dd(3) + t68 * t64 * dd(2) + t72 * t28 - t68 * t74 * dd(1) + t77 * t79 - t56 * t81 * bb(1) * dd(2) - t85 * t87 + t59 * t89 * bb(1) * dd(3) - t71 * t93 * bb(2) * dd(3) + t97 * t87 - t77 * t100 + t56 * t81 * bb(2) * dd(1) - t59 * t89 * bb(3) * dd(1) + t85 * t100 - t97 * t79 + t71 * t93 * bb(3) * dd(2);
t126 = t2 * t38;
t127 = cc(1) * dd(3);
t129 = t2 * t34;
t130 = cc(1) * dd(2);
t132 = t2 * t48;
t133 = cc(2) * dd(3);
t135 = cc(3) * dd(2);
t137 = cc(3) * dd(1);
t139 = cc(2) * dd(1);
t141 = t77 * t7 + t77 * t10 - t85 * t15 - t85 * t18 - t77 * t15 + t97 * t24 + t97 * t10 - t77 * t28 - t97 * t18 - t97 * t28 + t85 * t7 + t85 * t24 + t126 * t127 - t129 * t130 - t132 * t133 + t132 * t135 - t126 * t137 + t129 * t139;
t142 = t113 / 0.6e1 + t141 / 0.12e2;
t143 = t3 * t81;
t150 = t3 * t89;
t155 = t3 * t93;
t166 = -t143 * t18 / 0.18e2 + t143 * t10 / 0.18e2 + t143 * t24 / 0.18e2 - t150 * t28 / 0.18e2 - t150 * t18 / 0.18e2 + t155 * t24 / 0.18e2 + t155 * t7 / 0.18e2 + t8 / 0.72e2 + t11 / 0.72e2 - t16 / 0.72e2 - t19 / 0.72e2 - t20 / 0.72e2 + t25 / 0.72e2;
t173 = t56 * t81;
t182 = t59 * t89;
t187 = t71 * t93;
t208 = t26 / 0.72e2 - t29 / 0.72e2 - t30 / 0.72e2 - t31 / 0.72e2 + t32 / 0.72e2 + t33 / 0.72e2 - t173 * bb(1) * dd(2) * t1 / 0.18e2 - t85 * cc(1) * t133 * t1 / 0.18e2 + t182 * bb(1) * dd(3) * t1 / 0.18e2 - t187 * bb(2) * dd(3) * t1 / 0.18e2 + t97 * cc(2) * t127 * t1 / 0.18e2 - t77 * cc(2) * t137 * t1 / 0.18e2 + t173 * bb(2) * dd(1) * t1 / 0.18e2 + t77 * cc(1) * t135 * t1 / 0.18e2;
t227 = dd(3) * t1;
t230 = dd(2) * t1;
t237 = -t182 * bb(3) * dd(1) * t1 + t85 * cc(3) * t139 * t1 - t97 * cc(3) * t130 * t1 + t187 * bb(3) * dd(2) * t1 + t150 * t7 + t150 * t10 - t143 * t15 - t155 * t15 - t155 * t28 + t57 * t14 * t227 - t60 * t6 * t230 - t63 * t64 * t227 + t68 * t64 * t230;
t238 = dd(1) * t1;
t251 = t81 * cc(3);
t256 = t93 * cc(1);
t261 = t89 * cc(2);
t278 = t72 * t27 * t238 / 0.18e2 - t68 * t74 * t238 / 0.18e2 + t37 / 0.72e2 + t41 / 0.72e2 - t44 / 0.72e2 - t47 / 0.72e2 - t51 / 0.72e2 + t54 / 0.72e2 + t3 * t251 * bb(1) * dd(2) / 0.18e2 - t3 * t256 * bb(3) * dd(2) / 0.18e2 - t3 * t261 * bb(1) * dd(3) / 0.18e2 + t3 * t256 * bb(2) * dd(3) / 0.18e2 - t3 * t251 * bb(2) * dd(1) / 0.18e2 + t3 * t261 * bb(3) * dd(1) / 0.18e2;
%
cm(1) = t55 / 0.24e2;
cm(2) = t142;
cm(3) = -t55 / 0.12e2;
cm(4) = t55 / 0.24e2;
cm(5) = t166 + t208 + t237 / 0.18e2 + t278;
cm(6) = -t142;
%
t1 = (ko ^ 2);
t3 = (j ^ 2);
t5 = (0.1e1 / t1 / t3);
t7 = 3 * t5 * cm(2);
t10 = (0.1e1 / ko / j);
t19 = 8 / t1 / ko / t3 / j * cm(5);
t21 = 2 * t10 * cm(1);
t23 = 8 * t5 * cm(5);
t25 = 2 * t10 * cm(3);
t27 = 4 * t10 * cm(5);
t29 = 2 * t10 * cm(4);
t31 = 3 * t5 * cm(6);
%
coefm(1) = t7;
coefm(2) = 3 * t10 * cm(6);
coefm(3) = cm(4);
coefm(4) = cm(3);
coefm(5) = -t19;
coefm(6) = t19;
coefm(7) = -t21;
coefm(8) = cm(6);
coefm(9) = -t23;
coefm(10) = t25;
coefm(11) = t27;
coefm(12) = -t27;
coefm(13) = t29;
coefm(14) = t31;
coefm(15) = -t7;
coefm(16) = t23;
coefm(17) = cm(1);
coefm(18) = cm(5);
coefm(19) = -cm(5);
coefm(20) = cm(2);
coefm(21) = 3 * t10 * cm(2);
coefm(22) = -t31;
coefm(23) = -t25;
coefm(24) = -t29;
coefm(25) = -t19;
coefm(26) = t19;
coefm(27) = t21;
function [coef] = coefficients_f3_f3(r1,r2,r3,r4,r5)
%% coefficients_f3_f3

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

% OUTPUT DATA
% coef   
%

coef  = zeros(8);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = sqrt(0.3e1);
t2 = t1 * cc(2);
t4 = t2 * dd(3) * ee(1);
t5 = t1 * cc(1);
t7 = t5 * ee(2) * dd(3);
t9 = t5 * dd(2) * ee(3);
t11 = t2 * ee(3) * dd(1);
t12 = t1 * cc(3);
t14 = t12 * ee(1) * dd(2);
t16 = t12 * dd(1) * ee(2);
t19 = t5 * bb(3) * ee(2);
t20 = bb(1) * ee(2);
t21 = t12 * t20;
t22 = bb(2) * ee(3);
t23 = t5 * t22;
t24 = bb(3) * ee(1);
t25 = t2 * t24;
t27 = t2 * bb(1) * ee(3);
t29 = t12 * bb(2) * ee(1);
t32 = cc(1) * ee(2) * dd(3);
t34 = cc(1) * dd(2) * ee(3);
t36 = cc(2) * ee(3) * dd(1);
t38 = cc(2) * dd(3) * ee(1);
t40 = cc(3) * ee(1) * dd(2);
t42 = cc(3) * dd(1) * ee(2);
t44 = cc(3) * bb(2);
t45 = t44 * ee(1);
t46 = cc(1) * bb(3);
t47 = t46 * ee(2);
t48 = cc(2) * bb(3);
t49 = t48 * ee(1);
t50 = cc(2) * bb(1);
t51 = t50 * ee(3);
t52 = cc(3) * bb(1);
t53 = t52 * ee(2);
t54 = cc(1) * bb(2);
t55 = t54 * ee(3);
t57 = dd(3) * t1;
t59 = t20 * t57 / 0.12e2;
t60 = bb(3) * dd(1);
t63 = t60 * ee(2) * t1 / 0.12e2;
t64 = bb(1) * dd(2);
t67 = t64 * ee(3) * t1 / 0.12e2;
t70 = t54 * t57 / 0.12e2;
t71 = dd(1) * t1;
t73 = t44 * t71 / 0.12e2;
t76 = t50 * t57 / 0.12e2;
t79 = t22 * t71 / 0.12e2;
t81 = dd(2) * t1;
t83 = t24 * t81 / 0.12e2;
t86 = t46 * t81 / 0.12e2;
t88 = bb(2) * dd(3);
t91 = t88 * ee(1) * t1 / 0.12e2;
t93 = t48 * t71 / 0.12e2;
t95 = t52 * t81 / 0.12e2;
t96 = t59 - t63 - t67 + t19 / 0.6e1 + t70 - t73 + t27 / 0.6e1 - t76 + t29 / 0.6e1 + t79 - t23 / 0.6e1 + t83 - t25 / 0.6e1 - t86 - t21 / 0.6e1 - t91 + t93 + t95;
t103 = t79 - t63 - t67 - t86 - t76 + t95 - t91 + t93 + t83 - t73 + t59 - t11 / 0.6e1 - t14 / 0.6e1 + t70 - t7 / 0.6e1 + t9 / 0.6e1 + t16 / 0.6e1 + t4 / 0.6e1;
t104 = t52 * dd(2);
t105 = t22 * dd(1);
t106 = t88 * ee(1);
t107 = t46 * dd(2);
t108 = t64 * ee(3);
t109 = t20 * dd(3);
t110 = t54 * dd(3);
t111 = t48 * dd(1);
t112 = t44 * dd(1);
t113 = t50 * dd(3);
t114 = t24 * dd(2);
t115 = t60 * ee(2);
t116 = -t104 - t105 + t106 + t107 + t108 - t109 - t110 - t111 + t112 + t113 - t114 + t115;
t129 = -t105 / 0.12e2 - t114 / 0.12e2 + t113 / 0.12e2 + t108 / 0.12e2 + t115 / 0.12e2 + t106 / 0.12e2 + t107 / 0.12e2 - t110 / 0.12e2 - t109 / 0.12e2 - t111 / 0.12e2 - t104 / 0.12e2 + t53 / 0.6e1;
t142 = t112 / 0.12e2 + t32 / 0.6e1 - t45 / 0.6e1 - t34 / 0.6e1 - t51 / 0.6e1 - t42 / 0.6e1 + t49 / 0.6e1 - t38 / 0.6e1 + t55 / 0.6e1 + t36 / 0.6e1 - t47 / 0.6e1 + t40 / 0.6e1;
%
coef(1) = t4 / 0.6e1 - t7 / 0.6e1 + t9 / 0.6e1 - t11 / 0.6e1 - t14 / 0.6e1 + t16 / 0.6e1;
coef(2) = t19 / 0.6e1 - t21 / 0.6e1 - t23 / 0.6e1 - t25 / 0.6e1 + t27 / 0.6e1 + t29 / 0.6e1;
coef(3) = t32 / 0.2e1 - t34 / 0.2e1 + t36 / 0.2e1 - t38 / 0.2e1 + t40 / 0.2e1 - t42 / 0.2e1;
coef(4) = -t45 / 0.2e1 - t47 / 0.2e1 + t49 / 0.2e1 - t51 / 0.2e1 + t53 / 0.2e1 + t55 / 0.2e1;
coef(5) = t96;
coef(6) = t103;
coef(7) = t116 / 0.4e1;
coef(8) = t129 + t142;
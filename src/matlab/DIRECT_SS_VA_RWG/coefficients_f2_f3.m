function [coef] = coefficients_f2_f3(r1,r2,r3,r4,r5)
%% coefficients_f2_f3

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

coef  = zeros(7);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(2) * ee(3);
t2 = sqrt(0.3e1);
t3 = dd(1) * t2;
t4 = t1 * t3;
t5 = bb(1) * dd(2);
t6 = ee(3) * t2;
t7 = t5 * t6;
t8 = bb(1) * ee(2);
t9 = dd(3) * t2;
t10 = t8 * t9;
t11 = bb(2) * dd(3);
t12 = ee(1) * t2;
t13 = t11 * t12;
t14 = bb(3) * ee(1);
t15 = dd(2) * t2;
t16 = t14 * t15;
t17 = bb(3) * dd(1);
t18 = ee(2) * t2;
t19 = t17 * t18;
t20 = -t4 + t7 - t10 + t13 - t16 + t19;
t21 = t2 * bb(2);
t22 = cc(3) * ee(1);
t24 = t2 * bb(1);
t27 = t2 * bb(3);
t30 = cc(2) * ee(3);
t34 = cc(1) * ee(2);
t37 = t1 * dd(1);
t38 = t14 * dd(2);
t39 = t5 * ee(3);
t40 = t8 * dd(3);
t41 = t11 * ee(1);
t42 = t17 * ee(2);
t43 = t37 + t38 - t39 + t40 - t41 - t42;
t45 = bb(3) * cc(1);
t48 = cc(1) * dd(2);
t54 = bb(1) * cc(3);
t57 = bb(3) * cc(2);
t63 = bb(2) * cc(3);
t68 = bb(2) * cc(1);
t71 = -t40 / 0.12e2 + t45 * dd(2) / 0.6e1 - t48 * ee(3) / 0.6e1 + t34 * dd(3) / 0.6e1 + t41 / 0.12e2 - t54 * dd(2) / 0.6e1 - t57 * dd(1) / 0.6e1 + t22 * dd(2) / 0.6e1 - t38 / 0.12e2 - t63 * ee(1) / 0.3e1 + t30 * dd(1) / 0.6e1 - t68 * dd(3) / 0.6e1;
t76 = bb(1) * cc(2);
t82 = cc(3) * dd(1);
t91 = cc(2) * dd(3);
t96 = t39 / 0.12e2 - t37 / 0.12e2 + t54 * ee(2) / 0.3e1 - t76 * ee(3) / 0.3e1 + t63 * dd(1) / 0.6e1 + t42 / 0.12e2 - t82 * ee(2) / 0.6e1 + t76 * dd(3) / 0.6e1 - t45 * ee(2) / 0.3e1 + t57 * ee(1) / 0.3e1 - t91 * ee(1) / 0.6e1 + t68 * ee(3) / 0.3e1;
t128 = -t63 * t3 / 0.6e1 + t54 * t15 / 0.6e1 - t45 * t15 / 0.6e1 - t34 * t9 / 0.6e1 - t76 * t9 / 0.6e1 + t91 * t12 / 0.6e1 + t68 * t9 / 0.6e1 + t10 / 0.12e2 - t30 * t3 / 0.6e1 - t19 / 0.12e2 - t22 * t15 / 0.6e1 - t7 / 0.12e2 + t48 * t6 / 0.6e1 + t16 / 0.12e2 + t82 * t18 / 0.6e1 + t4 / 0.12e2 + t57 * t3 / 0.6e1 - t13 / 0.12e2;
%
coef(1) = t20 / 0.6e1;
coef(2) = t21 * t22 / 0.3e1 - t24 * cc(3) * ee(2) / 0.3e1 - t27 * cc(2) * ee(1) / 0.3e1 + t24 * t30 / 0.3e1 - t21 * cc(1) * ee(3) / 0.3e1 + t27 * t34 / 0.3e1;
coef(3) = t43 / 0.2e1;
coef(4) = t71 + t96;
coef(5) = -t43 / 0.4e1;
coef(6) = -t20 / 0.12e2;
coef(7) = t128;
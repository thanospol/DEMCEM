function [coef] = coefficients_f1_f3(r1,r2,r3,r4,r5)
%% coefficients_f1_f3

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

coef  = zeros(4);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(1) * ee(2);
t2 = sqrt(0.3e1);
t3 = dd(3) * t2;
t4 = t1 * t3;
t5 = bb(2) * ee(3);
t6 = dd(1) * t2;
t7 = t5 * t6;
t8 = bb(1) * dd(2);
t9 = ee(3) * t2;
t10 = t8 * t9;
t11 = bb(3) * ee(1);
t12 = dd(2) * t2;
t13 = t11 * t12;
t14 = bb(2) * dd(3);
t15 = ee(1) * t2;
t16 = t14 * t15;
t17 = bb(3) * dd(1);
t18 = ee(2) * t2;
t19 = t17 * t18;
t21 = t5 * dd(1);
t22 = t17 * ee(2);
t23 = t1 * dd(3);
t24 = t8 * ee(3);
t25 = t14 * ee(1);
t26 = t11 * dd(2);
t28 = cc(2) * ee(3);
t33 = cc(1) * ee(2);
t36 = cc(3) * ee(1);
t40 = cc(2) * dd(3);
t44 = cc(1) * dd(2);
t47 = cc(3) * dd(1);
t52 = t28 * dd(1) / 0.6e1 + t25 / 0.12e2 + t22 / 0.12e2 + t33 * dd(3) / 0.6e1 + t36 * dd(2) / 0.6e1 - t21 / 0.12e2 - t40 * ee(1) / 0.6e1 - t23 / 0.12e2 - t44 * ee(3) / 0.6e1 - t47 * ee(2) / 0.6e1 - t26 / 0.12e2 + t24 / 0.12e2;
t71 = t4 / 0.12e2 + t13 / 0.12e2 - t19 / 0.12e2 - t33 * t3 / 0.6e1 + t44 * t9 / 0.6e1 - t28 * t6 / 0.6e1 + t40 * t15 / 0.6e1 - t10 / 0.12e2 - t36 * t12 / 0.6e1 + t7 / 0.12e2 + t47 * t18 / 0.6e1 - t16 / 0.12e2;
%
coef(1) = t4 / 0.12e2 + t7 / 0.12e2 - t10 / 0.12e2 + t13 / 0.12e2 - t16 / 0.12e2 - t19 / 0.12e2;
coef(2) = -t21 / 0.4e1 + t22 / 0.4e1 - t23 / 0.4e1 + t24 / 0.4e1 + t25 / 0.4e1 - t26 / 0.4e1;
coef(3) = t52;
coef(4) = t71;
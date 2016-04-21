function [coef] = coefficients_f3_f1(r1,r2,r3,r4,r5)
%% coefficients_f3_f1

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
t1 = cc(1) * bb(2);
t2 = sqrt(0.3e1);
t3 = dd(3) * t2;
t4 = t1 * t3;
t6 = cc(1) * bb(3);
t7 = ee(2) * t2;
t10 = cc(3) * bb(2);
t11 = dd(1) * t2;
t12 = t10 * t11;
t14 = dd(2) * t2;
t15 = t6 * t14;
t17 = ee(3) * t2;
t20 = cc(2) * bb(3);
t21 = t20 * t11;
t23 = ee(1) * t2;
t26 = cc(2) * bb(1);
t27 = t26 * t3;
t29 = cc(3) * bb(1);
t32 = t29 * t14;
t38 = t4 / 0.12e2 + t6 * t7 / 0.6e1 - t12 / 0.12e2 - t15 / 0.12e2 - t1 * t17 / 0.6e1 + t21 / 0.12e2 - t20 * t23 / 0.6e1 - t27 / 0.12e2 - t29 * t7 / 0.6e1 + t32 / 0.12e2 + t10 * t23 / 0.6e1 + t26 * t17 / 0.6e1;
t41 = t1 * dd(3);
t45 = t10 * dd(1);
t47 = t20 * dd(1);
t49 = t6 * dd(2);
t57 = t29 * dd(2);
t59 = t26 * dd(3);
t63 = t1 * ee(3) / 0.6e1 - t41 / 0.12e2 - t26 * ee(3) / 0.6e1 + t45 / 0.12e2 - t47 / 0.12e2 + t49 / 0.12e2 + t29 * ee(2) / 0.6e1 + t20 * ee(1) / 0.6e1 - t10 * ee(1) / 0.6e1 - t57 / 0.12e2 + t59 / 0.12e2 - t6 * ee(2) / 0.6e1;
%
coef(1) = t38;
coef(2) = t63;
coef(3) = t59 / 0.4e1 + t45 / 0.4e1 - t57 / 0.4e1 - t41 / 0.4e1 - t47 / 0.4e1 + t49 / 0.4e1;
coef(4) = -t27 / 0.12e2 + t21 / 0.12e2 + t4 / 0.12e2 - t12 / 0.12e2 - t15 / 0.12e2 + t32 / 0.12e2;
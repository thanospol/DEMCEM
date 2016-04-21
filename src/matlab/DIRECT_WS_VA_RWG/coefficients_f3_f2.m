function [coef] = coefficients_f3_f2(r1,r2,r3,r4,r5)
%% coefficients_f3_f2

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

coef  = zeros(6);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = cc(1) * dd(1);
t2 = cc(2) * dd(2);
t3 = cc(3) * dd(3);
t4 = t1 + t2 + t3;
t5 = bb(2) * dd(2);
t6 = bb(1) * dd(1);
t7 = bb(3) * dd(3);
t8 = -t5 - t6 - t7;
t9 = sqrt(0.3e1);
t10 = t2 * t9;
t11 = t10 / 0.6e1;
t12 = t1 * t9;
t13 = t12 / 0.6e1;
t17 = t3 * t9;
t18 = t17 / 0.6e1;
t26 = t6 * t9;
t30 = t7 * t9;
t33 = t9 * bb(2) * dd(2);
t39 = bb(1) * ee(1);
t47 = bb(3) * ee(3);
t49 = bb(2) * ee(2);
t55 = t5 / 0.12e2 + t7 / 0.12e2 - t39 / 0.6e1 + t6 / 0.12e2 + cc(2) * ee(2) / 0.3e1 + cc(3) * ee(3) / 0.3e1 - t2 / 0.6e1 - t47 / 0.6e1 - t49 / 0.6e1 - t1 / 0.6e1 - t3 / 0.6e1 + cc(1) * ee(1) / 0.3e1;
t56 = t33 / 0.12e2;
t57 = t26 / 0.12e2;
t58 = t30 / 0.12e2;
%
coef(1) = t4;
coef(2) = t8 / 0.2e1;
coef(3) = -t4 / 0.2e1;
coef(4) = t11 + t13 - t9 * cc(1) * ee(1) / 0.3e1 + t18 - t9 * cc(3) * ee(3) / 0.3e1 - t9 * cc(2) * ee(2) / 0.3e1;
coef(5) = t26 / 0.6e1 - t10 / 0.3e1 - t12 / 0.3e1 + t30 / 0.6e1 + t33 / 0.6e1 - t17 / 0.3e1;
coef(6) = t55;
coef(7) = -t8 / 0.4e1;
coef(8) = t18 + t13 + t11 - t56 - t57 - t58;
coef(9) = -t56 - t57 + t49 * t9 / 0.6e1 + t39 * t9 / 0.6e1 - t58 + t47 * t9 / 0.6e1;
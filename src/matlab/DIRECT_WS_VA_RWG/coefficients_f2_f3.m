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

coef  = zeros(9);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(3) * ee(3);
t2 = bb(1) * ee(1);
t3 = bb(2) * ee(2);
t4 = t1 + t2 + t3;
t5 = sqrt(0.3e1);
t6 = t5 * bb(3);
t7 = t6 * dd(3);
t8 = t7 / 0.12e2;
t9 = cc(1) * dd(1);
t12 = t5 * bb(2);
t13 = t12 * dd(2);
t14 = t13 / 0.12e2;
t15 = cc(3) * dd(3);
t18 = t5 * bb(1);
t19 = t18 * dd(1);
t20 = t19 / 0.12e2;
t21 = cc(2) * dd(2);
t25 = t6 * ee(3);
t26 = t25 / 0.6e1;
t27 = t18 * ee(1);
t28 = t27 / 0.6e1;
t29 = t12 * ee(2);
t30 = t29 / 0.6e1;
t32 = bb(1) * dd(1);
t33 = bb(3) * dd(3);
t34 = bb(2) * dd(2);
t35 = t32 + t33 + t34;
t68 = cc(1) * ee(1) / 0.3e1 + t32 / 0.12e2 + t34 / 0.12e2 + t33 / 0.12e2 - t1 / 0.6e1 - t2 / 0.6e1 - t3 / 0.6e1 - t15 / 0.6e1 - t21 / 0.6e1 - t9 / 0.6e1 + cc(2) * ee(2) / 0.3e1 + cc(3) * ee(3) / 0.3e1;
%
coef(1) = t4;
coef(2) = -t8 + t9 * t5 / 0.6e1 - t14 + t15 * t5 / 0.6e1 - t20 + t21 * t5 / 0.6e1;
coef(3) = t26 + t28 - t8 - t20 + t30 - t14;
coef(4) = t35 / 0.4e1;
coef(5) = -t4 / 0.2e1;
coef(6) = -t35 / 0.2e1;
coef(7) = t7 / 0.6e1 + t13 / 0.6e1 + t19 / 0.6e1 - t25 / 0.3e1 - t27 / 0.3e1 - t29 / 0.3e1;
coef(8) = -t5 * cc(1) * ee(1) / 0.3e1 - t5 * cc(2) * ee(2) / 0.3e1 + t28 - t5 * cc(3) * ee(3) / 0.3e1 + t30 + t26;
coef(9) = t68;
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

coef  = zeros(6);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = cc(1) * ee(1);
t2 = cc(2) * ee(2);
t3 = cc(3) * ee(3);
t5 = cc(1) * dd(1);
t6 = cc(2) * dd(2);
t7 = cc(3) * dd(3);
t9 = bb(3) * ee(3);
t10 = bb(1) * ee(1);
t11 = bb(2) * ee(2);
t13 = sqrt(0.3e1);
t16 = t13 * bb(2) * ee(2) / 0.6e1;
t19 = t13 * bb(3) * ee(3) / 0.6e1;
t20 = t13 * cc(1);
t22 = t20 * ee(1) / 0.3e1;
t23 = t13 * cc(2);
t25 = t23 * ee(2) / 0.3e1;
t28 = t13 * bb(1) * ee(1) / 0.6e1;
t29 = t13 * cc(3);
t31 = t29 * ee(3) / 0.3e1;
t33 = bb(1) * dd(1);
t35 = bb(3) * dd(3);
t40 = bb(2) * dd(2);
t48 = t33 / 0.12e2 + t35 / 0.12e2 - t10 / 0.6e1 - t5 / 0.6e1 - t9 / 0.6e1 + t40 / 0.12e2 + t3 / 0.3e1 - t6 / 0.6e1 - t7 / 0.6e1 + t1 / 0.3e1 + t2 / 0.3e1 - t11 / 0.6e1;
t51 = t29 * dd(3) / 0.6e1;
t53 = t23 * dd(2) / 0.6e1;
t55 = t35 * t13 / 0.12e2;
t57 = t33 * t13 / 0.12e2;
t59 = t40 * t13 / 0.12e2;
t61 = t20 * dd(1) / 0.6e1;
%
coef(1) = t1 + t2 + t3;
coef(2) = -t5 / 0.2e1 - t6 / 0.2e1 - t7 / 0.2e1;
coef(3) = -t9 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1;
coef(4) = t16 + t19 - t22 - t25 + t28 - t31;
coef(5) = t48;
coef(6) = t35 / 0.4e1 + t33 / 0.4e1 + t40 / 0.4e1;
coef(7) = t51 + t53 - t55 - t57 - t59 + t61;
coef(8) = -t59 + t16 - t55 + t19 + t28 - t57;
coef(9) = t51 - t25 + t61 + t53 - t31 - t22;
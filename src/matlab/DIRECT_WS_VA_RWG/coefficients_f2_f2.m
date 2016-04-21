function [coef] = coefficients_f2_f2(r1,r2,r3,r4,r5)
%% coefficients_f2_f2

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
t1 = bb(1) * dd(1);
t2 = bb(2) * dd(2);
t3 = bb(3) * dd(3);
t4 = t1 + t2 + t3;
t26 = t1 / 0.12e2 - bb(1) * ee(1) / 0.6e1 - bb(3) * ee(3) / 0.6e1 + t3 / 0.12e2 - cc(1) * dd(1) / 0.6e1 - bb(2) * ee(2) / 0.6e1 + cc(2) * ee(2) / 0.3e1 - cc(2) * dd(2) / 0.6e1 + cc(3) * ee(3) / 0.3e1 - cc(3) * dd(3) / 0.6e1 + t2 / 0.12e2 + cc(1) * ee(1) / 0.3e1;
t27 = sqrt(0.3e1);
t29 = t27 * cc(3) * dd(3);
t31 = t27 * bb(1);
t32 = t31 * dd(1);
t33 = t32 / 0.6e1;
t34 = t27 * bb(3);
t35 = t34 * dd(3);
t36 = t35 / 0.6e1;
t37 = t27 * bb(2);
t38 = t37 * dd(2);
t39 = t38 / 0.6e1;
t41 = t27 * cc(1) * dd(1);
t44 = t27 * cc(2) * dd(2);
t47 = t31 * ee(1);
t49 = t37 * ee(2);
t51 = t34 * ee(3);
t57 = t32 / 0.12e2;
t58 = t35 / 0.12e2;
t59 = t38 / 0.12e2;
%
coef(1) = t4;
coef(2) = t26;
coef(3) = -t29 / 0.3e1 + t33 + t36 + t39 - t41 / 0.3e1 - t44 / 0.3e1;
coef(4) = -t4 / 0.2e1;
coef(5) = t39 - t47 / 0.3e1 + t33 + t36 - t49 / 0.3e1 - t51 / 0.3e1;
coef(6) = -t4 / 0.2e1;
coef(7) = t4 / 0.4e1;
coef(8) = t29 / 0.6e1 + t44 / 0.6e1 + t41 / 0.6e1 - t57 - t58 - t59;
coef(9) = -t57 + t51 / 0.6e1 - t58 + t49 / 0.6e1 + t47 / 0.6e1 - t59;
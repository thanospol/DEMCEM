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

coef  = zeros(6);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(1) * dd(1);
t2 = sqrt(0.3e1);
t4 = t1 * t2 / 0.12e2;
t5 = bb(1) * ee(1);
t8 = bb(3) * dd(3);
t10 = t8 * t2 / 0.12e2;
t11 = bb(3) * ee(3);
t14 = bb(2) * dd(2);
t16 = t14 * t2 / 0.12e2;
t17 = bb(2) * ee(2);
t21 = t2 * cc(1);
t23 = t21 * dd(1) / 0.6e1;
t24 = t2 * cc(2);
t26 = t24 * dd(2) / 0.6e1;
t27 = t2 * cc(3);
t29 = t27 * dd(3) / 0.6e1;
t38 = cc(2) * dd(2);
t42 = cc(1) * dd(1);
t46 = cc(3) * dd(3);
%
t50 = t8 / 0.12e2 - t5 / 0.6e1 - t17 / 0.6e1 + cc(2) * ee(2) / 0.3e1 + t1 / 0.12e2 - t38 / 0.6e1 + t14 / 0.12e2 - t11 / 0.6e1 - t42 / 0.6e1 + cc(3) * ee(3) / 0.3e1 - t46 / 0.6e1 + cc(1) * ee(1) / 0.3e1;
coef(1) = -t4 + t5 * t2 / 0.6e1 - t10 + t11 * t2 / 0.6e1 - t16 + t17 * t2 / 0.6e1;
coef(2) = -t4 - t10 + t23 + t26 - t16 + t29;
coef(3) = t1 / 0.4e1 + t14 / 0.4e1 + t8 / 0.4e1;
coef(4) = t50;
coef(5) = t29 - t21 * ee(1) / 0.3e1 - t24 * ee(2) / 0.3e1 + t26 + t23 - t27 * ee(3) / 0.3e1;
coef(6) = -t46 / 0.2e1 - t42 / 0.2e1 - t38 / 0.2e1;
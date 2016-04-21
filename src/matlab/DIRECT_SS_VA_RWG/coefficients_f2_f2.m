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

coef  = zeros(5);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(2) * dd(3);
t2 = sqrt(0.3e1);
t5 = bb(1) * ee(2);
t6 = dd(3) * t2;
t8 = bb(2) * ee(3);
t9 = dd(1) * t2;
t11 = bb(3) * ee(1);
t12 = dd(2) * t2;
t14 = bb(1) * dd(2);
t17 = bb(3) * dd(1);
t20 = -t1 * ee(1) * t2 + t5 * t6 + t8 * t9 + t11 * t12 - t14 * ee(3) * t2 - t17 * ee(2) * t2;
t21 = bb(2) * cc(3);
t23 = bb(1) * cc(2);
t25 = bb(2) * cc(1);
t27 = bb(3) * cc(2);
t29 = bb(3) * cc(1);
t31 = bb(1) * cc(3);
t33 = -t21 * t9 - t23 * t6 + t25 * t6 + t27 * t9 - t29 * t12 + t31 * t12;
t58 = -t5 * dd(3) / 0.6e1 - t23 * ee(3) / 0.3e1 + t1 * ee(1) / 0.6e1 - t11 * dd(2) / 0.6e1 - t25 * dd(3) / 0.6e1 - t27 * dd(1) / 0.6e1 + t29 * dd(2) / 0.6e1 - t31 * dd(2) / 0.6e1 + t31 * ee(2) / 0.3e1 + t21 * dd(1) / 0.6e1 + t17 * ee(2) / 0.6e1 - t8 * dd(1) / 0.6e1;
t89 = t14 * ee(3) / 0.6e1 + t23 * dd(3) / 0.6e1 + cc(3) * ee(1) * dd(2) / 0.3e1 - cc(2) * dd(3) * ee(1) / 0.3e1 - t29 * ee(2) / 0.3e1 + t27 * ee(1) / 0.3e1 + cc(2) * ee(3) * dd(1) / 0.3e1 - t21 * ee(1) / 0.3e1 + t25 * ee(3) / 0.3e1 + cc(1) * ee(2) * dd(3) / 0.3e1 - cc(3) * dd(1) * ee(2) / 0.3e1 - cc(1) * dd(2) * ee(3) / 0.3e1;
%
coef(1) = t20 / 0.6e1;
coef(2) = t33 / 0.6e1;
coef(3) = t58 + t89;
coef(4) = -t33 / 0.3e1;
coef(5) = -t20 / 0.3e1;
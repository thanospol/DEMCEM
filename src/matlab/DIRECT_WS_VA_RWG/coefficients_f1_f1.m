function [coef] = coefficients_f1_f1(r1,r2,r3,r4,r5)
%% coefficients_f1_f1

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
t1 = bb(1) * dd(1);
t3 = bb(2) * dd(2);
t5 = bb(3) * dd(3);
t7 = bb(1) * ee(1);
t9 = bb(2) * ee(2);
t11 = bb(3) * ee(3);
t13 = cc(1) * dd(1);
t15 = cc(2) * dd(2);
t17 = cc(3) * dd(3);
t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc(1) * ee(1) / 0.3e1 + cc(2) * ee(2) / 0.3e1 + cc(3) * ee(3) / 0.3e1;
t27 = sqrt(0.3e1);
t29 = t1 * t27 / 0.12e2;
t31 = t3 * t27 / 0.12e2;
t35 = t5 * t27 / 0.12e2;
%
coef(1) = t25;
coef(2) = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
coef(3) = -t29 - t31 + t9 * t27 / 0.6e1 - t35 + t7 * t27 / 0.6e1 + t11 * t27 / 0.6e1;
coef(4) = -t29 - t35 - t31 + t13 * t27 / 0.6e1 + t15 * t27 / 0.6e1 + t17 * t27 / 0.6e1;
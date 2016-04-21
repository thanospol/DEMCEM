function [coef] = coefficients_f2_f1(r1,r2,r3,r4,r5)
%% coefficients_f2_f1

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

coef  = zeros(2);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(3) * cc(2);
t4 = bb(3) * cc(1);
t7 = bb(1) * cc(3);
t10 = bb(1) * cc(2);
t13 = bb(2) * cc(1);
t16 = bb(2) * cc(3);
t31 = -t1 * dd(1) / 0.6e1 + t4 * dd(2) / 0.6e1 - t7 * dd(2) / 0.6e1 + t10 * dd(3) / 0.6e1 - t13 * dd(3) / 0.6e1 + t16 * dd(1) / 0.6e1 - t4 * ee(2) / 0.3e1 + t1 * ee(1) / 0.3e1 - t10 * ee(3) / 0.3e1 + t7 * ee(2) / 0.3e1 - t16 * ee(1) / 0.3e1 + t13 * ee(3) / 0.3e1;
t32 = sqrt(0.3e1);
t33 = dd(1) * t32;
t35 = dd(2) * t32;
t38 = dd(3) * t32;
%
coef(1) = t31;
coef(2) = t1 * t33 / 0.6e1 - t4 * t35 / 0.6e1 + t7 * t35 / 0.6e1 - t10 * t38 / 0.6e1 + t13 * t38 / 0.6e1 - t16 * t33 / 0.6e1;
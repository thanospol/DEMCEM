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

coef  = zeros(7);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
ee = r5-r1;
%
t1 = bb(2) * ee(1);
t4 = cc(1) * bb(3);
t7 = cc(3) * bb(2);
t8 = t7 * dd(1);
t10 = bb(2) * ee(3);
t13 = cc(2) * bb(1);
t16 = t4 * dd(2);
t18 = cc(3) * bb(1);
t21 = t18 * dd(2);
t23 = t13 * dd(3);
t25 = cc(1) * bb(2);
t26 = t25 * dd(3);
t31 = bb(1) * ee(2);
t34 = t1 * dd(3) / 0.6e1 - t4 * ee(2) / 0.6e1 + t8 / 0.12e2 - t10 * dd(1) / 0.6e1 - t13 * ee(3) / 0.6e1 + t16 / 0.12e2 + t18 * ee(2) / 0.6e1 - t21 / 0.12e2 + t23 / 0.12e2 - t26 / 0.12e2 + cc(1) * ee(2) * dd(3) / 0.3e1 - t31 * dd(3) / 0.6e1;
t35 = bb(3) * ee(1);
t41 = bb(3) * ee(2);
t44 = cc(2) * bb(3);
t45 = t44 * dd(1);
t47 = bb(1) * ee(3);
t68 = -t35 * dd(2) / 0.6e1 + cc(2) * ee(3) * dd(1) / 0.3e1 + t41 * dd(1) / 0.6e1 - t45 / 0.12e2 + t47 * dd(2) / 0.6e1 - t7 * ee(1) / 0.6e1 + t25 * ee(3) / 0.6e1 - cc(2) * ee(1) * dd(3) / 0.3e1 + t44 * ee(1) / 0.6e1 + cc(3) * ee(1) * dd(2) / 0.3e1 - cc(1) * ee(3) * dd(2) / 0.3e1 - cc(3) * ee(2) * dd(1) / 0.3e1;
t70 = t23 + t8 - t21 + t16 - t26 - t45;
t71 = sqrt(0.3e1);
t72 = dd(3) * t71;
t75 = dd(2) * t71;
t78 = t71 * cc(2);
t80 = t78 * bb(1) * dd(3);
t82 = t71 * cc(1);
t84 = t82 * bb(2) * dd(3);
t88 = t71 * cc(3);
t90 = t88 * bb(2) * dd(1);
t93 = t88 * bb(1) * dd(2);
t95 = ee(1) * t71;
t99 = t78 * bb(3) * dd(1);
t101 = ee(2) * t71;
t105 = t82 * bb(3) * dd(2);
t107 = ee(3) * t71;
t110 = dd(1) * t71;
t123 = t31 * t72 / 0.6e1 + t35 * t75 / 0.6e1 - t80 / 0.12e2 + t84 / 0.12e2 - t47 * t75 / 0.6e1 - t90 / 0.12e2 + t93 / 0.12e2 + t7 * t95 / 0.6e1 + t99 / 0.12e2 - t18 * t101 / 0.6e1 - t105 / 0.12e2 + t13 * t107 / 0.6e1 + t10 * t110 / 0.6e1 - t1 * t72 / 0.6e1 - t44 * t95 / 0.6e1 - t41 * t110 / 0.6e1 + t4 * t101 / 0.6e1 - t25 * t107 / 0.6e1;
t124 = -t80 + t84 - t90 + t93 + t99 - t105;
%
coef(1) = t34 + t68;
coef(2) = t70 / 0.4e1;
coef(3) = t123;
coef(4) = t124 / 0.12e2;
coef(5) = t78 * ee(1) * dd(3) / 0.3e1 - t88 * ee(1) * dd(2) / 0.3e1 + t82 * ee(3) * dd(2) / 0.3e1 + t88 * ee(2) * dd(1) / 0.3e1 - t82 * ee(2) * dd(3) / 0.3e1 - t78 * ee(3) * dd(1) / 0.3e1;
coef(6) = -t124 / 0.6e1;
coef(7) = -t70 / 0.2e1;
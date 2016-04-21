function [coef,coefm] = coefficients_f3_f2(r1,r2,r3,r4,ko)
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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% ko = wavenumber

% OUTPUT DATA
% coef   
% coefm
%

c  = zeros(3);
cm = zeros(3);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
% 
t1 = cc(1) * bb(2);
t2 = sqrt(0.3e1);
t3 = dd(3) * t2;
t5 = cc(2) * bb(1);
t7 = cc(3) * bb(2);
t8 = dd(1) * t2;
t10 = cc(1) * bb(3);
t11 = dd(2) * t2;
t13 = cc(3) * bb(1);
t15 = cc(2) * bb(3);
t17 = -t1 * t3 + t5 * t3 + t7 * t8 + t10 * t11 - t13 * t11 - t15 * t8;
%
c(1) = t17 / 0.6e1;
c(2) = -t5 * dd(3) / 0.6e1 - t10 * dd(2) / 0.6e1 + t13 * dd(2) / 0.6e1 + t1 * dd(3) / 0.6e1 + t15 * dd(1) / 0.6e1 - t7 * dd(1) / 0.6e1;
c(3) = t17 / 0.6e1;
%
t3 = (0.1e1 / ko / j);
t7 = 2 * t3 * c(1);
t8 = (ko ^ 2);
t10 = (j ^ 2);
t14 = 3 / t8 / t10 * c(2);
t16 = 2 * t3 * c(3);
%
coef(1) = c(2);
coef(2) = c(1);
coef(3) = 3 * t3 * c(2);
coef(4) = -t7;
coef(5) = -t14;
coef(6) = t16;
coef(7) = c(3);
coef(8) = t7;
coef(9) = t14;
coef(10) = -t16;
%
t1 = cc(1) * bb(2);
t2 = sqrt(0.3e1);
t3 = dd(3) * t2;
t5 = cc(2) * bb(1);
t7 = cc(3) * bb(2);
t8 = dd(1) * t2;
t10 = cc(1) * bb(3);
t11 = dd(2) * t2;
t13 = cc(3) * bb(1);
t15 = cc(2) * bb(3);
t17 = -t1 * t3 + t5 * t3 + t7 * t8 + t10 * t11 - t13 * t11 - t15 * t8;
%
cm(1) = t17 / 0.6e1;
cm(2) = -t17 / 0.6e1;
cm(3) = -t5 * dd(3) / 0.6e1 - t10 * dd(2) / 0.6e1 + t13 * dd(2) / 0.6e1 + t1 * dd(3) / 0.6e1 + t15 * dd(1) / 0.6e1 - t7 * dd(1) / 0.6e1;
%
t3 = (0.1e1 / ko / j);
t5 = 2 * t3 * cm(2);
t6 = (ko ^ 2);
t8 = (j ^ 2);
t12 = 3 / t6 / t8 * cm(3);
t14 = 2 * t3 * cm(1);
%
coefm(1) = t5;
coefm(2) = t12;
coefm(3) = cm(1);
coefm(4) = -t14;
coefm(5) = cm(3);
coefm(6) = cm(2);
coefm(7) = 3 * t3 * cm(3);
coefm(8) = -t5;
coefm(9) = t14;
coefm(10) = -t12;
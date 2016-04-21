function [coef,coefm] = coefficients_f2_f2(r1,r2,r3,r4,ko)
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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% ko = wavenumber

% OUTPUT DATA
% coef   
% coefm
%

c  = zeros(10);
cm = zeros(10);
%
coef  = zeros(40);
coefm = zeros(40);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
%
t1 = bb(1) ^ 2;
t2 = bb(2) ^ 2;
t3 = bb(3) ^ 2;
t4 = -t1 - t2 - t3;
t5 = bb(1) * dd(1);
t6 = sqrt(0.3e1);
t9 = bb(3) * dd(3);
t13 = t2 * t6 / 0.12e2;
t15 = t1 * t6 / 0.12e2;
t17 = t3 * t6 / 0.12e2;
t18 = bb(2) * dd(2);
t21 = -t5 * t6 / 0.6e1 - t9 * t6 / 0.6e1 + t13 + t15 + t17 - t18 * t6 / 0.6e1;
t22 = cc(2) * bb(2);
t25 = cc(3) * bb(3);
t28 = cc(1) * bb(1);
t31 = t22 * t6 / 0.6e1 + t25 * t6 / 0.6e1 - t13 - t17 - t15 + t28 * t6 / 0.6e1;
t47 = cc(3) * dd(3) / 0.3e1 + cc(1) * dd(1) / 0.3e1 - t9 / 0.6e1 + t2 / 0.12e2 - t22 / 0.6e1 + t3 / 0.12e2 - t18 / 0.6e1 + t1 / 0.12e2 + cc(2) * dd(2) / 0.3e1 - t5 / 0.6e1 - t25 / 0.6e1 - t28 / 0.6e1;
%
c(1) = t4 / 0.4e1;
c(2) = -t4 / 0.4e1;
c(3) = t21;
c(4) = t31;
c(5) = t47;
c(6) = -t4 / 0.4e1;
c(7) = -t31;
c(8) = -t21;
c(9) = t4 / 0.4e1;
c(10) = t31;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * c(10);
t11 = 2 * t7 * c(2);
t14 = (0.1e1 / t1 / t4);
t15 = t14 * c(1);
t17 = 2 * t7 * c(3);
t20 = (0.1e1 / ko / j);
t41 = t14 * c(6);
t42 = t1 ^ 2;
t44 = t4 ^ 2;
t46 = (0.1e1 / t42 / t44);
t48 = 6 * t46 * c(5);
t50 = 6 * t46 * c(7);
t52 = 2 * t7 * c(4);
t54 = 2 * t7 * c(8);
t56 = 2 * t7 * c(9);
%
coef(1) = t9;
coef(2) = t11;
coef(3) = -t15;
coef(4) = t17;
coef(5) = -t20 * c(9);
coef(6) = -t20 * c(5);
coef(7) = -2 * t14 * c(4);
coef(8) = -t20 * c(4);
coef(9) = -t20 * c(8);
coef(10) = -t20 * c(7);
coef(11) = -2 * t14 * c(8);
coef(12) = -2 * t14 * c(9);
coef(13) = -6 * t7 * c(5);
coef(14) = -6 * t7 * c(7);
coef(15) = -3 * t14 * c(5);
coef(16) = -3 * t14 * c(7);
coef(17) = -t9;
coef(18) = -t20 * c(1);
coef(19) = -t41;
coef(20) = -t17;
coef(21) = -t11;
coef(22) = t48;
coef(23) = t50;
coef(24) = t52;
coef(25) = t54;
coef(26) = t56;
coef(27) = t15;
coef(28) = -2 * t14 * c(2);
coef(29) = -t50;
coef(30) = -t54;
coef(31) = -t20 * c(3);
coef(32) = -t20 * c(2);
coef(33) = -2 * t14 * c(10);
coef(34) = -t52;
coef(35) = -2 * t14 * c(3);
coef(36) = -t48;
coef(37) = -t56;
coef(38) = -t20 * c(10);
coef(39) = -t20 * c(6);
coef(40) = t41;
%
t1 = bb(1) ^ 2;
t2 = bb(2) ^ 2;
t3 = bb(3) ^ 2;
t4 = -t1 - t2 - t3;
t5 = sqrt(0.3e1);
t7 = t3 * t5 / 0.12e2;
t8 = bb(3) * dd(3);
t12 = t2 * t5 / 0.12e2;
t14 = t1 * t5 / 0.12e2;
t15 = bb(2) * dd(2);
t18 = bb(1) * dd(1);
t21 = t7 - t8 * t5 / 0.6e1 + t12 + t14 - t15 * t5 / 0.6e1 - t18 * t5 / 0.6e1;
t22 = cc(1) * bb(1);
t25 = cc(3) * bb(3);
t28 = cc(2) * bb(2);
t31 = t22 * t5 / 0.6e1 + t25 * t5 / 0.6e1 - t7 - t14 - t12 + t28 * t5 / 0.6e1;
t47 = -t15 / 0.6e1 - t28 / 0.6e1 - t18 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + t3 / 0.12e2 + cc(3) * dd(3) / 0.3e1 - t25 / 0.6e1 + t1 / 0.12e2 - t22 / 0.6e1 - t8 / 0.6e1 + cc(1) * dd(1) / 0.3e1 + t2 / 0.12e2;
%
cm(1) = t4 / 0.4e1;
cm(2) = t21;
cm(3) = t4 / 0.4e1;
cm(4) = t21;
cm(5) = -t4 / 0.4e1;
cm(6) = t31;
cm(7) = t4 / 0.4e1;
cm(8) = t31;
cm(9) = -t31;
cm(10) = t47;
%
t1 = ko ^ 2;
t2 = t1 ^ 2;
t4 = j ^ 2;
t5 = t4 ^ 2;
t7 = (0.1e1 / t2 / t5);
t9 = 6 * t7 * cm(8);
t14 = (0.1e1 / t1 / ko / t4 / j);
t16 = 2 * t14 * cm(9);
t18 = 2 * t14 * cm(2);
t20 = 2 * t14 * cm(7);
t23 = (0.1e1 / t1 / t4);
t24 = t23 * cm(5);
t27 = (0.1e1 / ko / j);
t30 = 2 * t14 * cm(6);
t32 = 2 * t14 * cm(4);
t34 = 2 * t14 * cm(3);
t36 = 6 * t7 * cm(10);
t66 = t23 * cm(1);
%
coefm(1) = t9;
coefm(2) = t16;
coefm(3) = t18;
coefm(4) = t20;
coefm(5) = -t24;
coefm(6) = -t27 * cm(1);
coefm(7) = -t30;
coefm(8) = -t32;
coefm(9) = -t34;
coefm(10) = t36;
coefm(11) = -t20;
coefm(12) = -t9;
coefm(13) = -t27 * cm(5);
coefm(14) = -t36;
coefm(15) = -t27 * cm(4);
coefm(16) = -t27 * cm(6);
coefm(17) = -2 * t23 * cm(6);
coefm(18) = -t27 * cm(3);
coefm(19) = -2 * t23 * cm(3);
coefm(20) = -2 * t23 * cm(4);
coefm(21) = -t18;
coefm(22) = -t16;
coefm(23) = t32;
coefm(24) = -t27 * cm(2);
coefm(25) = -t27 * cm(7);
coefm(26) = -2 * t23 * cm(9);
coefm(27) = -t27 * cm(9);
coefm(28) = -t27 * cm(10);
coefm(29) = -t27 * cm(8);
coefm(30) = -6 * t14 * cm(10);
coefm(31) = -6 * t14 * cm(8);
coefm(32) = -3 * t23 * cm(10);
coefm(33) = -3 * t23 * cm(8);
coefm(34) = -2 * t23 * cm(2);
coefm(35) = -2 * t23 * cm(7);
coefm(36) = t66;
coefm(37) = -t66;
coefm(38) = t30;
coefm(39) = t34;
coefm(40) = t24;

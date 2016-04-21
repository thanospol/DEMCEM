function [coef,coefm] = coefficients_f1_f1(r1,r2,r3,r4,ko)
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
t5 = bb(3) * dd(3);
t6 = sqrt(0.3e1);
t10 = t1 * t6 / 0.12e2;
t12 = t3 * t6 / 0.12e2;
t13 = bb(1) * dd(1);
t16 = bb(2) * dd(2);
t20 = t2 * t6 / 0.12e2;
t21 = t5 * t6 / 0.6e1 - t10 - t12 + t13 * t6 / 0.6e1 + t16 * t6 / 0.6e1 - t20;
t22 = cc(2) * bb(2);
t25 = cc(1) * bb(1);
t28 = cc(3) * bb(3);
t31 = -t20 - t10 - t12 + t22 * t6 / 0.6e1 + t25 * t6 / 0.6e1 + t28 * t6 / 0.6e1;
t47 = -t25 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + cc(1) * dd(1) / 0.3e1 + t2 / 0.12e2 + cc(3) * dd(3) / 0.3e1 - t28 / 0.6e1 - t16 / 0.6e1 + t3 / 0.12e2 - t13 / 0.6e1 - t22 / 0.6e1 + t1 / 0.12e2 - t5 / 0.6e1;
%
c(1) = t4 / 0.4e1;
c(2) = t4 / 0.4e1;
c(3) = -t4 / 0.4e1;
c(4) = t21;
c(5) = t31;
c(6) = -t31;
c(7) = t21;
c(8) = t4 / 0.4e1;
c(9) = -t31;
c(10) = t47;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * c(4);
t10 = t1 ^ 2;
t12 = t4 ^ 2;
t14 = (0.1e1 / t10 / t12);
t16 = 6 * t14 * c(10);
t19 = (0.1e1 / ko / j);
t23 = (0.1e1 / t1 / t4);
t32 = 2 * t7 * c(7);
t34 = 2 * t7 * c(5);
t36 = 6 * t14 * c(6);
t38 = 2 * t7 * c(8);
t42 = 2 * t7 * c(2);
t45 = 2 * t7 * c(9);
t46 = t23 * c(3);
t47 = t23 * c(1);
%
coef(1) = t9;
coef(2) = -t16;
coef(3) = -t19 * c(4);
coef(4) = -2 * t23 * c(4);
coef(5) = -t19 * c(3);
coef(6) = -t19 * c(9);
coef(7) = -t19 * c(2);
coef(8) = -2 * t23 * c(9);
coef(9) = -t32;
coef(10) = -t34;
coef(11) = -t36;
coef(12) = -t38;
coef(13) = -2 * t23 * c(2);
coef(14) = t32;
coef(15) = t38;
coef(16) = -t9;
coef(17) = -t42;
coef(18) = t36;
coef(19) = t16;
coef(20) = -t19 * c(1);
coef(21) = -t45;
coef(22) = -t46;
coef(23) = t34;
coef(24) = -t47;
coef(25) = t45;
coef(26) = -t19 * c(7);
coef(27) = -t19 * c(8);
coef(28) = -t19 * c(5);
coef(29) = -2 * t23 * c(5);
coef(30) = -t19 * c(6);
coef(31) = -t19 * c(10);
coef(32) = -2 * t23 * c(7);
coef(33) = -2 * t23 * c(8);
coef(34) = -6 * t7 * c(6);
coef(35) = -6 * t7 * c(10);
coef(36) = -3 * t23 * c(6);
coef(37) = -3 * t23 * c(10);
coef(38) = t42;
coef(39) = t47;
coef(40) = t46;
%
t1 = bb(1) ^ 2;
t2 = bb(2) ^ 2;
t3 = bb(3) ^ 2;
t4 = -t1 - t2 - t3;
t5 = cc(1) * bb(1);
t6 = sqrt(0.3e1);
t10 = t3 * t6 / 0.12e2;
t11 = cc(2) * bb(2);
t15 = t2 * t6 / 0.12e2;
t16 = cc(3) * bb(3);
t20 = t1 * t6 / 0.12e2;
t21 = -t5 * t6 / 0.6e1 + t10 - t11 * t6 / 0.6e1 + t15 - t16 * t6 / 0.6e1 + t20;
t22 = bb(1) * dd(1);
t25 = bb(2) * dd(2);
t28 = bb(3) * dd(3);
t31 = t15 - t22 * t6 / 0.6e1 - t25 * t6 / 0.6e1 - t28 * t6 / 0.6e1 + t20 + t10;
t47 = -t5 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + cc(1) * dd(1) / 0.3e1 + t2 / 0.12e2 + cc(3) * dd(3) / 0.3e1 - t16 / 0.6e1 - t25 / 0.6e1 + t3 / 0.12e2 - t22 / 0.6e1 - t11 / 0.6e1 + t1 / 0.12e2 - t28 / 0.6e1;
%
cm(1) = t4 / 0.4e1;
cm(2) = t21;
cm(3) = -t4 / 0.4e1;
cm(4) = -t4 / 0.4e1;
cm(5) = t31;
cm(6) = t4 / 0.4e1;
cm(7) = -t31;
cm(8) = -t21;
cm(9) = t21;
cm(10) = t47;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t12 = (0.1e1 / t1 / ko / t3 / j);
t14 = 2 * t12 * cm(5);
t15 = t1 ^ 2;
t17 = t3 ^ 2;
t19 = (0.1e1 / t15 / t17);
t21 = 6 * t19 * cm(8);
t23 = 2 * t12 * cm(6);
t26 = (0.1e1 / ko / j);
t30 = 6 * t19 * cm(10);
t35 = 2 * t12 * cm(2);
t40 = 2 * t12 * cm(9);
t60 = t5 * cm(1);
t62 = 2 * t12 * cm(7);
t64 = 2 * t12 * cm(4);
t65 = t5 * cm(3);
%
coefm(1) = -2 * t5 * cm(4);
coefm(2) = -t14;
coefm(3) = -t21;
coefm(4) = -t23;
coefm(5) = -t26 * cm(9);
coefm(6) = -t26 * cm(4);
coefm(7) = -t30;
coefm(8) = -t26 * cm(3);
coefm(9) = -2 * t5 * cm(7);
coefm(10) = -t35;
coefm(11) = -t26 * cm(7);
coefm(12) = -2 * t5 * cm(9);
coefm(13) = t40;
coefm(14) = -t26 * cm(8);
coefm(15) = -2 * t5 * cm(5);
coefm(16) = -2 * t5 * cm(6);
coefm(17) = -2 * t5 * cm(2);
coefm(18) = -t26 * cm(2);
coefm(19) = -6 * t12 * cm(10);
coefm(20) = -6 * t12 * cm(8);
coefm(21) = -3 * t5 * cm(10);
coefm(22) = -3 * t5 * cm(8);
coefm(23) = -t26 * cm(10);
coefm(24) = -t26 * cm(5);
coefm(25) = -t26 * cm(6);
coefm(26) = t60;
coefm(27) = -t60;
coefm(28) = t62;
coefm(29) = t64;
coefm(30) = t65;
coefm(31) = -t65;
coefm(32) = -t40;
coefm(33) = -t26 * cm(1);
coefm(34) = t14;
coefm(35) = -t64;
coefm(36) = -t62;
coefm(37) = t23;
coefm(38) = t35;
coefm(39) = t30;
coefm(40) = t21;
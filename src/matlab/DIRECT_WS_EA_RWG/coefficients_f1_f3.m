function [coef,coefm] = coefficients_f1_f3(r1,r2,r3,r4,ko)
%% coefficients_f1_f3

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

c  = zeros(11);
cm = zeros(11);
%
coef  = zeros(43);
coefm = zeros(43);
%
bb = r2-r1;
cc = r3-r1;
dd = r4-r1;
j = sqrt(-1);
%
t1 = bb(1) * dd(1);
t3 = bb(2) * dd(2);
t5 = bb(3) * dd(3);
t7 = bb(2) ^ 2;
t9 = bb(3) ^ 2;
t11 = bb(1) ^ 2;
t14 = t7 + t11 + t9;
t15 = cc(2) * bb(2);
t16 = sqrt(0.3e1);
t18 = t15 * t16 / 0.6e1;
t19 = cc(3) * bb(3);
t21 = t19 * t16 / 0.6e1;
t23 = t9 * t16 / 0.12e2;
t25 = t11 * t16 / 0.12e2;
t26 = cc(1) * bb(1);
t28 = t26 * t16 / 0.6e1;
t30 = t7 * t16 / 0.12e2;
t31 = t18 + t21 - t23 - t25 + t28 - t30;
t33 = t3 * t16 / 0.6e1;
t35 = t1 * t16 / 0.6e1;
t37 = t5 * t16 / 0.6e1;
t38 = -t25 + t33 + t35 - t23 - t30 + t37;
t48 = -t25 - t30 - t23 + t18 - t16 * cc(3) * dd(3) / 0.3e1 + t21 + t33 + t28 + t35 + t37 - t16 * cc(1) * dd(1) / 0.3e1 - t16 * cc(2) * dd(2) / 0.3e1;
t64 = cc(1) * dd(1) / 0.3e1 - t1 / 0.6e1 + t9 / 0.12e2 - t15 / 0.6e1 - t5 / 0.6e1 - t3 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + t7 / 0.12e2 - t19 / 0.6e1 + t11 / 0.12e2 - t26 / 0.6e1 + cc(3) * dd(3) / 0.3e1;
%
c(1) = -t1 / 0.2e1 - t3 / 0.2e1 - t5 / 0.2e1 + t7 / 0.4e1 + t9 / 0.4e1 + t11 / 0.4e1;
c(2) = t14 / 0.4e1;
c(3) = t31;
c(4) = -t14 / 0.4e1;
c(5) = t38;
c(6) = -t14 / 0.4e1;
c(7) = t38;
c(8) = -t31;
c(9) = t48;
c(10) = t64;
c(11) = -t1 / 0.2e1 - t3 / 0.2e1 + t9 / 0.2e1 + t11 / 0.2e1 - t5 / 0.2e1 + t7 / 0.2e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * c(1));
t9 = (0.1e1 / ko / j);
t25 = (0.1e1 / t1 / ko / t3 / j);
t35 = 2 * t25 * c(7);
t37 = 2 * t25 * c(5);
t38 = t5 * c(11);
t40 = 2 * t25 * c(9);
t42 = 2 * t25 * c(4);
t43 = t5 * c(2);
t46 = 2 * t25 * c(3);
t47 = t1 ^ 2;
t49 = t3 ^ 2;
t51 = (0.1e1 / t47 / t49);
t53 = 6 * t51 * c(10);
t55 = 2 * t25 * c(6);
t63 = 6 * t51 * c(8);
%
coef(1) = t6;
coef(2) = -t9 * c(6);
coef(3) = -t9 * c(7);
coef(4) = -2 * t5 * c(3);
coef(5) = -t9 * c(3);
coef(6) = -t9 * c(8);
coef(7) = -t9 * c(10);
coef(8) = -2 * t5 * c(6);
coef(9) = -2 * t5 * c(7);
coef(10) = -6 * t25 * c(8);
coef(11) = -6 * t25 * c(10);
coef(12) = -3 * t5 * c(8);
coef(13) = -3 * t5 * c(10);
coef(14) = t35;
coef(15) = t37;
coef(16) = t38;
coef(17) = t40;
coef(18) = -t6;
coef(19) = t42;
coef(20) = t43;
coef(21) = -t38;
coef(22) = -t9 * c(2);
coef(23) = -t46;
coef(24) = -t53;
coef(25) = -t55;
coef(26) = -t35;
coef(27) = -t9 * c(11);
coef(28) = -2 * t5 * c(9);
coef(29) = -t9 * c(5);
coef(30) = -t9 * c(9);
coef(31) = -t9 * c(4);
coef(32) = -t63;
coef(33) = -2 * t5 * c(5);
coef(34) = -2 * t5 * c(4);
coef(35) = t46;
coef(36) = t53;
coef(37) = t55;
coef(38) = -t42;
coef(39) = -t37;
coef(40) = t63;
coef(41) = -t43;
coef(42) = -t40;
coef(43) = -t9 * c(1);
%
t1 = bb(1) ^ 2;
t3 = bb(2) ^ 2;
t5 = bb(3) ^ 2;
t7 = bb(2) * dd(2);
t9 = bb(3) * dd(3);
t11 = bb(1) * dd(1);
t14 = cc(2) * bb(2);
t16 = cc(1) * bb(1);
t26 = cc(3) * bb(3);
t32 = -t14 / 0.6e1 - t16 / 0.6e1 + t1 / 0.12e2 + t3 / 0.12e2 - t7 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + cc(3) * dd(3) / 0.3e1 + t5 / 0.12e2 - t26 / 0.6e1 + cc(1) * dd(1) / 0.3e1 - t11 / 0.6e1 - t9 / 0.6e1;
t33 = sqrt(0.3e1);
t35 = t14 * t33 / 0.6e1;
t37 = t26 * t33 / 0.6e1;
t39 = t3 * t33 / 0.12e2;
t41 = t16 * t33 / 0.6e1;
t43 = t1 * t33 / 0.12e2;
t45 = t5 * t33 / 0.12e2;
t46 = -t35 - t37 + t39 - t41 + t43 + t45;
t47 = t5 + t1 + t3;
t49 = t7 * t33 / 0.6e1;
t51 = t9 * t33 / 0.6e1;
t53 = t11 * t33 / 0.6e1;
t54 = -t49 + t45 - t51 - t53 + t43 + t39;
t64 = -t45 + t51 - t43 - t39 + t53 - t33 * cc(1) * dd(1) / 0.3e1 - t33 * cc(2) * dd(2) / 0.3e1 - t33 * cc(3) * dd(3) / 0.3e1 + t37 + t41 + t35 + t49;
%
cm(1) = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1 - t7 / 0.2e1 - t9 / 0.2e1 - t11 / 0.2e1;
cm(2) = t32;
cm(3) = t46;
cm(4) = t47 / 0.4e1;
cm(5) = -t47 / 0.4e1;
cm(6) = t54;
cm(7) = -t46;
cm(8) = -t54;
cm(9) = t47 / 0.4e1;
cm(10) = t64;
cm(11) = t9 / 0.2e1 - t5 / 0.2e1 + t11 / 0.2e1 + t7 / 0.2e1 - t1 / 0.2e1 - t3 / 0.2e1;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * cm(10);
t12 = (0.1e1 / t1 / t4);
t13 = t12 * cm(1);
t15 = 2 * t7 * cm(9);
t18 = (0.1e1 / ko / j);
t39 = 2 * t7 * cm(3);
t40 = t1 ^ 2;
t42 = t4 ^ 2;
t44 = (0.1e1 / t40 / t42);
t46 = 6 * t44 * cm(2);
t48 = 6 * t44 * cm(7);
t50 = 2 * t7 * cm(5);
t52 = 2 * t7 * cm(6);
t53 = t12 * cm(4);
t54 = t12 * cm(11);
t57 = 2 * t7 * cm(8);
%
coefm(1) = t9;
coefm(2) = -t13;
coefm(3) = t15;
coefm(4) = -t18 * cm(5);
coefm(5) = -t18 * cm(6);
coefm(6) = -2 * t12 * cm(3);
coefm(7) = -t18 * cm(3);
coefm(8) = -t18 * cm(2);
coefm(9) = -t18 * cm(7);
coefm(10) = -2 * t12 * cm(5);
coefm(11) = -2 * t12 * cm(6);
coefm(12) = -6 * t7 * cm(2);
coefm(13) = -6 * t7 * cm(7);
coefm(14) = -3 * t12 * cm(2);
coefm(15) = -3 * t12 * cm(7);
coefm(16) = t39;
coefm(17) = t46;
coefm(18) = t48;
coefm(19) = t50;
coefm(20) = t52;
coefm(21) = -t9;
coefm(22) = -t53;
coefm(23) = -t54;
coefm(24) = -t18 * cm(1);
coefm(25) = -t57;
coefm(26) = -t15;
coefm(27) = t54;
coefm(28) = t53;
coefm(29) = t57;
coefm(30) = t13;
coefm(31) = -t18 * cm(4);
coefm(32) = -t50;
coefm(33) = -2 * t12 * cm(8);
coefm(34) = -2 * t12 * cm(9);
coefm(35) = -t46;
coefm(36) = -t52;
coefm(37) = -t48;
coefm(38) = -t18 * cm(11);
coefm(39) = -t39;
coefm(40) = -t18 * cm(8);
coefm(41) = -2 * t12 * cm(10);
coefm(42) = -t18 * cm(9);
coefm(43) = -t18 * cm(10);
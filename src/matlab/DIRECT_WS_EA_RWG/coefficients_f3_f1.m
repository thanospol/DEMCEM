function [coef,coefm] = coefficients_f3_f1(r1,r2,r3,r4,ko)
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
t1 = cc(2) * bb(2);
t3 = cc(3) * bb(3);
t5 = bb(3) ^ 2;
t7 = bb(1) ^ 2;
t9 = bb(2) ^ 2;
t11 = cc(1) * bb(1);
t13 = t1 / 0.2e1 + t3 / 0.2e1 - t5 / 0.4e1 - t7 / 0.4e1 - t9 / 0.4e1 + t11 / 0.2e1;
t14 = sqrt(0.3e1);
t16 = t11 * t14 / 0.6e1;
t18 = t3 * t14 / 0.6e1;
t20 = t5 * t14 / 0.12e2;
t22 = t1 * t14 / 0.6e1;
t24 = t9 * t14 / 0.12e2;
t26 = t7 * t14 / 0.12e2;
t27 = -t16 - t18 + t20 - t22 + t24 + t26;
t28 = -t7 - t9 - t5;
t32 = bb(1) * dd(1);
t36 = bb(3) * dd(3);
t46 = bb(2) * dd(2);
t48 = cc(2) * dd(2) / 0.3e1 - t32 / 0.6e1 + cc(3) * dd(3) / 0.3e1 - t36 / 0.6e1 - t3 / 0.6e1 + t5 / 0.12e2 + t7 / 0.12e2 - t11 / 0.6e1 - t1 / 0.6e1 + t9 / 0.12e2 + cc(1) * dd(1) / 0.3e1 - t46 / 0.6e1;
t50 = t32 * t14 / 0.6e1;
t52 = t46 * t14 / 0.6e1;
t57 = t36 * t14 / 0.6e1;
t64 = t16 + t18 + t50 + t52 + t22 - t14 * cc(3) * dd(3) / 0.3e1 - t20 + t57 - t26 - t14 * cc(2) * dd(2) / 0.3e1 - t24 - t14 * cc(1) * dd(1) / 0.3e1;
%
c(1) = t13;
c(2) = t27;
c(3) = t28 / 0.4e1;
c(4) = -t3 / 0.2e1 - t11 / 0.2e1 - t1 / 0.2e1;
c(5) = t27;
c(6) = t48;
c(7) = -t28 / 0.4e1;
c(8) = t64;
c(9) = t13;
c(10) = -t27;
c(11) = t57 - t20 - t26 + t52 - t24 + t50;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * c(1));
t7 = (t5 * c(4));
t12 = (0.1e1 / t1 / ko / t3 / j);
t14 = 2 * t12 * c(9);
t16 = 2 * t12 * c(5);
t17 = t5 * c(7);
t19 = 2 * t12 * c(8);
t20 = t1 ^ 2;
t22 = t3 ^ 2;
t24 = (0.1e1 / t20 / t22);
t26 = 6 * t24 * c(2);
t28 = 6 * t24 * c(6);
t30 = 2 * t12 * c(3);
t32 = 2 * t12 * c(11);
t34 = 2 * t12 * c(10);
t37 = (0.1e1 / ko / j);
%
coef(1) = -t6;
coef(2) = t7;
coef(3) = t14;
coef(4) = t16;
coef(5) = t17;
coef(6) = t19;
coef(7) = t26;
coef(8) = t28;
coef(9) = t30;
coef(10) = t32;
coef(11) = t34;
coef(12) = -t19;
coef(13) = -t14;
coef(14) = -t16;
coef(15) = -t37 * c(1);
coef(16) = -t17;
coef(17) = -t7;
coef(18) = -t37 * c(3);
coef(19) = -t37 * c(11);
coef(20) = -2 * t5 * c(10);
coef(21) = -t37 * c(10);
coef(22) = -t37 * c(2);
coef(23) = -t37 * c(6);
coef(24) = -2 * t5 * c(3);
coef(25) = -2 * t5 * c(11);
coef(26) = -6 * t12 * c(2);
coef(27) = -6 * t12 * c(6);
coef(28) = -3 * t5 * c(2);
coef(29) = -3 * t5 * c(6);
coef(30) = t6;
coef(31) = -t34;
coef(32) = -t26;
coef(33) = -2 * t5 * c(5);
coef(34) = -t37 * c(8);
coef(35) = -t28;
coef(36) = -t37 * c(4);
coef(37) = -t37 * c(5);
coef(38) = -t37 * c(7);
coef(39) = -t32;
coef(40) = -t37 * c(9);
coef(41) = -2 * t5 * c(9);
coef(42) = -t30;
coef(43) = -2 * t5 * c(8);
%
t1 = cc(1) * bb(1);
t3 = cc(2) * bb(2);
t5 = cc(3) * bb(3);
t7 = bb(2) ^ 2;
t9 = bb(3) ^ 2;
t11 = bb(1) ^ 2;
t13 = t1 / 0.2e1 + t3 / 0.2e1 + t5 / 0.2e1 - t7 / 0.4e1 - t9 / 0.4e1 - t11 / 0.4e1;
t15 = -t7 - t9 - t11;
t16 = sqrt(0.3e1);
t18 = t1 * t16 / 0.6e1;
t20 = t9 * t16 / 0.12e2;
t22 = t3 * t16 / 0.6e1;
t24 = t5 * t16 / 0.6e1;
t26 = t7 * t16 / 0.12e2;
t28 = t11 * t16 / 0.12e2;
t29 = t18 - t20 + t22 + t24 - t26 - t28;
t30 = bb(2) * dd(2);
t32 = t30 * t16 / 0.6e1;
t33 = bb(1) * dd(1);
t35 = t33 * t16 / 0.6e1;
t36 = bb(3) * dd(3);
t38 = t36 * t16 / 0.6e1;
t49 = t22 - t28 + t18 + t35 - t26 - t16 * cc(2) * dd(2) / 0.3e1 - t16 * cc(1) * dd(1) / 0.3e1 + t32 - t20 + t38 + t24 - t16 * cc(3) * dd(3) / 0.3e1;
t65 = cc(1) * dd(1) / 0.3e1 - t33 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + t7 / 0.12e2 + cc(3) * dd(3) / 0.3e1 - t30 / 0.6e1 - t1 / 0.6e1 + t11 / 0.12e2 - t3 / 0.6e1 - t36 / 0.6e1 - t5 / 0.6e1 + t9 / 0.12e2;
%
cm(1) = t13;
cm(2) = t5 / 0.2e1 + t3 / 0.2e1 + t1 / 0.2e1;
cm(3) = t15 / 0.4e1;
cm(4) = t29;
cm(5) = -t32 - t35 + t20 + t28 + t26 - t38;
cm(6) = -t15 / 0.4e1;
cm(7) = -t29;
cm(8) = t49;
cm(9) = t65;
cm(10) = -t13;
cm(11) = -t29;
%
t3 = (0.1e1 / ko / j);
t7 = ko ^ 2;
t9 = j ^ 2;
t11 = (0.1e1 / t7 / t9);
t24 = (0.1e1 / t7 / ko / t9 / j);
t33 = (t11 * cm(1));
t35 = 2 * t24 * cm(7);
t36 = t11 * cm(6);
t37 = t11 * cm(2);
t40 = 2 * t24 * cm(3);
t42 = 2 * t24 * cm(5);
t44 = 2 * t24 * cm(11);
t46 = 2 * t24 * cm(8);
t48 = 2 * t24 * cm(10);
t49 = t7 ^ 2;
t51 = t9 ^ 2;
t53 = (0.1e1 / t49 / t51);
t55 = 6 * t53 * cm(9);
t57 = 6 * t53 * cm(4);
%
coefm(1) = -t3 * cm(3);
coefm(2) = -t3 * cm(5);
coefm(3) = -t3 * cm(11);
coefm(4) = -2 * t11 * cm(11);
coefm(5) = -t3 * cm(9);
coefm(6) = -t3 * cm(4);
coefm(7) = -2 * t11 * cm(3);
coefm(8) = -2 * t11 * cm(5);
coefm(9) = -6 * t24 * cm(9);
coefm(10) = -6 * t24 * cm(4);
coefm(11) = -3 * t11 * cm(9);
coefm(12) = -3 * t11 * cm(4);
coefm(13) = t33;
coefm(14) = -t35;
coefm(15) = -t36;
coefm(16) = -t37;
coefm(17) = -t3 * cm(1);
coefm(18) = t40;
coefm(19) = t42;
coefm(20) = t44;
coefm(21) = -t46;
coefm(22) = -t48;
coefm(23) = t55;
coefm(24) = t57;
coefm(25) = t35;
coefm(26) = -t33;
coefm(27) = t36;
coefm(28) = t46;
coefm(29) = t48;
coefm(30) = t37;
coefm(31) = -t57;
coefm(32) = -2 * t11 * cm(10);
coefm(33) = -t42;
coefm(34) = -t3 * cm(8);
coefm(35) = -t40;
coefm(36) = -t44;
coefm(37) = -t3 * cm(7);
coefm(38) = -2 * t11 * cm(7);
coefm(39) = -t3 * cm(2);
coefm(40) = -t3 * cm(10);
coefm(41) = -t55;
coefm(42) = -2 * t11 * cm(8);
coefm(43) = -t3 * cm(6);

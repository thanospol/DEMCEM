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
t5 = bb(1) ^ 2;
t7 = bb(2) ^ 2;
t9 = bb(3) ^ 2;
t11 = cc(1) * bb(1);
t13 = -t1 / 0.2e1 - t3 / 0.2e1 + t5 / 0.4e1 + t7 / 0.4e1 + t9 / 0.4e1 - t11 / 0.2e1;
t14 = sqrt(0.3e1);
t16 = t3 * t14 / 0.6e1;
t18 = t5 * t14 / 0.12e2;
t20 = t1 * t14 / 0.6e1;
t22 = t9 * t14 / 0.12e2;
t24 = t7 * t14 / 0.12e2;
t26 = t11 * t14 / 0.6e1;
t27 = -t16 + t18 - t20 + t22 + t24 - t26;
t28 = -t5 - t7 - t9;
t29 = bb(2) * dd(2);
t31 = t29 * t14 / 0.6e1;
t32 = bb(3) * dd(3);
t34 = t32 * t14 / 0.6e1;
t35 = bb(1) * dd(1);
t37 = t35 * t14 / 0.6e1;
t49 = t34 + t20 + t31 + t26 - t14 * cc(3) * dd(3) / 0.3e1 - t18 + t37 - t24 - t14 * cc(2) * dd(2) / 0.3e1 + t16 - t22 - t14 * cc(1) * dd(1) / 0.3e1;
t65 = t7 / 0.12e2 - t35 / 0.6e1 - t29 / 0.6e1 + t9 / 0.12e2 + t5 / 0.12e2 - t3 / 0.6e1 - t11 / 0.6e1 + cc(2) * dd(2) / 0.3e1 - t1 / 0.6e1 + cc(3) * dd(3) / 0.3e1 + cc(1) * dd(1) / 0.3e1 - t32 / 0.6e1;
%
c(1) = t13;
c(2) = t27;
c(3) = t28 / 0.4e1;
c(4) = -t24 - t22 + t31 + t34 + t37 - t18;
c(5) = -t3 / 0.2e1 - t1 / 0.2e1 + t5 / 0.2e1 - t11 / 0.2e1 + t7 / 0.2e1 + t9 / 0.2e1;
c(6) = -t13;
c(7) = t49;
c(8) = -t27;
c(9) = -t27;
c(10) = t65;
c(11) = -t28 / 0.4e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * c(1));
t7 = (t5 * c(11));
t12 = (0.1e1 / t1 / ko / t3 / j);
t14 = 2 * t12 * c(9);
t15 = t1 ^ 2;
t17 = t3 ^ 2;
t19 = (0.1e1 / t15 / t17);
t21 = 6 * t19 * c(2);
t23 = 6 * t19 * c(10);
t25 = 2 * t12 * c(8);
t27 = 2 * t12 * c(3);
t29 = 2 * t12 * c(4);
t31 = 2 * t12 * c(7);
t33 = 2 * t12 * c(6);
t34 = t5 * c(5);
t37 = (0.1e1 / ko / j);
%
coef(1) = -t6;
coef(2) = t7;
coef(3) = t14;
coef(4) = t21;
coef(5) = t23;
coef(6) = t25;
coef(7) = t27;
coef(8) = t29;
coef(9) = -t31;
coef(10) = -t33;
coef(11) = -t14;
coef(12) = -t7;
coef(13) = -t34;
coef(14) = -t37 * c(1);
coef(15) = -t37 * c(3);
coef(16) = -t37 * c(4);
coef(17) = -2 * t5 * c(8);
coef(18) = -t37 * c(8);
coef(19) = -t37 * c(2);
coef(20) = -t37 * c(10);
coef(21) = -2 * t5 * c(3);
coef(22) = -2 * t5 * c(4);
coef(23) = -6 * t12 * c(2);
coef(24) = -6 * t12 * c(10);
coef(25) = -3 * t5 * c(2);
coef(26) = -3 * t5 * c(10);
coef(27) = -t37 * c(6);
coef(28) = -t37 * c(9);
coef(29) = -2 * t5 * c(9);
coef(30) = -t25;
coef(31) = -t21;
coef(32) = -t23;
coef(33) = -t37 * c(11);
coef(34) = -t27;
coef(35) = -t29;
coef(36) = -2 * t5 * c(6);
coef(37) = -2 * t5 * c(7);
coef(38) = -t37 * c(7);
coef(39) = t6;
coef(40) = -t37 * c(5);
coef(41) = t31;
coef(42) = t33;
coef(43) = t34;
%
t1 = bb(1) ^ 2;
t3 = bb(2) ^ 2;
t5 = bb(3) ^ 2;
t7 = cc(1) * bb(1);
t9 = cc(2) * bb(2);
t11 = cc(3) * bb(3);
t13 = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1 - t7 / 0.2e1 - t9 / 0.2e1 - t11 / 0.2e1;
t15 = bb(1) * dd(1);
t16 = sqrt(0.3e1);
t18 = t15 * t16 / 0.6e1;
t20 = t1 * t16 / 0.12e2;
t25 = t3 * t16 / 0.12e2;
t27 = t9 * t16 / 0.6e1;
t32 = t5 * t16 / 0.12e2;
t34 = t11 * t16 / 0.6e1;
t38 = bb(3) * dd(3);
t40 = t38 * t16 / 0.6e1;
t41 = bb(2) * dd(2);
t43 = t41 * t16 / 0.6e1;
t45 = t7 * t16 / 0.6e1;
t46 = t18 - t20 - t16 * cc(1) * dd(1) / 0.3e1 - t25 + t27 - t16 * cc(3) * dd(3) / 0.3e1 - t32 + t34 - t16 * cc(2) * dd(2) / 0.3e1 + t40 + t43 + t45;
t47 = -t1 - t3 - t5;
t48 = -t25 - t20 + t34 + t45 - t32 + t27;
t65 = t5 / 0.12e2 - t15 / 0.6e1 + cc(2) * dd(2) / 0.3e1 - t11 / 0.6e1 - t41 / 0.6e1 - t9 / 0.6e1 - t38 / 0.6e1 + t1 / 0.12e2 - t7 / 0.6e1 + cc(1) * dd(1) / 0.3e1 + t3 / 0.12e2 + cc(3) * dd(3) / 0.3e1;
%
cm(1) = t13;
cm(2) = -t1 / 0.2e1 - t5 / 0.2e1 + t11 / 0.2e1 - t3 / 0.2e1 + t9 / 0.2e1 + t7 / 0.2e1;
cm(3) = t13;
cm(4) = t46;
cm(5) = t47 / 0.4e1;
cm(6) = t48;
cm(7) = t25 + t20 + t32 - t43 - t40 - t18;
cm(8) = t65;
cm(9) = t48;
cm(10) = -t48;
cm(11) = -t47 / 0.4e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * cm(2));
t11 = (0.1e1 / t1 / ko / t3 / j);
t13 = 2 * t11 * cm(4);
t15 = 2 * t11 * cm(3);
t16 = t1 ^ 2;
t18 = t3 ^ 2;
t20 = (0.1e1 / t16 / t18);
t22 = 6 * t20 * cm(8);
t24 = 6 * t20 * cm(6);
t26 = 2 * t11 * cm(10);
t28 = 2 * t11 * cm(5);
t30 = 2 * t11 * cm(7);
t32 = 2 * t11 * cm(9);
t33 = t5 * cm(11);
t36 = (0.1e1 / ko / j);
t38 = t5 * cm(1);
%
coefm(1) = t6;
coefm(2) = -t13;
coefm(3) = -t15;
coefm(4) = t22;
coefm(5) = t24;
coefm(6) = t26;
coefm(7) = t28;
coefm(8) = t30;
coefm(9) = -t32;
coefm(10) = -t33;
coefm(11) = -t6;
coefm(12) = -t36 * cm(1);
coefm(13) = t13;
coefm(14) = -t38;
coefm(15) = t15;
coefm(16) = t33;
coefm(17) = -t36 * cm(11);
coefm(18) = -t36 * cm(9);
coefm(19) = -t24;
coefm(20) = -t22;
coefm(21) = -t36 * cm(2);
coefm(22) = -t36 * cm(3);
coefm(23) = -t26;
coefm(24) = -t36 * cm(4);
coefm(25) = -2 * t5 * cm(4);
coefm(26) = -t28;
coefm(27) = -t30;
coefm(28) = -2 * t5 * cm(3);
coefm(29) = -2 * t5 * cm(9);
coefm(30) = t38;
coefm(31) = t32;
coefm(32) = -t36 * cm(6);
coefm(33) = -2 * t5 * cm(5);
coefm(34) = -2 * t5 * cm(7);
coefm(35) = -6 * t11 * cm(8);
coefm(36) = -6 * t11 * cm(6);
coefm(37) = -3 * t5 * cm(8);
coefm(38) = -3 * t5 * cm(6);
coefm(39) = -t36 * cm(7);
coefm(40) = -t36 * cm(8);
coefm(41) = -t36 * cm(10);
coefm(42) = -2 * t5 * cm(10);
coefm(43) = -t36 * cm(5);

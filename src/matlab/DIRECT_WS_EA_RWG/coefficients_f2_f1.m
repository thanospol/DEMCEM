function [coef,coefm] = coefficients_f2_f1(r1,r2,r3,r4,ko)
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
t1 = bb(3) ^ 2;
t2 = bb(1) ^ 2;
t3 = bb(2) ^ 2;
t4 = t1 + t2 + t3;
t5 = sqrt(0.3e1);
t7 = t2 * t5 / 0.12e2;
t9 = t3 * t5 / 0.12e2;
t10 = cc(1) * bb(1);
t14 = t1 * t5 / 0.12e2;
t15 = cc(3) * bb(3);
t18 = cc(2) * bb(2);
t21 = t7 + t9 - t10 * t5 / 0.6e1 + t14 - t15 * t5 / 0.6e1 - t18 * t5 / 0.6e1;
t22 = bb(2) * dd(2);
t25 = bb(3) * dd(3);
t28 = bb(1) * dd(1);
t31 = t22 * t5 / 0.6e1 + t25 * t5 / 0.6e1 - t9 + t28 * t5 / 0.6e1 - t7 - t14;
t47 = -t18 / 0.6e1 + t3 / 0.12e2 - t10 / 0.6e1 - t28 / 0.6e1 + cc(1) * dd(1) / 0.3e1 + t2 / 0.12e2 - t15 / 0.6e1 - t25 / 0.6e1 - t22 / 0.6e1 + t1 / 0.12e2 + cc(2) * dd(2) / 0.3e1 + cc(3) * dd(3) / 0.3e1;
%
c(1) = t4 / 0.4e1;
c(2) = t21;
c(3) = t31;
c(4) = -t4 / 0.4e1;
c(5) = -t4 / 0.2e1;
c(6) = t21;
c(7) = -t31;
c(8) = t4 / 0.4e1;
c(9) = -t21;
c(10) = t47;
c(11) = t4 / 0.4e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * c(11));
t11 = (0.1e1 / t1 / ko / t3 / j);
t13 = 2 * t11 * c(8);
t14 = t5 * c(5);
t16 = 2 * t11 * c(7);
t17 = t5 * c(1);
t18 = t1 ^ 2;
t20 = t3 ^ 2;
t22 = (0.1e1 / t18 / t20);
t24 = 6 * t22 * c(2);
t26 = 6 * t22 * c(10);
t28 = 2 * t11 * c(9);
t30 = 2 * t11 * c(4);
t32 = 2 * t11 * c(3);
t34 = 2 * t11 * c(6);
t37 = (0.1e1 / ko / j);
%
coef(1) = t6;
coef(2) = t13;
coef(3) = t14;
coef(4) = t16;
coef(5) = -t17;
coef(6) = t24;
coef(7) = t26;
coef(8) = t28;
coef(9) = t30;
coef(10) = t32;
coef(11) = -t16;
coef(12) = -t13;
coef(13) = -t34;
coef(14) = -t6;
coef(15) = -t14;
coef(16) = -t37 * c(1);
coef(17) = -t37 * c(4);
coef(18) = -t37 * c(3);
coef(19) = -2 * t5 * c(9);
coef(20) = -t37 * c(9);
coef(21) = -t37 * c(2);
coef(22) = -t37 * c(10);
coef(23) = -2 * t5 * c(4);
coef(24) = -2 * t5 * c(3);
coef(25) = -6 * t11 * c(2);
coef(26) = -6 * t11 * c(10);
coef(27) = -3 * t5 * c(2);
coef(28) = -3 * t5 * c(10);
coef(29) = -t37 * c(8);
coef(30) = -t30;
coef(31) = -2 * t5 * c(7);
coef(32) = -t37 * c(7);
coef(33) = -t32;
coef(34) = -2 * t5 * c(8);
coef(35) = -t28;
coef(36) = -t37 * c(11);
coef(37) = -t37 * c(5);
coef(38) = -t37 * c(6);
coef(39) = -2 * t5 * c(6);
coef(40) = -t24;
coef(41) = -t26;
coef(42) = t17;
coef(43) = t34;
%
t1 = bb(1) ^ 2;
t2 = bb(2) ^ 2;
t3 = bb(3) ^ 2;
t4 = t1 + t2 + t3;
t5 = cc(3) * bb(3);
t9 = cc(2) * bb(2);
t13 = cc(1) * bb(1);
t16 = bb(2) * dd(2);
t18 = bb(1) * dd(1);
t22 = bb(3) * dd(3);
t26 = -t5 / 0.6e1 + cc(3) * dd(3) / 0.3e1 - t9 / 0.6e1 + t2 / 0.12e2 + t3 / 0.12e2 - t13 / 0.6e1 + t1 / 0.12e2 - t16 / 0.6e1 - t18 / 0.6e1 + cc(1) * dd(1) / 0.3e1 - t22 / 0.6e1 + cc(2) * dd(2) / 0.3e1;
t27 = sqrt(0.3e1);
t29 = t1 * t27 / 0.12e2;
t31 = t3 * t27 / 0.12e2;
t37 = t2 * t27 / 0.12e2;
t40 = t29 + t31 - t18 * t27 / 0.6e1 - t16 * t27 / 0.6e1 + t37 - t22 * t27 / 0.6e1;
t47 = t13 * t27 / 0.6e1 + t5 * t27 / 0.6e1 - t29 - t31 + t9 * t27 / 0.6e1 - t37;
%
cm(1) = t4 / 0.4e1;
cm(2) = t26;
cm(3) = t4 / 0.4e1;
cm(4) = t4 / 0.2e1;
cm(5) = t40;
cm(6) = t47;
cm(7) = t40;
cm(8) = -t4 / 0.4e1;
cm(9) = -t47;
cm(10) = -t47;
cm(11) = -t4 / 0.4e1;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * cm(5);
t12 = (0.1e1 / t1 / t4);
t13 = t12 * cm(3);
t14 = t12 * cm(4);
t15 = t1 ^ 2;
t17 = t4 ^ 2;
t19 = (0.1e1 / t15 / t17);
t21 = 6 * t19 * cm(2);
t23 = 6 * t19 * cm(6);
t26 = (0.1e1 / ko / j);
t29 = 2 * t7 * cm(9);
t31 = 2 * t7 * cm(10);
t33 = 2 * t7 * cm(7);
t35 = 2 * t7 * cm(8);
t37 = 2 * t7 * cm(11);
t38 = t12 * cm(1);
%
coefm(1) = t9;
coefm(2) = t13;
coefm(3) = t14;
coefm(4) = t21;
coefm(5) = t23;
coefm(6) = -t14;
coefm(7) = -t13;
coefm(8) = -t26 * cm(1);
coefm(9) = -t29;
coefm(10) = t31;
coefm(11) = t33;
coefm(12) = t35;
coefm(13) = -t37;
coefm(14) = -t9;
coefm(15) = t38;
coefm(16) = -t31;
coefm(17) = -2 * t12 * cm(9);
coefm(18) = -t26 * cm(9);
coefm(19) = -t26 * cm(4);
coefm(20) = -t33;
coefm(21) = -t35;
coefm(22) = -2 * t12 * cm(5);
coefm(23) = -2 * t12 * cm(11);
coefm(24) = -t23;
coefm(25) = -t26 * cm(3);
coefm(26) = -t26 * cm(5);
coefm(27) = -t26 * cm(11);
coefm(28) = -t21;
coefm(29) = -t38;
coefm(30) = t37;
coefm(31) = -t26 * cm(8);
coefm(32) = -t26 * cm(7);
coefm(33) = -2 * t12 * cm(10);
coefm(34) = -t26 * cm(10);
coefm(35) = -t26 * cm(2);
coefm(36) = -t26 * cm(6);
coefm(37) = -2 * t12 * cm(8);
coefm(38) = -2 * t12 * cm(7);
coefm(39) = -6 * t7 * cm(2);
coefm(40) = -6 * t7 * cm(6);
coefm(41) = -3 * t12 * cm(2);
coefm(42) = -3 * t12 * cm(6);
coefm(43) = t29;

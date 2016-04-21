function [coef,coefm] = coefficients_f3_f3(r1,r2,r3,r4,ko)
%% coefficients_f3_f3

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
t3 = bb(3) * dd(3);
t5 = bb(3) ^ 2;
t6 = t5 / 0.4e1;
t7 = bb(2) ^ 2;
t8 = t7 / 0.4e1;
t9 = bb(1) ^ 2;
t10 = t9 / 0.4e1;
t11 = cc(1) * dd(1);
t12 = cc(2) * bb(2);
t13 = t12 / 0.2e1;
t14 = cc(3) * bb(3);
t15 = t14 / 0.2e1;
t16 = bb(2) * dd(2);
t18 = cc(2) * dd(2);
t19 = cc(3) * dd(3);
t20 = cc(1) * bb(1);
t21 = t20 / 0.2e1;
t22 = -t1 / 0.2e1 - t3 / 0.2e1 + t6 + t8 + t10 + t11 - t13 - t15 - t16 / 0.2e1 + t18 + t19 - t21;
t25 = sqrt(0.3e1);
t27 = t7 * t25 / 0.12e2;
t29 = t1 * t25 / 0.6e1;
t31 = t9 * t25 / 0.12e2;
t33 = t16 * t25 / 0.6e1;
t35 = t3 * t25 / 0.6e1;
t37 = t5 * t25 / 0.12e2;
t39 = -t9 - t5 - t7;
t41 = t14 * t25 / 0.6e1;
t43 = t12 * t25 / 0.6e1;
t45 = t20 * t25 / 0.6e1;
t46 = t27 - t41 - t43 + t37 - t45 + t31;
t56 = t41 - t31 - t37 + t45 + t43 - t25 * cc(2) * dd(2) / 0.3e1 + t35 - t25 * cc(1) * dd(1) / 0.3e1 + t29 - t27 + t33 - t25 * cc(3) * dd(3) / 0.3e1;
t69 = t7 / 0.12e2 - t1 / 0.6e1 - t3 / 0.6e1 + t11 / 0.3e1 - t20 / 0.6e1 + t18 / 0.3e1 + t19 / 0.3e1 + t5 / 0.12e2 - t14 / 0.6e1 + t9 / 0.12e2 - t12 / 0.6e1 - t16 / 0.6e1;
%
c(1) = t22;
c(2) = t21 + t13 - t10 + t15 - t6 - t8;
c(3) = -t14 / 0.2e1 - t20 / 0.2e1 + t9 / 0.2e1 - t16 / 0.2e1 - t12 / 0.2e1 - t1 / 0.2e1 - t3 / 0.2e1 + t7 / 0.2e1 + t5 / 0.2e1;
c(4) = -t27 + t29 - t31 + t33 + t35 - t37;
c(5) = t39 / 0.4e1;
c(6) = t46;
c(7) = t56;
c(8) = t69;
c(9) = t56;
c(10) = -t46;
c(11) = -t39 / 0.4e1;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * c(7);
t12 = (0.1e1 / t1 / t4);
t13 = t12 * c(1);
t14 = t12 * c(11);
t15 = t12 * c(3);
t17 = 2 * t7 * c(9);
t19 = 2 * t7 * c(2);
t22 = (0.1e1 / ko / j);
t42 = t1 ^ 2;
t44 = t4 ^ 2;
t46 = (0.1e1 / t42 / t44);
t48 = 6 * t46 * c(6);
t50 = 6 * t46 * c(8);
t52 = 2 * t7 * c(10);
t54 = 2 * t7 * c(5);
t56 = 2 * t7 * c(4);
%
coef(1) = t9;
coef(2) = -t13;
coef(3) = t14;
coef(4) = t15;
coef(5) = t17;
coef(6) = t19;
coef(7) = -t22 * c(5);
coef(8) = -t22 * c(4);
coef(9) = -t22 * c(10);
coef(10) = -2 * t12 * c(10);
coef(11) = -t22 * c(6);
coef(12) = -t22 * c(8);
coef(13) = -2 * t12 * c(5);
coef(14) = -2 * t12 * c(4);
coef(15) = -6 * t7 * c(6);
coef(16) = -6 * t7 * c(8);
coef(17) = -3 * t12 * c(6);
coef(18) = -3 * t12 * c(8);
coef(19) = t48;
coef(20) = t50;
coef(21) = t52;
coef(22) = t54;
coef(23) = t56;
coef(24) = -t17;
coef(25) = -t19;
coef(26) = -t9;
coef(27) = -t14;
coef(28) = -t15;
coef(29) = -t22 * c(1);
coef(30) = t13;
coef(31) = -t22 * c(2);
coef(32) = -t22 * c(11);
coef(33) = -t52;
coef(34) = -2 * t12 * c(2);
coef(35) = -t22 * c(7);
coef(36) = -t22 * c(9);
coef(37) = -t22 * c(3);
coef(38) = -2 * t12 * c(9);
coef(39) = -t56;
coef(40) = -t54;
coef(41) = -2 * t12 * c(7);
coef(42) = -t48;
coef(43) = -t50;
%
t1 = cc(1) * bb(1);
t2 = t1 / 0.2e1;
t3 = bb(2) ^ 2;
t4 = t3 / 0.4e1;
t5 = bb(1) ^ 2;
t6 = t5 / 0.4e1;
t7 = bb(3) ^ 2;
t8 = t7 / 0.4e1;
t9 = cc(2) * bb(2);
t10 = t9 / 0.2e1;
t11 = bb(2) * dd(2);
t13 = cc(3) * bb(3);
t14 = t13 / 0.2e1;
t15 = bb(3) * dd(3);
t17 = cc(1) * dd(1);
t18 = cc(2) * dd(2);
t19 = cc(3) * dd(3);
t20 = bb(1) * dd(1);
t22 = -t2 + t4 + t6 + t8 - t10 - t11 / 0.2e1 - t14 - t15 / 0.2e1 + t17 + t18 + t19 - t20 / 0.2e1;
t23 = sqrt(0.3e1);
t25 = t3 * t23 / 0.12e2;
t27 = t13 * t23 / 0.6e1;
t29 = t9 * t23 / 0.6e1;
t31 = t1 * t23 / 0.6e1;
t33 = t5 * t23 / 0.12e2;
t35 = t7 * t23 / 0.12e2;
t36 = -t25 + t27 + t29 + t31 - t33 - t35;
t38 = t20 * t23 / 0.6e1;
t40 = t15 * t23 / 0.6e1;
t42 = t11 * t23 / 0.6e1;
t44 = -t3 - t7 - t5;
t57 = t7 / 0.12e2 + t18 / 0.3e1 + t19 / 0.3e1 - t20 / 0.6e1 + t17 / 0.3e1 - t11 / 0.6e1 + t5 / 0.12e2 - t1 / 0.6e1 - t15 / 0.6e1 - t9 / 0.6e1 - t13 / 0.6e1 + t3 / 0.12e2;
t67 = -t23 * cc(1) * dd(1) / 0.3e1 + t29 - t35 - t23 * cc(3) * dd(3) / 0.3e1 + t38 - t23 * cc(2) * dd(2) / 0.3e1 - t33 + t42 - t25 + t27 + t31 + t40;
%
cm(1) = t22;
cm(2) = t36;
cm(3) = t25 - t38 + t35 + t33 - t40 - t42;
cm(4) = t44 / 0.4e1;
cm(5) = t57;
cm(6) = t67;
cm(7) = t11 / 0.2e1 - t5 / 0.2e1 - t7 / 0.2e1 - t3 / 0.2e1 + t13 / 0.2e1 + t9 / 0.2e1 + t1 / 0.2e1 + t15 / 0.2e1 + t20 / 0.2e1;
cm(8) = -t44 / 0.4e1;
cm(9) = -t36;
cm(10) = t67;
cm(11) = t6 + t8 + t4 - t14 - t2 - t10;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t10 = (0.1e1 / ko / j);
t18 = (0.1e1 / t1 / ko / t3 / j);
t20 = 2 * t18 * cm(4);
t22 = 2 * t18 * cm(3);
t23 = t1 ^ 2;
t25 = t3 ^ 2;
t27 = (0.1e1 / t23 / t25);
t29 = 6 * t27 * cm(2);
t34 = 6 * t27 * cm(5);
t36 = 2 * t18 * cm(9);
t40 = t5 * cm(7);
t41 = t5 * cm(8);
t44 = 2 * t18 * cm(6);
t46 = 2 * t18 * cm(10);
t48 = 2 * t18 * cm(11);
t68 = t5 * cm(1);
%
coefm(1) = -2 * t5 * cm(6);
coefm(2) = -t10 * cm(8);
coefm(3) = -2 * t5 * cm(10);
coefm(4) = -t20;
coefm(5) = -t22;
coefm(6) = -t29;
coefm(7) = -2 * t5 * cm(11);
coefm(8) = -t10 * cm(7);
coefm(9) = -t34;
coefm(10) = -t36;
coefm(11) = -t10 * cm(6);
coefm(12) = -t10 * cm(11);
coefm(13) = -t10 * cm(10);
coefm(14) = t40;
coefm(15) = t20;
coefm(16) = t22;
coefm(17) = -t41;
coefm(18) = t34;
coefm(19) = t29;
coefm(20) = -t40;
coefm(21) = -t10 * cm(1);
coefm(22) = -t44;
coefm(23) = t36;
coefm(24) = -t46;
coefm(25) = -t48;
coefm(26) = -t10 * cm(4);
coefm(27) = -t10 * cm(3);
coefm(28) = -t10 * cm(9);
coefm(29) = -2 * t5 * cm(9);
coefm(30) = -t10 * cm(5);
coefm(31) = -t10 * cm(2);
coefm(32) = -2 * t5 * cm(4);
coefm(33) = -2 * t5 * cm(3);
coefm(34) = -6 * t18 * cm(5);
coefm(35) = -6 * t18 * cm(2);
coefm(36) = -3 * t5 * cm(5);
coefm(37) = -3 * t5 * cm(2);
coefm(38) = t48;
coefm(39) = t68;
coefm(40) = -t68;
coefm(41) = t44;
coefm(42) = t41;
coefm(43) = t46;

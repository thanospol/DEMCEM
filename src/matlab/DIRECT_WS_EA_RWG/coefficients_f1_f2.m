function [coef,coefm] = coefficients_f1_f2(r1,r2,r3,r4,ko)
%% coefficients_f1_f2

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
t1 = bb(2) ^ 2;
t2 = bb(3) ^ 2;
t3 = bb(1) ^ 2;
t4 = t1 + t2 + t3;
t5 = bb(3) * dd(3);
t6 = sqrt(0.3e1);
t9 = bb(1) * dd(1);
t12 = bb(2) * dd(2);
t16 = t3 * t6 / 0.12e2;
t18 = t1 * t6 / 0.12e2;
t20 = t2 * t6 / 0.12e2;
t21 = t5 * t6 / 0.6e1 + t9 * t6 / 0.6e1 + t12 * t6 / 0.6e1 - t16 - t18 - t20;
t22 = cc(2) * bb(2);
t25 = cc(1) * bb(1);
t28 = cc(3) * bb(3);
t31 = -t22 * t6 / 0.6e1 + t20 - t25 * t6 / 0.6e1 + t18 - t28 * t6 / 0.6e1 + t16;
t47 = -t25 / 0.6e1 - t9 / 0.6e1 - t28 / 0.6e1 + cc(1) * dd(1) / 0.3e1 - t5 / 0.6e1 + t2 / 0.12e2 + t3 / 0.12e2 - t22 / 0.6e1 + cc(3) * dd(3) / 0.3e1 + t1 / 0.12e2 - t12 / 0.6e1 + cc(2) * dd(2) / 0.3e1;
%
c(1) = t4 / 0.4e1;
c(2) = t21;
c(3) = t31;
c(4) = -t4 / 0.4e1;
c(5) = -t4 / 0.4e1;
c(6) = t21;
c(7) = t4 / 0.2e1;
c(8) = -t31;
c(9) = t47;
c(10) = t4 / 0.4e1;
c(11) = -t31;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * c(7));
t11 = (0.1e1 / t1 / ko / t3 / j);
t13 = 2 * t11 * c(8);
t15 = 2 * t11 * c(6);
t17 = 2 * t11 * c(5);
t18 = t5 * c(10);
t19 = t1 ^ 2;
t21 = t3 ^ 2;
t23 = (0.1e1 / t19 / t21);
t25 = 6 * t23 * c(3);
t27 = 6 * t23 * c(9);
t29 = 2 * t11 * c(4);
t31 = 2 * t11 * c(2);
t33 = 2 * t11 * c(11);
t36 = (0.1e1 / ko / j);
t38 = t5 * c(1);
%
coef(1) = t6;
coef(2) = t13;
coef(3) = t15;
coef(4) = t17;
coef(5) = t18;
coef(6) = t25;
coef(7) = t27;
coef(8) = t29;
coef(9) = t31;
coef(10) = t33;
coef(11) = -t15;
coef(12) = -t17;
coef(13) = -t6;
coef(14) = -t36 * c(1);
coef(15) = -t13;
coef(16) = -t18;
coef(17) = -t38;
coef(18) = -t36 * c(4);
coef(19) = -t36 * c(2);
coef(20) = -t36 * c(11);
coef(21) = -2 * t5 * c(11);
coef(22) = -t36 * c(3);
coef(23) = -t36 * c(9);
coef(24) = -2 * t5 * c(4);
coef(25) = -2 * t5 * c(2);
coef(26) = -6 * t11 * c(3);
coef(27) = -6 * t11 * c(9);
coef(28) = -3 * t5 * c(3);
coef(29) = -3 * t5 * c(9);
coef(30) = t38;
coef(31) = -t36 * c(6);
coef(32) = -t33;
coef(33) = -t31;
coef(34) = -t36 * c(7);
coef(35) = -t25;
coef(36) = -t29;
coef(37) = -2 * t5 * c(5);
coef(38) = -t36 * c(8);
coef(39) = -2 * t5 * c(6);
coef(40) = -2 * t5 * c(8);
coef(41) = -t36 * c(10);
coef(42) = -t27;
coef(43) = -t36 * c(5);
%
t1 = bb(1) ^ 2;
t2 = bb(2) ^ 2;
t3 = bb(3) ^ 2;
t4 = t1 + t2 + t3;
t5 = bb(2) * dd(2);
t6 = sqrt(0.3e1);
t10 = t2 * t6 / 0.12e2;
t12 = t1 * t6 / 0.12e2;
t14 = t3 * t6 / 0.12e2;
t15 = bb(1) * dd(1);
t18 = bb(3) * dd(3);
t21 = -t5 * t6 / 0.6e1 + t10 + t12 + t14 - t15 * t6 / 0.6e1 - t18 * t6 / 0.6e1;
t22 = cc(1) * bb(1);
t25 = cc(2) * bb(2);
t28 = cc(3) * bb(3);
t31 = t22 * t6 / 0.6e1 - t14 + t25 * t6 / 0.6e1 + t28 * t6 / 0.6e1 - t10 - t12;
t47 = -t25 / 0.6e1 - t28 / 0.6e1 - t5 / 0.6e1 + cc(1) * dd(1) / 0.3e1 - t22 / 0.6e1 + cc(2) * dd(2) / 0.3e1 - t18 / 0.6e1 + cc(3) * dd(3) / 0.3e1 + t1 / 0.12e2 - t15 / 0.6e1 + t3 / 0.12e2 + t2 / 0.12e2;
%
cm(1) = t4 / 0.4e1;
cm(2) = -t4 / 0.2e1;
cm(3) = -t4 / 0.4e1;
cm(4) = t21;
cm(5) = t31;
cm(6) = -t31;
cm(7) = t31;
cm(8) = t4 / 0.4e1;
cm(9) = t47;
cm(10) = -t21;
cm(11) = t4 / 0.4e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t6 = (t5 * cm(1));
t7 = (t5 * cm(2));
t10 = (0.1e1 / ko / j);
t16 = (0.1e1 / t1 / ko / t3 / j);
t18 = 2 * t16 * cm(6);
t20 = 2 * t16 * cm(3);
t22 = 2 * t16 * cm(4);
t24 = 2 * t16 * cm(10);
t25 = t1 ^ 2;
t27 = t3 ^ 2;
t29 = (0.1e1 / t25 / t27);
t31 = 6 * t29 * cm(5);
t33 = 6 * t29 * cm(9);
t35 = 2 * t16 * cm(11);
t37 = 2 * t16 * cm(7);
t38 = t5 * cm(8);
%
coefm(1) = t6;
coefm(2) = -t7;
coefm(3) = -t10 * cm(1);
coefm(4) = t18;
coefm(5) = t20;
coefm(6) = t22;
coefm(7) = -t24;
coefm(8) = t31;
coefm(9) = t33;
coefm(10) = -t35;
coefm(11) = -t37;
coefm(12) = -t38;
coefm(13) = t7;
coefm(14) = t24;
coefm(15) = t35;
coefm(16) = -t10 * cm(3);
coefm(17) = -t10 * cm(4);
coefm(18) = -2 * t5 * cm(6);
coefm(19) = -t10 * cm(6);
coefm(20) = -t10 * cm(9);
coefm(21) = -t10 * cm(5);
coefm(22) = -2 * t5 * cm(3);
coefm(23) = -2 * t5 * cm(4);
coefm(24) = -6 * t16 * cm(9);
coefm(25) = -6 * t16 * cm(5);
coefm(26) = -3 * t5 * cm(9);
coefm(27) = -3 * t5 * cm(5);
coefm(28) = -t6;
coefm(29) = t38;
coefm(30) = t37;
coefm(31) = -t10 * cm(11);
coefm(32) = -t10 * cm(2);
coefm(33) = -2 * t5 * cm(11);
coefm(34) = -t10 * cm(8);
coefm(35) = -t22;
coefm(36) = -2 * t5 * cm(7);
coefm(37) = -t18;
coefm(38) = -2 * t5 * cm(10);
coefm(39) = -t10 * cm(10);
coefm(40) = -t10 * cm(7);
coefm(41) = -t20;
coefm(42) = -t33;
coefm(43) = -t31;

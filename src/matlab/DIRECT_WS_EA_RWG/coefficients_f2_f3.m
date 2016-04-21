function [coef,coefm] = coefficients_f2_f3(r1,r2,r3,r4,ko)
%% coefficients_f2_f3

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
t1 = bb(1) ^ 2;
t3 = bb(2) ^ 2;
t5 = bb(3) ^ 2;
t7 = bb(3) * dd(3);
t9 = bb(1) * dd(1);
t11 = bb(2) * dd(2);
t14 = -t1 - t3 - t5;
t15 = sqrt(0.3e1);
t17 = t9 * t15 / 0.6e1;
t19 = t11 * t15 / 0.6e1;
t21 = t1 * t15 / 0.12e2;
t23 = t3 * t15 / 0.12e2;
t25 = t5 * t15 / 0.12e2;
t27 = t7 * t15 / 0.6e1;
t28 = t17 + t19 - t21 - t23 - t25 + t27;
t29 = cc(2) * bb(2);
t31 = t29 * t15 / 0.6e1;
t32 = cc(1) * bb(1);
t34 = t32 * t15 / 0.6e1;
t35 = cc(3) * bb(3);
t37 = t35 * t15 / 0.6e1;
t38 = -t23 + t31 + t34 - t25 - t21 + t37;
t49 = -t21 - t25 - t15 * cc(1) * dd(1) / 0.3e1 - t15 * cc(3) * dd(3) / 0.3e1 - t23 + t37 + t34 - t15 * cc(2) * dd(2) / 0.3e1 + t31 + t17 + t19 + t27;
t65 = t5 / 0.12e2 + cc(3) * dd(3) / 0.3e1 - t7 / 0.6e1 + t3 / 0.12e2 - t9 / 0.6e1 - t32 / 0.6e1 - t11 / 0.6e1 + cc(2) * dd(2) / 0.3e1 - t29 / 0.6e1 + t1 / 0.12e2 - t35 / 0.6e1 + cc(1) * dd(1) / 0.3e1;
%
c(1) = -t1 / 0.4e1 - t3 / 0.4e1 - t5 / 0.4e1 + t7 / 0.2e1 + t9 / 0.2e1 + t11 / 0.2e1;
c(2) = t14 / 0.4e1;
c(3) = t28;
c(4) = -t14 / 0.4e1;
c(5) = -t14 / 0.4e1;
c(6) = -t28;
c(7) = t38;
c(8) = -t11 / 0.2e1 - t7 / 0.2e1 - t9 / 0.2e1;
c(9) = -t38;
c(10) = t49;
c(11) = t65;
%
t1 = ko ^ 2;
t4 = j ^ 2;
t7 = (0.1e1 / t1 / ko / t4 / j);
t9 = 2 * t7 * c(2);
t12 = (0.1e1 / t1 / t4);
t15 = t1 ^ 2;
t17 = t4 ^ 2;
t19 = (0.1e1 / t15 / t17);
t21 = 6 * t19 * c(11);
t25 = 2 * t7 * c(3);
t28 = (0.1e1 / ko / j);
t31 = 6 * t19 * c(9);
t34 = 2 * t7 * c(7);
t40 = t12 * c(1);
t41 = t12 * c(4);
t43 = 2 * t7 * c(5);
t45 = 2 * t7 * c(6);
t47 = 2 * t7 * c(10);
t48 = t12 * c(8);
%
coef(1) = -t9;
coef(2) = -2 * t12 * c(6);
coef(3) = -t21;
coef(4) = -2 * t12 * c(5);
coef(5) = -t25;
coef(6) = -t28 * c(6);
coef(7) = -t31;
coef(8) = -t28 * c(5);
coef(9) = -t34;
coef(10) = -t28 * c(4);
coef(11) = -t28 * c(8);
coef(12) = -t28 * c(10);
coef(13) = -2 * t12 * c(10);
coef(14) = t40;
coef(15) = t41;
coef(16) = t43;
coef(17) = t45;
coef(18) = -t40;
coef(19) = t47;
coef(20) = t48;
coef(21) = -t28 * c(2);
coef(22) = -t28 * c(3);
coef(23) = -t28 * c(7);
coef(24) = -2 * t12 * c(7);
coef(25) = -t28 * c(9);
coef(26) = -t28 * c(11);
coef(27) = -2 * t12 * c(2);
coef(28) = -2 * t12 * c(3);
coef(29) = -6 * t7 * c(9);
coef(30) = -6 * t7 * c(11);
coef(31) = -3 * t12 * c(9);
coef(32) = -3 * t12 * c(11);
coef(33) = -t47;
coef(34) = -t41;
coef(35) = -t48;
coef(36) = -t28 * c(1);
coef(37) = t31;
coef(38) = t21;
coef(39) = t34;
coef(40) = t9;
coef(41) = t25;
coef(42) = -t45;
coef(43) = -t43;
%
t1 = bb(3) * dd(3);
t3 = bb(2) ^ 2;
t5 = bb(3) ^ 2;
t7 = bb(1) ^ 2;
t9 = bb(1) * dd(1);
t11 = bb(2) * dd(2);
t14 = sqrt(0.3e1);
t16 = t5 * t14 / 0.12e2;
t18 = t3 * t14 / 0.12e2;
t19 = cc(1) * bb(1);
t21 = t19 * t14 / 0.6e1;
t22 = cc(2) * bb(2);
t24 = t22 * t14 / 0.6e1;
t26 = t7 * t14 / 0.12e2;
t27 = cc(3) * bb(3);
t29 = t27 * t14 / 0.6e1;
t30 = -t16 - t18 + t21 + t24 - t26 + t29;
t32 = t1 * t14 / 0.6e1;
t34 = t11 * t14 / 0.6e1;
t36 = t9 * t14 / 0.6e1;
t37 = t26 + t18 + t16 - t32 - t34 - t36;
t38 = -t3 - t7 - t5;
t54 = -t1 / 0.6e1 - t27 / 0.6e1 - t11 / 0.6e1 + t7 / 0.12e2 - t9 / 0.6e1 - t19 / 0.6e1 + t3 / 0.12e2 + t5 / 0.12e2 + cc(1) * dd(1) / 0.3e1 - t22 / 0.6e1 + cc(2) * dd(2) / 0.3e1 + cc(3) * dd(3) / 0.3e1;
t64 = -t16 - t26 - t14 * cc(2) * dd(2) / 0.3e1 + t24 - t18 + t34 + t29 + t36 - t14 * cc(3) * dd(3) / 0.3e1 + t32 + t21 - t14 * cc(1) * dd(1) / 0.3e1;
%
cm(1) = t1 / 0.2e1 - t3 / 0.4e1 - t5 / 0.4e1 - t7 / 0.4e1 + t9 / 0.2e1 + t11 / 0.2e1;
cm(2) = t30;
cm(3) = t37;
cm(4) = t38 / 0.4e1;
cm(5) = t37;
cm(6) = t54;
cm(7) = t38 / 0.4e1;
cm(8) = -t30;
cm(9) = -t38 / 0.4e1;
cm(10) = t64;
cm(11) = t9 / 0.2e1 + t11 / 0.2e1 + t1 / 0.2e1;
%
t1 = ko ^ 2;
t3 = j ^ 2;
t5 = (0.1e1 / t1 / t3);
t10 = (0.1e1 / ko / j);
t22 = (0.1e1 / t1 / ko / t3 / j);
t33 = (t5 * cm(1));
t34 = (t5 * cm(9));
t36 = 2 * t22 * cm(7);
t41 = t1 ^ 2;
t43 = t3 ^ 2;
t45 = (0.1e1 / t41 / t43);
t47 = 6 * t45 * cm(6);
t50 = 2 * t22 * cm(4);
t53 = 2 * t22 * cm(3);
t57 = 6 * t45 * cm(2);
t62 = 2 * t22 * cm(8);
t63 = t5 * cm(11);
t65 = 2 * t22 * cm(10);
t68 = 2 * t22 * cm(5);
%
coefm(1) = -2 * t5 * cm(8);
coefm(2) = -t10 * cm(8);
coefm(3) = -t10 * cm(6);
coefm(4) = -t10 * cm(2);
coefm(5) = -2 * t5 * cm(4);
coefm(6) = -2 * t5 * cm(3);
coefm(7) = -6 * t22 * cm(6);
coefm(8) = -6 * t22 * cm(2);
coefm(9) = -3 * t5 * cm(6);
coefm(10) = -3 * t5 * cm(2);
coefm(11) = -t10 * cm(4);
coefm(12) = -t10 * cm(3);
coefm(13) = -t33;
coefm(14) = t34;
coefm(15) = t36;
coefm(16) = -2 * t5 * cm(5);
coefm(17) = -t10 * cm(5);
coefm(18) = -t10 * cm(9);
coefm(19) = -t47;
coefm(20) = -t10 * cm(10);
coefm(21) = -t50;
coefm(22) = -t10 * cm(7);
coefm(23) = -t53;
coefm(24) = -2 * t5 * cm(10);
coefm(25) = -t57;
coefm(26) = -t10 * cm(11);
coefm(27) = -2 * t5 * cm(7);
coefm(28) = -t62;
coefm(29) = t63;
coefm(30) = t47;
coefm(31) = t62;
coefm(32) = -t65;
coefm(33) = -t10 * cm(1);
coefm(34) = -t63;
coefm(35) = -t34;
coefm(36) = t53;
coefm(37) = t50;
coefm(38) = -t36;
coefm(39) = -t68;
coefm(40) = t57;
coefm(41) = t68;
coefm(42) = t33;
coefm(43) = t65;

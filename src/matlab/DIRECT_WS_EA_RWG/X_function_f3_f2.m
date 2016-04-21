function X = X_function_f3_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f3_f2

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
% theta, Psi, B, Bm, coef, coefm, N, Nm

% OUTPUT DATA
% X   
%
t2 = B ^ 2;
t3 = t2 ^ 2;
t5 = 0.1e1 / t3 / B;
t8 = cos(Psi);
t9 = t8 ^ 2;
t10 = cos(theta);
t11 = t9 * t10;
t12 = sin(Psi);
t13 = t11 * t12;
t15 = 0.1e1 / t3;
t18 = t8 * t12;
t21 = 0.1e1 / t2 / B;
t31 = sin(theta);
t32 = t9 * t12 * t31;
t36 = t9 * t31;
t56 = 0.1e1 / t2;
t182 = N(1) * (t5 * coef(4) * t13 + t15 * coef(3) * t18 + t21 * coef(39) * t8 + t15 * coef(42) * t11 + t5 * coef(5) * t32 + t15 * coef(41) * t36) + N(8) * (t15 * coef(34) * t11 + t15 * coef(35) * t36 + t21 * coef(13) * t8 + t15 * coef(30) * t18) + N(11) * t56 * t8 * coef(33) + N(2) * (t15 * coef(9) * t36 + t15 * coef(10) * t11 + t5 * coef(31) * t13 + t5 * coef(32) * t32 + t21 * coef(1) * t8 + t15 * coef(11) * t18) + N(12) * (t56 * coef(15) * t11 + t56 * coef(16) * t36 + t56 * coef(18) * t18) + N(10) * (t21 * coef(21) * t11 + t21 * coef(22) * t36 + t56 * coef(40) * t8 + t21 * coef(17) * t18) + N(5) * (t56 * coef(19) * t13 + t56 * coef(20) * t32) + N(3) * (t21 * coef(36) * t11 + t21 * coef(37) * t36 + t15 * coef(23) * t13 + t15 * coef(24) * t32 + t56 * coef(14) * t8 + t21 * coef(29) * t18) + N(4) * (t56 * coef(27) * t11 + t56 * coef(38) * t36 + t21 * coef(25) * t13 + t21 * coef(26) * t32 + t56 * coef(28) * t18) + N(9) * t21 * t8 * coef(12) + N(6) * (t15 * coef(7) * t11 + t15 * coef(8) * t36 + t21 * coef(43) * t8 + t15 * coef(6) * t18) + N(7) * t21 * t8 * coef(2);
t184 = Bm ^ 2;
t186 = 0.1e1 / t184 / Bm;
t190 = t184 ^ 2;
t191 = 0.1e1 / t190;
t196 = 0.1e1 / t190 / Bm;
t253 = 0.1e1 / t184;
t353 = Nm(1) * (t186 * coefm(30) * t8 + t191 * coefm(15) * t11 + t196 * coefm(4) * t32 + t196 * coefm(5) * t13 + t191 * coefm(31) * t18 + t191 * coefm(13) * t36) + Nm(6) * (t191 * coefm(7) * t11 + t191 * coefm(8) * t36 + t186 * coefm(1) * t8 + t191 * coefm(6) * t18) + Nm(9) * t186 * t8 * coefm(10) + Nm(8) * (t191 * coefm(26) * t11 + t191 * coefm(27) * t36 + t186 * coefm(11) * t8 + t191 * coefm(23) * t18) + Nm(10) * (t186 * coefm(33) * t11 + t186 * coefm(34) * t36 + t253 * coefm(21) * t8 + t186 * coefm(42) * t18) + Nm(12) * (t253 * coefm(39) * t36 + t253 * coefm(43) * t11 + t253 * coefm(41) * t18) + Nm(2) * (t191 * coefm(2) * t36 + t191 * coefm(3) * t11 + t196 * coefm(19) * t13 + t196 * coefm(20) * t32 + t186 * coefm(14) * t8 + t191 * coefm(9) * t18) + Nm(5) * (t253 * coefm(32) * t13 + t253 * coefm(40) * t32) + Nm(3) * (t186 * coefm(25) * t36 + t186 * coefm(28) * t11 + t191 * coefm(35) * t32 + t191 * coefm(36) * t13 + t253 * coefm(12) * t8 + t186 * coefm(29) * t18) + Nm(4) * (t253 * coefm(22) * t11 + t253 * coefm(24) * t36 + t186 * coefm(37) * t32 + t186 * coefm(38) * t13 + t253 * coefm(18) * t18) + Nm(11) * t253 * t8 * coefm(17) + Nm(7) * t186 * t8 * coefm(16);
% Output
X = t353 + t182;
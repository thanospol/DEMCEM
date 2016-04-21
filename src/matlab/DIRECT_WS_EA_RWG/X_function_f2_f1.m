function X = X_function_f2_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f2_f1

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
t10 = sin(Psi);
t12 = sin(theta);
t13 = t9 * t10 * t12;
t15 = 0.1e1 / t3;
t18 = t8 * t10;
t21 = 0.1e1 / t2 / B;
t27 = cos(theta);
t28 = t9 * t27;
t32 = t28 * t10;
t36 = t9 * t12;
t41 = 0.1e1 / t2;
t182 = N(1) * (t5 * coef(7) * t13 + t15 * coef(43) * t18 + t21 * coef(42) * t8 + t15 * coef(2) * t28 + t5 * coef(6) * t32 + t15 * coef(4) * t36) + N(12) * (t41 * coef(17) * t28 + t41 * coef(18) * t36 + t41 * coef(20) * t18) + N(10) * (t21 * coef(23) * t28 + t21 * coef(24) * t36 + t41 * coef(37) * t8 + t21 * coef(19) * t18) + N(5) * (t41 * coef(21) * t32 + t41 * coef(22) * t13) + N(3) * (t21 * coef(31) * t36 + t21 * coef(34) * t28 + t15 * coef(25) * t32 + t15 * coef(26) * t13 + t41 * coef(16) * t8 + t21 * coef(39) * t18) + N(4) * (t41 * coef(29) * t28 + t41 * coef(32) * t36 + t21 * coef(27) * t32 + t21 * coef(28) * t13 + t41 * coef(38) * t18) + N(9) * t21 * t8 * coef(14) + N(8) * (t15 * coef(30) * t28 + t15 * coef(33) * t36 + t21 * coef(15) * t8 + t15 * coef(35) * t18) + N(11) * t41 * t8 * coef(36) + N(2) * (t15 * coef(11) * t36 + t15 * coef(12) * t28 + t5 * coef(40) * t32 + t5 * coef(41) * t13 + t21 * coef(5) * t8 + t15 * coef(13) * t18) + N(6) * (t15 * coef(9) * t28 + t15 * coef(10) * t36 + t21 * coef(3) * t8 + t15 * coef(8) * t18) + N(7) * t21 * t8 * coef(1);
t184 = Bm ^ 2;
t186 = 0.1e1 / t184 / Bm;
t190 = t184 ^ 2;
t192 = 0.1e1 / t190 / Bm;
t196 = 0.1e1 / t190;
t217 = 0.1e1 / t184;
t353 = Nm(1) * (t186 * coefm(15) * t8 + t192 * coefm(4) * t13 + t196 * coefm(30) * t28 + t196 * coefm(43) * t18 + t196 * coefm(1) * t36 + t192 * coefm(5) * t32) + Nm(7) * t186 * t8 * coefm(2) + Nm(5) * (t217 * coefm(35) * t13 + t217 * coefm(36) * t32) + Nm(3) * (t186 * coefm(22) * t36 + t186 * coefm(23) * t28 + t196 * coefm(39) * t13 + t196 * coefm(40) * t32 + t217 * coefm(8) * t8 + t186 * coefm(17) * t18) + Nm(4) * (t217 * coefm(26) * t36 + t217 * coefm(27) * t28 + t186 * coefm(41) * t13 + t186 * coefm(42) * t32 + t217 * coefm(18) * t18) + Nm(6) * (t196 * coefm(11) * t36 + t196 * coefm(12) * t28 + t186 * coefm(3) * t8 + t196 * coefm(10) * t18) + Nm(2) * (t196 * coefm(13) * t28 + t196 * coefm(14) * t36 + t192 * coefm(24) * t32 + t192 * coefm(28) * t13 + t186 * coefm(29) * t8 + t196 * coefm(9) * t18) + Nm(10) * (t186 * coefm(37) * t28 + t186 * coefm(38) * t36 + t217 * coefm(19) * t8 + t186 * coefm(33) * t18) + Nm(8) * (t196 * coefm(20) * t36 + t196 * coefm(21) * t28 + t186 * coefm(6) * t8 + t196 * coefm(16) * t18) + Nm(11) * t217 * t8 * coefm(25) + Nm(9) * t186 * t8 * coefm(7) + Nm(12) * (t217 * coefm(31) * t28 + t217 * coefm(32) * t36 + t217 * coefm(34) * t18);
% Output
X = t353 + t182;
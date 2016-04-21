function X = X_function_f1_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f1_f2

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
t4 = 0.1e1 / t2 / B;
t7 = cos(Psi);
t9 = t2 ^ 2;
t10 = 0.1e1 / t9;
t13 = sin(Psi);
t14 = t7 * t13;
t18 = t7 ^ 2;
t19 = sin(theta);
t20 = t18 * t19;
t24 = cos(theta);
t25 = t18 * t24;
t28 = 0.1e1 / t9 / B;
t31 = t25 * t13;
t36 = t18 * t13 * t19;
t61 = 0.1e1 / t2;
t182 = N(1) * (t4 * coef(30) * t7 + t10 * coef(2) * t14 + t10 * coef(3) * t20 + t10 * coef(4) * t25 + t28 * coef(6) * t31 + t28 * coef(7) * t36) + N(6) * (t10 * coef(8) * t25 + t10 * coef(9) * t20 + t4 * coef(1) * t7 + t10 * coef(10) * t14) + N(7) * t4 * t7 * coef(5) + N(4) * (t61 * coef(31) * t20 + t61 * coef(43) * t25 + t4 * coef(28) * t31 + t4 * coef(29) * t36 + t61 * coef(38) * t14) + N(3) * (t4 * coef(37) * t25 + t4 * coef(39) * t20 + t10 * coef(26) * t31 + t10 * coef(27) * t36 + t61 * coef(14) * t7 + t4 * coef(40) * t14) + N(8) * (t10 * coef(33) * t20 + t10 * coef(36) * t25 + t4 * coef(13) * t7 + t10 * coef(32) * t14) + N(10) * (t4 * coef(24) * t25 + t4 * coef(25) * t20 + t61 * coef(34) * t7 + t4 * coef(21) * t14) + N(2) * (t10 * coef(11) * t20 + t10 * coef(12) * t25 + t28 * coef(35) * t31 + t28 * coef(42) * t36 + t4 * coef(17) * t7 + t10 * coef(15) * t14) + N(12) * (t61 * coef(18) * t25 + t61 * coef(19) * t20 + t61 * coef(20) * t14) + N(5) * (t61 * coef(22) * t31 + t61 * coef(23) * t36) + N(9) * t4 * t7 * coef(16) + N(11) * t61 * t7 * coef(41);
t184 = Bm ^ 2;
t186 = 0.1e1 / t184 / Bm;
t190 = t184 ^ 2;
t191 = 0.1e1 / t190;
t199 = 0.1e1 / t190 / Bm;
t233 = 0.1e1 / t184;
t353 = Nm(1) * (t186 * coefm(1) * t7 + t191 * coefm(15) * t25 + t191 * coefm(30) * t14 + t199 * coefm(9) * t36 + t191 * coefm(14) * t20 + t199 * coefm(8) * t31) + Nm(6) * (t191 * coefm(5) * t25 + t191 * coefm(6) * t20 + t186 * coefm(13) * t7 + t191 * coefm(4) * t14) + Nm(10) * (t186 * coefm(22) * t25 + t186 * coefm(23) * t20 + t233 * coefm(32) * t7 + t186 * coefm(18) * t14) + Nm(11) * t233 * t7 * coefm(34) + Nm(8) * (t191 * coefm(35) * t20 + t191 * coefm(41) * t25 + t186 * coefm(2) * t7 + t191 * coefm(37) * t14) + Nm(4) * (t233 * coefm(31) * t25 + t233 * coefm(39) * t20 + t186 * coefm(26) * t36 + t186 * coefm(27) * t31 + t233 * coefm(40) * t14) + Nm(3) * (t186 * coefm(33) * t25 + t186 * coefm(38) * t20 + t191 * coefm(24) * t36 + t191 * coefm(25) * t31 + t233 * coefm(3) * t7 + t186 * coefm(36) * t14) + Nm(9) * t186 * t7 * coefm(12) + Nm(7) * t186 * t7 * coefm(29) + Nm(2) * (t191 * coefm(7) * t20 + t191 * coefm(10) * t25 + t199 * coefm(42) * t36 + t199 * coefm(43) * t31 + t186 * coefm(28) * t7 + t191 * coefm(11) * t14) + Nm(5) * (t233 * coefm(20) * t36 + t233 * coefm(21) * t31) + Nm(12) * (t233 * coefm(16) * t25 + t233 * coefm(17) * t20 + t233 * coefm(19) * t14);
% Output
X = t353 + t182;
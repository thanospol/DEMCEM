function X = X_function_f3_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f3_f1

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
t13 = t7 ^ 2;
t14 = sin(theta);
t15 = t13 * t14;
t19 = cos(theta);
t20 = t13 * t19;
t24 = sin(Psi);
t25 = t7 * t24;
t28 = 0.1e1 / t9 / B;
t32 = t13 * t24 * t14;
t36 = t20 * t24;
t76 = 0.1e1 / t2;
t182 = N(1) * (t4 * coef(30) * t7 + t10 * coef(6) * t15 + t10 * coef(3) * t20 + t10 * coef(4) * t25 + t28 * coef(8) * t32 + t28 * coef(7) * t36) + N(6) * (t10 * coef(9) * t20 + t10 * coef(10) * t15 + t4 * coef(2) * t7 + t10 * coef(11) * t25) + N(7) * t4 * t7 * coef(5) + N(8) * (t10 * coef(39) * t15 + t10 * coef(42) * t20 + t4 * coef(17) * t7 + t10 * coef(31) * t25) + N(12) * (t76 * coef(18) * t20 + t76 * coef(19) * t15 + t76 * coef(21) * t25) + N(5) * (t76 * coef(22) * t36 + t76 * coef(23) * t32) + N(2) * (t10 * coef(12) * t15 + t10 * coef(13) * t20 + t28 * coef(32) * t36 + t28 * coef(35) * t32 + t4 * coef(1) * t7 + t10 * coef(14) * t25) + N(9) * t4 * t7 * coef(16) + N(3) * (t4 * coef(41) * t20 + t4 * coef(43) * t15 + t10 * coef(26) * t36 + t10 * coef(27) * t32 + t76 * coef(15) * t7 + t4 * coef(33) * t25) + N(4) * (t76 * coef(34) * t15 + t76 * coef(40) * t20 + t4 * coef(28) * t36 + t4 * coef(29) * t32 + t76 * coef(37) * t25) + N(10) * (t4 * coef(24) * t20 + t4 * coef(25) * t15 + t76 * coef(36) * t7 + t4 * coef(20) * t25) + N(11) * t76 * t7 * coef(38);
t184 = Bm ^ 2;
t185 = t184 ^ 2;
t187 = 0.1e1 / t185 / Bm;
t195 = 0.1e1 / t184 / Bm;
t199 = 0.1e1 / t185;
t232 = 0.1e1 / t184;
t353 = Nm(1) * (t187 * coefm(23) * t32 + t187 * coefm(24) * t36 + t195 * coefm(13) * t7 + t199 * coefm(28) * t15 + t199 * coefm(29) * t20 + t199 * coefm(25) * t25) + Nm(6) * (t199 * coefm(18) * t20 + t199 * coefm(19) * t15 + t195 * coefm(30) * t7 + t199 * coefm(20) * t25) + Nm(7) * t195 * t7 * coefm(27) + Nm(12) * (t232 * coefm(1) * t20 + t232 * coefm(2) * t15 + t232 * coefm(3) * t25) + Nm(10) * (t195 * coefm(7) * t20 + t195 * coefm(8) * t15 + t232 * coefm(39) * t7 + t195 * coefm(4) * t25) + Nm(9) * t195 * t7 * coefm(15) + Nm(8) * (t199 * coefm(33) * t15 + t199 * coefm(35) * t20 + t195 * coefm(16) * t7 + t199 * coefm(36) * t25) + Nm(11) * t232 * t7 * coefm(43) + Nm(5) * (t232 * coefm(5) * t32 + t232 * coefm(6) * t36) + Nm(2) * (t199 * coefm(21) * t15 + t199 * coefm(22) * t20 + t187 * coefm(31) * t36 + t187 * coefm(41) * t32 + t195 * coefm(26) * t7 + t199 * coefm(14) * t25) + Nm(3) * (t195 * coefm(32) * t20 + t195 * coefm(42) * t15 + t199 * coefm(9) * t32 + t199 * coefm(10) * t36 + t232 * coefm(17) * t7 + t195 * coefm(38) * t25) + Nm(4) * (t232 * coefm(34) * t15 + t232 * coefm(40) * t20 + t195 * coefm(11) * t32 + t195 * coefm(12) * t36 + t232 * coefm(37) * t25);
% Output
X = t353 + t182;
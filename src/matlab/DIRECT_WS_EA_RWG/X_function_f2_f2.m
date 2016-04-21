function X = X_function_f2_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f2_f2

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
t11 = 0.1e1 / t9 / B;
t14 = t7 ^ 2;
t15 = cos(theta);
t16 = t14 * t15;
t17 = sin(Psi);
t18 = t16 * t17;
t20 = 0.1e1 / t9;
t23 = t7 * t17;
t31 = sin(theta);
t32 = t14 * t17 * t31;
t36 = t14 * t31;
t58 = 0.1e1 / t2;
t173 = N(1) * (t4 * coef(27) * t7 + t11 * coef(23) * t18 + t20 * coef(1) * t23 + t20 * coef(2) * t16 + t11 * coef(22) * t32 + t20 * coef(4) * t36) + N(7) * t4 * t7 * coef(40) + N(8) * (t20 * coef(30) * t36 + t20 * coef(37) * t16 + t20 * coef(34) * t23) + N(4) * (t58 * coef(31) * t36 + t58 * coef(32) * t16 + t4 * coef(15) * t32 + t4 * coef(16) * t18 + t58 * coef(38) * t23) + N(3) * (t4 * coef(28) * t16 + t4 * coef(35) * t36 + t20 * coef(13) * t32 + t20 * coef(14) * t18 + t58 * coef(18) * t7 + t4 * coef(33) * t23) + N(11) * t58 * t7 * coef(39) + N(12) * (t58 * coef(5) * t16 + t58 * coef(9) * t36 + t58 * coef(8) * t23) + N(5) * (t58 * coef(6) * t32 + t58 * coef(10) * t18) + N(10) * (t4 * coef(11) * t36 + t4 * coef(12) * t16 + t4 * coef(7) * t23) + N(9) * t4 * t7 * coef(19) + N(2) * (t20 * coef(20) * t36 + t20 * coef(21) * t16 + t11 * coef(29) * t18 + t11 * coef(36) * t32 + t4 * coef(3) * t7 + t20 * coef(17) * t23) + N(6) * (t20 * coef(25) * t36 + t20 * coef(26) * t16 + t20 * coef(24) * t23);
t175 = Bm ^ 2;
t177 = 0.1e1 / t175 / Bm;
t181 = t175 ^ 2;
t182 = 0.1e1 / t181;
t190 = 0.1e1 / t181 / Bm;
t203 = 0.1e1 / t175;
t335 = Nm(1) * (t177 * coefm(36) * t7 + t182 * coefm(38) * t23 + t182 * coefm(39) * t16 + t190 * coefm(10) * t32 + t182 * coefm(23) * t36 + t190 * coefm(1) * t18) + Nm(12) * (t203 * coefm(24) * t36 + t203 * coefm(25) * t16 + t203 * coefm(27) * t23) + Nm(10) * (t177 * coefm(34) * t36 + t177 * coefm(35) * t16 + t177 * coefm(26) * t23) + Nm(4) * (t203 * coefm(15) * t36 + t203 * coefm(18) * t16 + t177 * coefm(32) * t32 + t177 * coefm(33) * t18 + t203 * coefm(16) * t23) + Nm(3) * (t177 * coefm(19) * t16 + t177 * coefm(20) * t36 + t182 * coefm(30) * t32 + t182 * coefm(31) * t18 + t203 * coefm(6) * t7 + t177 * coefm(17) * t23) + Nm(5) * (t203 * coefm(28) * t32 + t203 * coefm(29) * t18) + Nm(11) * t203 * t7 * coefm(13) + Nm(8) * (t182 * coefm(11) * t16 + t182 * coefm(21) * t36 + t182 * coefm(22) * t23) + Nm(6) * (t182 * coefm(3) * t36 + t182 * coefm(4) * t16 + t182 * coefm(2) * t23) + Nm(2) * (t182 * coefm(8) * t36 + t182 * coefm(9) * t16 + t190 * coefm(12) * t18 + t190 * coefm(14) * t32 + t177 * coefm(37) * t7 + t182 * coefm(7) * t23) + Nm(7) * t177 * t7 * coefm(40) + Nm(9) * t177 * t7 * coefm(5);
% Output
X = t335 + t173;

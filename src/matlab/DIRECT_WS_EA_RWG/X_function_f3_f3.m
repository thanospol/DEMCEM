function X = X_function_f3_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f3_f3

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
t22 = t9 * t12;
t26 = cos(theta);
t27 = t9 * t26;
t31 = t27 * t10;
t34 = 0.1e1 / t2 / B;
t56 = 0.1e1 / t2;
t182 = N(1) * (t5 * coef(20) * t13 + t15 * coef(1) * t18 + t15 * coef(5) * t22 + t15 * coef(6) * t27 + t5 * coef(19) * t31 + t34 * coef(30) * t8) + N(6) * (t15 * coef(22) * t27 + t15 * coef(23) * t22 + t34 * coef(4) * t8 + t15 * coef(21) * t18) + N(5) * (t56 * coef(11) * t31 + t56 * coef(12) * t13) + N(10) * (t34 * coef(13) * t27 + t34 * coef(14) * t22 + t56 * coef(37) * t8 + t34 * coef(10) * t18) + N(3) * (t34 * coef(34) * t27 + t34 * coef(38) * t22 + t15 * coef(15) * t31 + t15 * coef(16) * t13 + t56 * coef(29) * t8 + t34 * coef(41) * t18) + N(4) * (t56 * coef(31) * t27 + t56 * coef(36) * t22 + t34 * coef(17) * t31 + t34 * coef(18) * t13 + t56 * coef(35) * t18) + N(12) * (t56 * coef(7) * t27 + t56 * coef(8) * t22 + t56 * coef(9) * t18) + N(9) * t34 * t8 * coef(27) + N(8) * (t15 * coef(39) * t22 + t15 * coef(40) * t27 + t34 * coef(28) * t8 + t15 * coef(33) * t18) + N(7) * t34 * t8 * coef(3) + N(11) * t56 * t8 * coef(32) + N(2) * (t15 * coef(24) * t22 + t15 * coef(25) * t27 + t5 * coef(42) * t31 + t5 * coef(43) * t13 + t34 * coef(2) * t8 + t15 * coef(26) * t18);
t184 = Bm ^ 2;
t185 = t184 ^ 2;
t186 = 0.1e1 / t185;
t191 = 0.1e1 / t185 / Bm;
t199 = 0.1e1 / t184 / Bm;
t227 = 0.1e1 / t184;
t353 = Nm(1) * (t186 * coefm(41) * t18 + t191 * coefm(19) * t31 + t186 * coefm(43) * t22 + t199 * coefm(39) * t8 + t186 * coefm(38) * t27 + t191 * coefm(18) * t13) + Nm(6) * (t186 * coefm(15) * t27 + t186 * coefm(16) * t22 + t199 * coefm(14) * t8 + t186 * coefm(23) * t18) + Nm(5) * (t227 * coefm(30) * t13 + t227 * coefm(31) * t31) + Nm(3) * (t199 * coefm(3) * t22 + t199 * coefm(7) * t27 + t186 * coefm(34) * t13 + t186 * coefm(35) * t31 + t227 * coefm(21) * t8 + t199 * coefm(1) * t18) + Nm(4) * (t227 * coefm(12) * t27 + t227 * coefm(13) * t22 + t199 * coefm(36) * t13 + t199 * coefm(37) * t31 + t227 * coefm(11) * t18) + Nm(9) * t199 * t8 * coefm(17) + Nm(8) * (t186 * coefm(4) * t27 + t186 * coefm(5) * t22 + t199 * coefm(20) * t8 + t186 * coefm(10) * t18) + Nm(7) * t199 * t8 * coefm(42) + Nm(12) * (t227 * coefm(26) * t27 + t227 * coefm(27) * t22 + t227 * coefm(28) * t18) + Nm(10) * (t199 * coefm(32) * t27 + t199 * coefm(33) * t22 + t227 * coefm(8) * t8 + t199 * coefm(29) * t18) + Nm(11) * t227 * t8 * coefm(2) + Nm(2) * (t186 * coefm(24) * t22 + t186 * coefm(25) * t27 + t191 * coefm(6) * t31 + t191 * coefm(9) * t13 + t199 * coefm(40) * t8 + t186 * coefm(22) * t18);
% Output
X = t353 + t182;

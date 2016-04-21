function X = X_function_f2_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f2_f3

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
t4 = 0.1e1 / t3;
t7 = cos(Psi);
t8 = t7 ^ 2;
t9 = sin(theta);
t10 = t8 * t9;
t13 = 0.1e1 / t2 / B;
t19 = sin(Psi);
t20 = t7 * t19;
t23 = 0.1e1 / t3 / B;
t26 = cos(theta);
t27 = t8 * t26;
t28 = t27 * t19;
t33 = t8 * t19 * t9;
t67 = 0.1e1 / t2;
t182 = N(1) * (t4 * coef(17) * t10 + t13 * coef(14) * t7 + t4 * coef(19) * t20 + t23 * coef(37) * t28 + t23 * coef(38) * t33 + t4 * coef(16) * t27) + N(7) * t13 * t7 * coef(15) + N(2) * (t4 * coef(42) * t10 + t4 * coef(43) * t27 + t23 * coef(3) * t33 + t23 * coef(7) * t28 + t13 * coef(18) * t7 + t4 * coef(33) * t20) + N(12) * (t67 * coef(21) * t27 + t67 * coef(22) * t10 + t67 * coef(23) * t20) + N(10) * (t13 * coef(27) * t27 + t13 * coef(28) * t10 + t67 * coef(11) * t7 + t13 * coef(24) * t20) + N(5) * (t67 * coef(25) * t28 + t67 * coef(26) * t33) + N(6) * (t4 * coef(40) * t27 + t4 * coef(41) * t10 + t13 * coef(20) * t7 + t4 * coef(39) * t20) + N(8) * (t4 * coef(1) * t27 + t4 * coef(5) * t10 + t13 * coef(35) * t7 + t4 * coef(9) * t20) + N(11) * t67 * t7 * coef(10) + N(9) * t13 * t7 * coef(34) + N(3) * (t13 * coef(2) * t10 + t13 * coef(4) * t27 + t4 * coef(29) * t28 + t4 * coef(30) * t33 + t67 * coef(36) * t7 + t13 * coef(13) * t20) + N(4) * (t67 * coef(6) * t10 + t67 * coef(8) * t27 + t13 * coef(31) * t28 + t13 * coef(32) * t33 + t67 * coef(12) * t20);
t184 = Bm ^ 2;
t185 = t184 ^ 2;
t187 = 0.1e1 / t185 / Bm;
t191 = 0.1e1 / t185;
t199 = 0.1e1 / t184 / Bm;
t239 = 0.1e1 / t184;
t353 = Nm(1) * (t187 * coefm(30) * t33 + t191 * coefm(41) * t10 + t191 * coefm(43) * t20 + t199 * coefm(42) * t7 + t191 * coefm(15) * t27 + t187 * coefm(40) * t28) + Nm(2) * (t191 * coefm(38) * t27 + t191 * coefm(39) * t10 + t187 * coefm(19) * t33 + t187 * coefm(25) * t28 + t199 * coefm(13) * t7 + t191 * coefm(32) * t20) + Nm(10) * (t199 * coefm(5) * t27 + t199 * coefm(6) * t10 + t239 * coefm(26) * t7 + t199 * coefm(1) * t20) + Nm(5) * (t239 * coefm(3) * t33 + t239 * coefm(4) * t28) + Nm(8) * (t191 * coefm(21) * t27 + t191 * coefm(23) * t10 + t199 * coefm(34) * t7 + t191 * coefm(28) * t20) + Nm(11) * t239 * t7 * coefm(18) + Nm(7) * t199 * t7 * coefm(14) + Nm(9) * t199 * t7 * coefm(35) + Nm(12) * (t239 * coefm(11) * t27 + t239 * coefm(12) * t10 + t239 * coefm(2) * t20) + Nm(3) * (t199 * coefm(16) * t10 + t199 * coefm(27) * t27 + t191 * coefm(7) * t33 + t191 * coefm(8) * t28 + t239 * coefm(33) * t7 + t199 * coefm(24) * t20) + Nm(4) * (t239 * coefm(17) * t10 + t239 * coefm(22) * t27 + t199 * coefm(9) * t33 + t199 * coefm(10) * t28 + t239 * coefm(20) * t20) + Nm(6) * (t191 * coefm(36) * t10 + t191 * coefm(37) * t27 + t199 * coefm(29) * t7 + t191 * coefm(31) * t20);
% Output
X = t353 + t182;

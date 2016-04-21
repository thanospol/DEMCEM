function X = X_function_f1_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f1_f1

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
t18 = sin(theta);
t19 = t9 * t12 * t18;
t21 = 0.1e1 / t3;
t24 = t9 * t18;
t28 = t8 * t12;
t34 = 0.1e1 / t2 / B;
t53 = 0.1e1 / t2;
t173 = N(1) * (t5 * coef(18) * t13 + t5 * coef(19) * t19 + t21 * coef(1) * t24 + t21 * coef(25) * t28 + t21 * coef(38) * t11 + t34 * coef(39) * t8) + N(3) * (t34 * coef(4) * t24 + t34 * coef(13) * t11 + t21 * coef(34) * t13 + t21 * coef(35) * t19 + t53 * coef(20) * t8 + t34 * coef(8) * t28) + N(9) * t34 * t8 * coef(22) + N(8) * (t21 * coef(9) * t24 + t21 * coef(12) * t11 + t21 * coef(10) * t28) + N(6) * (t21 * coef(14) * t24 + t21 * coef(15) * t11 + t21 * coef(23) * t28) + N(7) * t34 * t8 * coef(40) + N(2) * (t21 * coef(16) * t24 + t21 * coef(17) * t11 + t5 * coef(2) * t19 + t5 * coef(11) * t13 + t34 * coef(24) * t8 + t21 * coef(21) * t28) + N(11) * t53 * t8 * coef(5) + N(12) * (t53 * coef(26) * t24 + t53 * coef(27) * t11 + t53 * coef(28) * t28) + N(10) * (t34 * coef(32) * t24 + t34 * coef(33) * t11 + t34 * coef(29) * t28) + N(5) * (t53 * coef(30) * t13 + t53 * coef(31) * t19) + N(4) * (t53 * coef(3) * t24 + t53 * coef(7) * t11 + t34 * coef(36) * t13 + t34 * coef(37) * t19 + t53 * coef(6) * t28);
t175 = Bm ^ 2;
t176 = t175 ^ 2;
t178 = 0.1e1 / t176 / Bm;
t182 = 0.1e1 / t176;
t193 = 0.1e1 / t175 / Bm;
t215 = 0.1e1 / t175;
t335 = Nm(1) * (t178 * coefm(39) * t19 + t182 * coefm(28) * t24 + t182 * coefm(29) * t11 + t182 * coefm(13) * t28 + t193 * coefm(26) * t8 + t178 * coefm(40) * t13) + Nm(6) * (t182 * coefm(34) * t24 + t182 * coefm(37) * t11 + t182 * coefm(38) * t28) + Nm(12) * (t215 * coefm(24) * t24 + t215 * coefm(25) * t11 + t215 * coefm(18) * t28) + Nm(7) * t193 * t8 * coefm(30) + Nm(8) * (t182 * coefm(2) * t24 + t182 * coefm(4) * t11 + t182 * coefm(10) * t28) + Nm(9) * t193 * t8 * coefm(31) + Nm(3) * (t193 * coefm(1) * t11 + t193 * coefm(9) * t24 + t182 * coefm(19) * t19 + t182 * coefm(20) * t13 + t215 * coefm(33) * t8 + t193 * coefm(12) * t28) + Nm(11) * t215 * t8 * coefm(8) + Nm(2) * (t182 * coefm(35) * t11 + t182 * coefm(36) * t24 + t178 * coefm(3) * t13 + t178 * coefm(7) * t19 + t193 * coefm(27) * t8 + t182 * coefm(32) * t28) + Nm(10) * (t193 * coefm(15) * t24 + t193 * coefm(16) * t11 + t193 * coefm(17) * t28) + Nm(5) * (t215 * coefm(14) * t13 + t215 * coefm(23) * t19) + Nm(4) * (t215 * coefm(6) * t11 + t215 * coefm(11) * t24 + t193 * coefm(21) * t19 + t193 * coefm(22) * t13 + t215 * coefm(5) * t28);
% Output
X = t335 + t173;
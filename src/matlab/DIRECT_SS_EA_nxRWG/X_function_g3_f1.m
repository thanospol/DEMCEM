function X = X_function_g3_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g3_f1

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
t16 = 0.1e1 / t3 / t2;
t19 = t9 * t12;
t23 = t9 ^ 2;
t24 = t23 * t12;
t26 = 0.1e1 / t3;
t53 = 0.1e1 / t2 / B;
t133 = N(1) * (t5 * coef(14) * t13 + t16 * coef(9) * t19 + t16 * coef(8) * t24 + t26 * coef(27) * t19) + N(6) * (t26 * coef(15) * t19 + t5 * coef(23) * t13) + N(7) * t26 * t19 * coef(10) + N(9) * t26 * t19 * coef(22) + N(4) * (t53 * coef(1) * t13 + t26 * coef(21) * t24 + t26 * coef(20) * t19) + N(5) * (t53 * coef(17) * t24 + t53 * coef(18) * t19) + N(3) * (t53 * coef(16) * t19 + t26 * coef(4) * t13 + t5 * coef(19) * t24 + t5 * coef(26) * t19) + N(8) * (t26 * coef(24) * t19 + t5 * coef(5) * t13) + N(2) * (t26 * coef(11) * t19 + t5 * coef(25) * t13 + t16 * coef(12) * t24 + t16 * coef(13) * t19) + N(10) * (t53 * coef(2) * t19 + t26 * coef(6) * t13) + N(12) * t10 * t53 * t19 * coef(7) + N(11) * t53 * t19 * coef(3);
t135 = Bm ^ 2;
t136 = t135 ^ 2;
t138 = 0.1e1 / t136 / t135;
t143 = 0.1e1 / t136 / Bm;
t150 = 0.1e1 / t136;
t167 = 0.1e1 / t135 / Bm;
t257 = Nm(1) * (t138 * coefm(7) * t19 + t143 * coefm(13) * t13 + t138 * coefm(6) * t24 + t150 * coefm(8) * t19) + Nm(6) * (t150 * coefm(14) * t19 + t143 * coefm(22) * t13) + Nm(4) * (t167 * coefm(1) * t13 + t150 * coefm(27) * t24 + t150 * coefm(24) * t19) + Nm(9) * t150 * t19 * coefm(25) + Nm(8) * (t150 * coefm(26) * t19 + t143 * coefm(4) * t13) + Nm(12) * t10 * t167 * t19 * coefm(9) + Nm(10) * (t167 * coefm(2) * t19 + t150 * coefm(10) * t13) + Nm(11) * t167 * t19 * coefm(3) + Nm(7) * t150 * t19 * coefm(15) + Nm(5) * (t167 * coefm(19) * t24 + t167 * coefm(17) * t19) + Nm(3) * (t167 * coefm(18) * t19 + t150 * coefm(5) * t13 + t143 * coefm(23) * t24 + t143 * coefm(20) * t19) + Nm(2) * (t150 * coefm(16) * t19 + t143 * coefm(21) * t13 + t138 * coefm(12) * t24 + t138 * coefm(11) * t19);
% Output
X = t257 + t133;
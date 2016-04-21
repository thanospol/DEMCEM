function X = X_function_g1_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g1_f2

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
t21 = 0.1e1 / t3;
t27 = t9 ^ 2;
t28 = t27 * t12;
t54 = 0.1e1 / t2 / B;
t133 = N(1) * (t5 * coef(3) * t13 + t16 * coef(9) * t19 + t21 * coef(10) * t19 + t16 * coef(11) * t28) + N(2) * (t21 * coef(8) * t19 + t5 * coef(21) * t13 + t16 * coef(7) * t28 + t16 * coef(4) * t19) + N(9) * t21 * t19 * coef(20) + N(5) * (t54 * coef(13) * t28 + t54 * coef(12) * t19) + N(4) * (t54 * coef(25) * t13 + t21 * coef(19) * t28 + t21 * coef(18) * t19) + N(12) * t10 * t54 * t19 * coef(1) + N(6) * (t21 * coef(5) * t19 + t5 * coef(14) * t13) + N(7) * t21 * t19 * coef(6) + N(10) * (t54 * coef(26) * t19 + t21 * coef(2) * t13) + N(11) * t54 * t19 * coef(22) + N(8) * (t21 * coef(17) * t19 + t5 * coef(23) * t13) + N(3) * (t54 * coef(27) * t19 + t21 * coef(24) * t13 + t5 * coef(16) * t28 + t5 * coef(15) * t19);
t135 = Bm ^ 2;
t136 = t135 ^ 2;
t137 = 0.1e1 / t136;
t142 = 0.1e1 / t136 / Bm;
t147 = 0.1e1 / t136 / t135;
t187 = 0.1e1 / t135 / Bm;
t257 = Nm(1) * (t137 * coefm(7) * t19 + t142 * coefm(22) * t13 + t147 * coefm(6) * t28 + t147 * coefm(5) * t19) + Nm(6) * (t137 * coefm(23) * t19 + t142 * coefm(15) * t13) + Nm(2) * (t137 * coefm(27) * t19 + t142 * coefm(14) * t13 + t147 * coefm(25) * t28 + t147 * coefm(26) * t19) + Nm(9) * t137 * t19 * coefm(13) + Nm(5) * (t187 * coefm(19) * t28 + t187 * coefm(18) * t19) + Nm(8) * (t137 * coefm(10) * t19 + t142 * coefm(1) * t13) + Nm(11) * t187 * t19 * coefm(3) + Nm(10) * (t187 * coefm(4) * t19 + t137 * coefm(21) * t13) + Nm(4) * (t187 * coefm(8) * t13 + t137 * coefm(12) * t28 + t137 * coefm(11) * t19) + Nm(12) * t10 * t187 * t19 * coefm(20) + Nm(7) * t137 * t19 * coefm(24) + Nm(3) * (t187 * coefm(17) * t19 + t137 * coefm(2) * t13 + t142 * coefm(9) * t28 + t142 * coefm(16) * t19);
% output
X = t257 + t133;
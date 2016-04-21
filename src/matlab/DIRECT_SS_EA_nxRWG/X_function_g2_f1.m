function X = X_function_g2_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g2_f1

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
t5 = 0.1e1 / t3 / t2;
t8 = cos(Psi);
t9 = t8 ^ 2;
t10 = t9 ^ 2;
t11 = sin(theta);
t12 = t10 * t11;
t14 = 0.1e1 / t3;
t17 = t9 * t11;
t20 = 0.1e1 / t3 / B;
t23 = sin(Psi);
t25 = t9 * t23 * t11;
t49 = 0.1e1 / t2 / B;
t133 = N(1) * (t5 * coef(9) * t12 + t14 * coef(10) * t17 + t20 * coef(6) * t25 + t5 * coef(11) * t17) + N(2) * (t14 * coef(7) * t17 + t20 * coef(27) * t25 + t5 * coef(8) * t12 + t5 * coef(3) * t17) + N(5) * (t49 * coef(17) * t12 + t49 * coef(18) * t17) + N(9) * t14 * t17 * coef(20) + N(6) * (t14 * coef(5) * t17 + t20 * coef(24) * t25) + N(7) * t14 * t17 * coef(4) + N(10) * (t49 * coef(13) * t17 + t14 * coef(2) * t25) + N(4) * (t49 * coef(14) * t25 + t14 * coef(21) * t12 + t14 * coef(25) * t17) + N(3) * (t49 * coef(19) * t17 + t14 * coef(15) * t25 + t20 * coef(22) * t12 + t20 * coef(26) * t17) + N(8) * (t14 * coef(23) * t17 + t20 * coef(16) * t25) + N(11) * t49 * t17 * coef(12) + N(12) * t23 * t49 * t17 * coef(1);
t135 = Bm ^ 2;
t136 = t135 ^ 2;
t137 = 0.1e1 / t136;
t142 = 0.1e1 / t136 / Bm;
t147 = 0.1e1 / t136 / t135;
t182 = 0.1e1 / t135 / Bm;
t257 = Nm(1) * (t137 * coefm(10) * t17 + t142 * coefm(13) * t25 + t147 * coefm(9) * t17 + t147 * coefm(8) * t12) + Nm(6) * (t137 * coefm(16) * t17 + t142 * coefm(26) * t25) + Nm(2) * (t137 * coefm(11) * t17 + t142 * coefm(20) * t25 + t147 * coefm(12) * t12 + t147 * coefm(14) * t17) + Nm(5) * (t182 * coefm(19) * t12 + t182 * coefm(18) * t17) + Nm(9) * t137 * t17 * coefm(24) + Nm(8) * (t137 * coefm(21) * t17 + t142 * coefm(6) * t25) + Nm(11) * t182 * t17 * coefm(4) + Nm(12) * t23 * t182 * t17 * coefm(2) + Nm(4) * (t182 * coefm(5) * t25 + t137 * coefm(22) * t12 + t137 * coefm(23) * t17) + Nm(3) * (t182 * coefm(17) * t17 + t137 * coefm(7) * t25 + t142 * coefm(25) * t12 + t142 * coefm(27) * t17) + Nm(10) * (t182 * coefm(3) * t17 + t137 * coefm(1) * t25) + Nm(7) * t137 * t17 * coefm(15);
% Output
X = t257 + t133;
function X = X_function_g3_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g3_f2

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
t15 = 0.1e1 / t3 / B;
t18 = sin(Psi);
t20 = t9 * t18 * t11;
t24 = t9 * t11;
t26 = 0.1e1 / t3;
t72 = 0.1e1 / t2 / B;
t133 = N(1) * (t5 * coef(18) * t12 + t15 * coef(27) * t20 + t5 * coef(17) * t24 + t26 * coef(19) * t24) + N(2) * (t26 * coef(23) * t24 + t15 * coef(13) * t20 + t5 * coef(22) * t12 + t5 * coef(26) * t24) + N(6) * (t26 * coef(25) * t24 + t15 * coef(8) * t20) + N(7) * t26 * t24 * coef(24) + N(8) * (t26 * coef(7) * t24 + t15 * coef(2) * t20) + N(4) * (t72 * coef(4) * t20 + t26 * coef(11) * t12 + t26 * coef(10) * t24) + N(11) * t72 * t24 * coef(5) + N(3) * (t72 * coef(16) * t24 + t26 * coef(3) * t20 + t15 * coef(9) * t12 + t15 * coef(6) * t24) + N(12) * t18 * t72 * t24 * coef(20) + N(10) * (t72 * coef(1) * t24 + t26 * coef(21) * t20) + N(9) * t26 * t24 * coef(12) + N(5) * (t72 * coef(14) * t12 + t72 * coef(15) * t24);
t135 = Bm ^ 2;
t136 = t135 ^ 2;
t138 = 0.1e1 / t136 / t135;
t143 = 0.1e1 / t136 / Bm;
t150 = 0.1e1 / t136;
t187 = 0.1e1 / t135 / Bm;
t257 = Nm(1) * (t138 * coefm(9) * t24 + t143 * coefm(3) * t20 + t138 * coefm(7) * t12 + t150 * coefm(8) * t24) + Nm(6) * (t150 * coefm(5) * t24 + t143 * coefm(22) * t20) + Nm(2) * (t150 * coefm(6) * t24 + t143 * coefm(24) * t20 + t138 * coefm(2) * t12 + t138 * coefm(1) * t24) + Nm(9) * t150 * t24 * coefm(18) + Nm(3) * (t187 * coefm(25) * t24 + t150 * coefm(12) * t20 + t143 * coefm(19) * t12 + t143 * coefm(23) * t24) + Nm(8) * (t150 * coefm(20) * t24 + t143 * coefm(13) * t20) + Nm(11) * t187 * t24 * coefm(14) + Nm(10) * (t187 * coefm(15) * t24 + t150 * coefm(11) * t20) + Nm(4) * (t187 * coefm(16) * t20 + t150 * coefm(17) * t12 + t150 * coefm(21) * t24) + Nm(12) * t18 * t187 * t24 * coefm(10) + Nm(7) * t150 * t24 * coefm(4) + Nm(5) * (t187 * coefm(27) * t12 + t187 * coefm(26) * t24);
% Output
X = t257 + t133;
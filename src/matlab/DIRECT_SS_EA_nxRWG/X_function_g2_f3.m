function X = X_function_g2_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g2_f3

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
t9 = cos(theta);
t10 = t8 * t9;
t13 = 0.1e1 / t3 / t2;
t16 = t8 ^ 2;
t17 = sin(theta);
t18 = t16 * t17;
t22 = t8 * t17;
t25 = 0.1e1 / t3 / B;
t28 = sin(Psi);
t29 = t10 * t28;
t33 = t7 * t28;
t38 = t8 * t28 * t17;
t45 = t8 * t7;
t66 = 0.1e1 / t2 / B;
t84 = t66 * t8;
t95 = t4 * t8;
t203 = N(1) * (t4 * coef(17) * t10 + t13 * coef(41) * t18 + t13 * coef(42) * t22 + t25 * coef(4) * t29 + t4 * coef(44) * t33 + t25 * coef(1) * t38 + t25 * coef(19) * t7 + t25 * coef(18) * t45) + N(6) * (t4 * coef(11) * t10 + t4 * coef(12) * t22 + t25 * coef(26) * t38 + t4 * coef(7) * t33) + N(10) * (t66 * coef(39) * t10 + t66 * coef(40) * t22 + t4 * coef(16) * t38 + t66 * coef(38) * t33) + N(12) * t17 * t28 * t84 * coef(15) + N(11) * t17 * t84 * coef(37) + N(7) * t17 * t95 * coef(10) + N(3) * (t66 * coef(20) * t10 + t4 * coef(45) * t29 + t4 * coef(46) * t38 + t25 * coef(27) * t18 + t66 * coef(22) * t33 + t4 * coef(9) * t45 + t25 * coef(28) * t22 + t4 * coef(2) * t7) + N(4) * (t66 * coef(35) * t29 + t66 * coef(36) * t38 + t4 * coef(24) * t18 + t66 * coef(14) * t45 + t4 * coef(33) * t22 + t66 * coef(13) * t7) + N(5) * (t66 * coef(23) * t18 + t66 * coef(21) * t22) + N(8) * (t4 * coef(29) * t22 + t4 * coef(30) * t10 + t25 * coef(47) * t38 + t4 * coef(32) * t33) + N(9) * t17 * t95 * coef(31) + N(2) * (t4 * coef(8) * t10 + t25 * coef(25) * t38 + t25 * coef(34) * t29 + t13 * coef(5) * t18 + t4 * coef(3) * t33 + t25 * coef(43) * t45 + t13 * coef(6) * t22 + t25 * coef(48) * t7);
t205 = Bm ^ 2;
t206 = t205 ^ 2;
t207 = 0.1e1 / t206;
t212 = 0.1e1 / t206 / t205;
t217 = 0.1e1 / t206 / Bm;
t256 = 0.1e1 / t205 / Bm;
t257 = t256 * t8;
t347 = t207 * t8;
t392 = Nm(1) * (t207 * coefm(5) * t10 + t212 * coefm(6) * t22 + t217 * coefm(35) * t38 + t217 * coefm(36) * t29 + t212 * coefm(9) * t18 + t207 * coefm(10) * t33 + t217 * coefm(1) * t45 + t217 * coefm(2) * t7) + Nm(6) * (t207 * coefm(37) * t22 + t207 * coefm(38) * t10 + t217 * coefm(14) * t38 + t207 * coefm(40) * t33) + Nm(11) * t17 * t257 * coefm(29) + Nm(12) * t17 * t28 * t257 * coefm(3) + Nm(4) * (t256 * coefm(30) * t38 + t256 * coefm(31) * t29 + t207 * coefm(18) * t18 + t256 * coefm(48) * t45 + t207 * coefm(15) * t22 + t256 * coefm(47) * t7) + Nm(3) * (t256 * coefm(22) * t10 + t207 * coefm(32) * t38 + t207 * coefm(33) * t29 + t217 * coefm(21) * t18 + t256 * coefm(23) * t33 + t207 * coefm(42) * t45 + t217 * coefm(13) * t22 + t207 * coefm(44) * t7) + Nm(10) * (t256 * coefm(26) * t10 + t256 * coefm(27) * t22 + t207 * coefm(4) * t38 + t256 * coefm(28) * t33) + Nm(8) * (t207 * coefm(16) * t10 + t207 * coefm(17) * t22 + t217 * coefm(34) * t38 + t207 * coefm(12) * t33) + Nm(9) * t17 * t347 * coefm(20) + Nm(5) * (t256 * coefm(24) * t18 + t256 * coefm(25) * t22) + Nm(7) * t17 * t347 * coefm(39) + Nm(2) * (t207 * coefm(41) * t10 + t217 * coefm(11) * t29 + t217 * coefm(19) * t38 + t212 * coefm(46) * t18 + t207 * coefm(43) * t33 + t217 * coefm(7) * t45 + t212 * coefm(45) * t22 + t217 * coefm(8) * t7);
% Output
X = t392 + t203;
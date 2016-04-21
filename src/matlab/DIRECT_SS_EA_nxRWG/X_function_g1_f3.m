function X = X_function_g1_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g1_f3

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
t12 = t8 ^ 2;
t13 = t12 * t8;
t16 = 0.1e1 / t3 / t2;
t19 = t12 ^ 2;
t20 = sin(theta);
t21 = t19 * t20;
t23 = 0.1e1 / t3;
t26 = cos(theta);
t27 = t12 * t26;
t31 = sin(Psi);
t32 = t8 * t31;
t36 = t27 * t31;
t40 = t12 * t20;
t45 = t12 * t31 * t20;
t51 = t23 * t12;
t57 = 0.1e1 / t2 / B;
t84 = t57 * t12;
t203 = N(1) * (t5 * coef(27) * t8 + t5 * coef(26) * t13 + t16 * coef(5) * t21 + t23 * coef(2) * t27 + t23 * coef(4) * t32 + t5 * coef(32) * t36 + t16 * coef(6) * t40 + t5 * coef(31) * t45) + N(9) * t20 * t51 * coef(15) + N(5) * (t57 * coef(21) * t21 + t57 * coef(20) * t40) + N(6) * (t23 * coef(35) * t27 + t23 * coef(36) * t40 + t5 * coef(17) * t45 + t23 * coef(30) * t32) + N(12) * t20 * t31 * t84 * coef(29) + N(4) * (t57 * coef(46) * t45 + t57 * coef(48) * t36 + t23 * coef(16) * t21 + t57 * coef(25) * t13 + t23 * coef(10) * t40 + t57 * coef(24) * t8) + N(3) * (t57 * coef(19) * t27 + t23 * coef(40) * t36 + t23 * coef(41) * t45 + t5 * coef(14) * t21 + t57 * coef(18) * t32 + t23 * coef(34) * t13 + t5 * coef(12) * t40 + t23 * coef(23) * t8) + N(10) * (t57 * coef(43) * t27 + t57 * coef(44) * t40 + t23 * coef(28) * t45 + t57 * coef(45) * t32) + N(11) * t20 * t84 * coef(47) + N(8) * (t23 * coef(8) * t27 + t23 * coef(13) * t40 + t5 * coef(42) * t45 + t23 * coef(11) * t32) + N(7) * t20 * t51 * coef(37) + N(2) * (t23 * coef(33) * t27 + t5 * coef(7) * t45 + t5 * coef(9) * t36 + t16 * coef(39) * t21 + t23 * coef(22) * t32 + t5 * coef(3) * t13 + t16 * coef(38) * t40 + t5 * coef(1) * t8);
t205 = Bm ^ 2;
t206 = t205 ^ 2;
t208 = 0.1e1 / t206 / t205;
t212 = 0.1e1 / t206;
t220 = 0.1e1 / t206 / Bm;
t255 = t212 * t12;
t261 = 0.1e1 / t205 / Bm;
t367 = t261 * t12;
t392 = Nm(1) * (t208 * coefm(37) * t40 + t212 * coefm(35) * t27 + t212 * coefm(36) * t32 + t220 * coefm(32) * t13 + t220 * coefm(24) * t36 + t220 * coefm(25) * t45 + t220 * coefm(33) * t8 + t208 * coefm(48) * t21) + Nm(6) * (t212 * coefm(27) * t40 + t212 * coefm(29) * t27 + t220 * coefm(15) * t45 + t212 * coefm(23) * t32) + Nm(7) * t20 * t255 * coefm(28) + Nm(5) * (t261 * coefm(4) * t21 + t261 * coefm(2) * t40) + Nm(9) * t20 * t255 * coefm(12) + Nm(2) * (t212 * coefm(19) * t27 + t220 * coefm(5) * t45 + t220 * coefm(9) * t36 + t208 * coefm(21) * t21 + t212 * coefm(26) * t32 + t220 * coefm(38) * t13 + t208 * coefm(22) * t40 + t220 * coefm(34) * t8) + Nm(3) * (t261 * coefm(3) * t27 + t212 * coefm(45) * t36 + t212 * coefm(46) * t45 + t220 * coefm(13) * t21 + t261 * coefm(1) * t32 + t212 * coefm(18) * t13 + t220 * coefm(14) * t40 + t212 * coefm(20) * t8) + Nm(4) * (t261 * coefm(42) * t36 + t261 * coefm(43) * t45 + t212 * coefm(7) * t21 + t261 * coefm(31) * t13 + t212 * coefm(10) * t40 + t261 * coefm(30) * t8) + Nm(10) * (t261 * coefm(41) * t27 + t261 * coefm(44) * t40 + t212 * coefm(17) * t45 + t261 * coefm(39) * t32) + Nm(11) * t20 * t367 * coefm(40) + Nm(8) * (t212 * coefm(6) * t40 + t212 * coefm(8) * t27 + t220 * coefm(47) * t45 + t212 * coefm(11) * t32) + Nm(12) * t20 * t31 * t367 * coefm(16);
% Output
X = t392 + t203;
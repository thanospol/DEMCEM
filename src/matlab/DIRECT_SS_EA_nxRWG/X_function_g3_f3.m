function X = X_function_g3_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g3_f3

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
t17 = cos(theta);
t18 = t9 * t17;
t19 = t18 * t10;
t21 = 0.1e1 / t3;
t26 = 0.1e1 / t3 / t2;
t29 = t9 ^ 2;
t30 = t29 * t12;
t34 = t9 * t12;
t38 = t8 * t10;
t45 = t9 * t8;
t51 = 0.1e1 / t2 / B;
t69 = t51 * t9;
t95 = t21 * t9;
t203 = N(1) * (t5 * coef(30) * t13 + t5 * coef(31) * t19 + t21 * coef(44) * t18 + t26 * coef(45) * t30 + t26 * coef(46) * t34 + t21 * coef(48) * t38 + t5 * coef(33) * t8 + t5 * coef(32) * t45) + N(10) * (t51 * coef(38) * t34 + t51 * coef(39) * t18 + t21 * coef(2) * t13 + t51 * coef(34) * t38) + N(12) * t12 * t10 * t69 * coef(1) + N(11) * t12 * t69 * coef(37) + N(6) * (t21 * coef(26) * t34 + t21 * coef(27) * t18 + t5 * coef(8) * t13 + t21 * coef(25) * t38) + N(7) * t12 * t95 * coef(20) + N(3) * (t51 * coef(4) * t18 + t21 * coef(41) * t19 + t21 * coef(42) * t13 + t5 * coef(12) * t30 + t51 * coef(5) * t38 + t21 * coef(28) * t45 + t5 * coef(10) * t34 + t21 * coef(24) * t8) + N(4) * (t51 * coef(35) * t13 + t51 * coef(36) * t19 + t21 * coef(14) * t30 + t51 * coef(19) * t45 + t21 * coef(7) * t34 + t51 * coef(18) * t8) + N(8) * (t21 * coef(15) * t18 + t21 * coef(17) * t34 + t5 * coef(40) * t13 + t21 * coef(9) * t38) + N(9) * t12 * t95 * coef(16) + N(5) * (t51 * coef(6) * t30 + t51 * coef(3) * t34) + N(2) * (t21 * coef(29) * t18 + t5 * coef(11) * t19 + t5 * coef(13) * t13 + t26 * coef(21) * t30 + t21 * coef(23) * t38 + t5 * coef(43) * t45 + t26 * coef(22) * t34 + t5 * coef(47) * t8);
t205 = Bm ^ 2;
t206 = t205 ^ 2;
t208 = 0.1e1 / t206 / Bm;
t212 = 0.1e1 / t206;
t220 = 0.1e1 / t206 / t205;
t255 = t212 * t9;
t261 = 0.1e1 / t205 / Bm;
t278 = t261 * t9;
t392 = Nm(1) * (t208 * coefm(8) * t19 + t212 * coefm(35) * t38 + t208 * coefm(18) * t45 + t220 * coefm(37) * t30 + t220 * coefm(38) * t34 + t212 * coefm(39) * t18 + t208 * coefm(9) * t13 + t208 * coefm(17) * t8) + Nm(6) * (t212 * coefm(5) * t34 + t212 * coefm(7) * t18 + t208 * coefm(24) * t13 + t212 * coefm(10) * t38) + Nm(7) * t12 * t255 * coefm(6) + Nm(10) * (t261 * coefm(40) * t34 + t261 * coefm(43) * t18 + t212 * coefm(2) * t13 + t261 * coefm(41) * t38) + Nm(11) * t12 * t278 * coefm(42) + Nm(12) * t12 * t10 * t278 * coefm(1) + Nm(2) * (t212 * coefm(3) * t18 + t208 * coefm(25) * t13 + t208 * coefm(26) * t19 + t220 * coefm(11) * t30 + t212 * coefm(13) * t38 + t208 * coefm(34) * t45 + t220 * coefm(12) * t34 + t208 * coefm(36) * t8) + Nm(4) * (t261 * coefm(44) * t19 + t261 * coefm(45) * t13 + t212 * coefm(28) * t30 + t261 * coefm(16) * t45 + t212 * coefm(27) * t34 + t261 * coefm(15) * t8) + Nm(3) * (t261 * coefm(21) * t18 + t212 * coefm(47) * t19 + t212 * coefm(48) * t13 + t208 * coefm(29) * t30 + t261 * coefm(19) * t38 + t212 * coefm(4) * t45 + t208 * coefm(33) * t34 + t212 * coefm(14) * t8) + Nm(8) * (t212 * coefm(30) * t34 + t212 * coefm(32) * t18 + t208 * coefm(46) * t13 + t212 * coefm(23) * t38) + Nm(9) * t12 * t255 * coefm(31) + Nm(5) * (t261 * coefm(22) * t30 + t261 * coefm(20) * t34);
% Output
X = t392 + t203;
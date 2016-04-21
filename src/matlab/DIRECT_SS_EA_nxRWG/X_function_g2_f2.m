function X = X_function_g2_f2(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g2_f2

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
t10 = sin(theta);
t11 = t9 * t10;
t13 = 0.1e1 / t3;
t19 = t9 ^ 2;
t20 = t19 * t10;
t38 = 0.1e1 / t2 / B;
t66 = 0.1e1 / t3 / B;
t76 = sin(Psi);
t110 = N(1) * (t5 * coef(9) * t11 + t13 * coef(3) * t11 + t5 * coef(8) * t20) + N(2) * (t13 * coef(4) * t11 + t5 * coef(6) * t20 + t5 * coef(7) * t11) + N(5) * (t38 * coef(12) * t20 + t38 * coef(13) * t11) + N(4) * (t13 * coef(20) * t20 + t13 * coef(19) * t11) + N(9) * t13 * t11 * coef(16) + N(3) * (t38 * coef(14) * t11 + t66 * coef(18) * t20 + t66 * coef(17) * t11) + N(10) * t76 * t13 * t11 * coef(1) + N(12) * t76 * t38 * t11 * coef(2) + N(11) * t38 * t11 * coef(10) + N(8) * t76 * t66 * t11 * coef(11) + N(6) * t76 * t66 * t11 * coef(15) + N(7) * t13 * t11 * coef(5);
t112 = Bm ^ 2;
t113 = t112 ^ 2;
t115 = 0.1e1 / t113 / t112;
t119 = 0.1e1 / t113;
t131 = 0.1e1 / t113 / Bm;
t150 = 0.1e1 / t112 / Bm;
t213 = Nm(1) * (t115 * coefm(5) * t20 + t119 * coefm(6) * t11 + t115 * coefm(4) * t11) + Nm(6) * t76 * t131 * t11 * coefm(11) + Nm(2) * (t119 * coefm(1) * t11 + t115 * coefm(17) * t20 + t115 * coefm(18) * t11) + Nm(3) * (t150 * coefm(8) * t11 + t131 * coefm(14) * t20 + t131 * coefm(13) * t11) + Nm(5) * (t150 * coefm(10) * t20 + t150 * coefm(9) * t11) + Nm(7) * t119 * t11 * coefm(7) + Nm(9) * t119 * t11 * coefm(12) + Nm(4) * (t119 * coefm(16) * t20 + t119 * coefm(15) * t11) + Nm(12) * t76 * t150 * t11 * coefm(3) + Nm(10) * t76 * t119 * t11 * coefm(2) + Nm(11) * t150 * t11 * coefm(19) + Nm(8) * t76 * t131 * t11 * coefm(20);
% Output
X = t213 + t110;
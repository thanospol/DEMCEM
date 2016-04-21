function X = X_function_g1_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_g1_f1

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
t16 = t9 * t11;
t18 = 0.1e1 / t3;
t38 = 0.1e1 / t2 / B;
t43 = 0.1e1 / t3 / B;
t76 = sin(Psi);
t110 = N(1) * (t5 * coef(9) * t12 + t5 * coef(7) * t16 + t18 * coef(8) * t16) + N(2) * (t18 * coef(5) * t16 + t5 * coef(3) * t12 + t5 * coef(4) * t16) + N(3) * (t38 * coef(14) * t16 + t43 * coef(18) * t12 + t43 * coef(17) * t16) + N(5) * (t38 * coef(12) * t12 + t38 * coef(13) * t16) + N(4) * (t18 * coef(20) * t12 + t18 * coef(19) * t16) + N(9) * t18 * t16 * coef(16) + N(6) * t76 * t43 * t16 * coef(15) + N(7) * t18 * t16 * coef(6) + N(10) * t76 * t18 * t16 * coef(2) + N(11) * t38 * t16 * coef(10) + N(8) * t76 * t43 * t16 * coef(11) + N(12) * t76 * t38 * t16 * coef(1);
t112 = Bm ^ 2;
t113 = t112 ^ 2;
t115 = 0.1e1 / t113 / t112;
t122 = 0.1e1 / t113;
t131 = 0.1e1 / t113 / Bm;
t155 = 0.1e1 / t112 / Bm;
% Output
t213 = Nm(1) * (t115 * coefm(5) * t12 + t115 * coefm(19) * t16 + t122 * coefm(20) * t16) + Nm(6) * t76 * t131 * t16 * coefm(12) + Nm(2) * (t122 * coefm(6) * t16 + t115 * coefm(8) * t12 + t115 * coefm(7) * t16) + Nm(7) * t122 * t16 * coefm(9) + Nm(5) * (t155 * coefm(11) * t12 + t155 * coefm(10) * t16) + Nm(9) * t122 * t16 * coefm(13) + Nm(3) * (t155 * coefm(18) * t16 + t131 * coefm(15) * t12 + t131 * coefm(14) * t16) + Nm(4) * (t122 * coefm(17) * t12 + t122 * coefm(16) * t16) + Nm(12) * t76 * t155 * t16 * coefm(4) + Nm(11) * t155 * t16 * coefm(2) + Nm(8) * t76 * t131 * t16 * coefm(1) + Nm(10) * t76 * t122 * t16 * coefm(3);
X = t213 + t110;
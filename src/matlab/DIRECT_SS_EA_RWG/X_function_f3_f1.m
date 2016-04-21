function X = X_function_f3_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f3_f1

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
t9 = sin(theta);
t10 = t8 * t9;
t13 = 0.1e1 / t3 / B;
t16 = sin(Psi);
t18 = t8 * t16 * t9;
t34 = 0.1e1 / t2 / B;
t64 = Bm ^ 2;
t65 = t64 ^ 2;
t67 = 0.1e1 / t65 / Bm;
t71 = 0.1e1 / t65;
t89 = 0.1e1 / t64 / Bm;
% Output
X = N(1) * (t4 * coef(1) * t10 + t13 * coef(5) * t18) + N(8) * t4 * t10 * coef(9) + N(6) * t4 * t10 * coef(6) + N(10) * t34 * t10 * coef(2) + N(3) * (t34 * coef(8) * t10 + t4 * coef(4) * t18) + N(4) * t16 * t34 * t10 * coef(3) + N(2) * (t4 * coef(7) * t10 + t13 * coef(10) * t18) + Nm(1) * (t67 * coefm(2) * t18 + t71 * coefm(4) * t10) + Nm(6) * t71 * t10 * coefm(1) + Nm(8) * t71 * t10 * coefm(8) + Nm(3) * (t89 * coefm(10) * t10 + t71 * coefm(7) * t18) + Nm(4) * t16 * t89 * t10 * coefm(5) + Nm(2) * (t71 * coefm(3) * t10 + t67 * coefm(9) * t18) + Nm(10) * t89 * t10 * coefm(6);
function X = X_function_f3_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f3_f3

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
t13 = 0.1e1 / t3 / B;
t16 = sin(Psi);
t17 = t10 * t16;
t23 = 0.1e1 / t2 / B;
t35 = sin(theta);
t36 = t8 * t35;
t40 = t7 * t16;
t79 = Bm ^ 2;
t80 = t79 ^ 2;
t81 = 0.1e1 / t80;
t86 = 0.1e1 / t80 / Bm;
t103 = 0.1e1 / t79 / Bm;
% Output
X = N(1) * (t4 * coef(1) * t10 + t13 * coef(5) * t17) + N(3) * (t23 * coef(6) * t10 + t4 * coef(12) * t17) + N(8) * (t4 * coef(7) * t36 + t4 * coef(8) * t40) + N(10) * (t23 * coef(13) * t36 + t23 * coef(11) * t40) + N(4) * t16 * t9 * t23 * t8 * coef(10) + N(6) * (t4 * coef(3) * t36 + t4 * coef(2) * t40) + N(2) * (t4 * coef(4) * t10 + t13 * coef(9) * t17) + Nm(1) * (t81 * coefm(9) * t10 + t86 * coefm(12) * t17) + Nm(6) * (t81 * coefm(10) * t36 + t81 * coefm(11) * t40) + Nm(10) * (t103 * coefm(4) * t36 + t103 * coefm(1) * t40) + Nm(3) * (t103 * coefm(8) * t10 + t81 * coefm(2) * t17) + Nm(4) * t16 * t9 * t103 * t8 * coefm(3) + Nm(8) * (t81 * coefm(5) * t36 + t81 * coefm(6) * t40) + Nm(2) * (t81 * coefm(13) * t10 + t86 * coefm(7) * t17);
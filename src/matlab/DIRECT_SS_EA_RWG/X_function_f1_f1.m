function X = X_function_f1_f1(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f1_f1

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
t9 = cos(Psi);
t10 = t9 ^ 2;
t11 = sin(Psi);
t13 = sin(theta);
t14 = t10 * t11 * t13;
t21 = t11 * t13;
t39 = Bm ^ 2;
t40 = t39 ^ 2;
t42 = 0.1e1 / t40 / Bm;
% Output
X = N(1) * t5 * coef(1) * t14 + N(4) / t2 / B * t10 * t21 * coef(3) + N(3) / t3 * t10 * t21 * coef(2) + N(2) * t5 * t10 * t21 * coef(4) + Nm(1) * t42 * coefm(1) * t14 + Nm(2) * t42 * t10 * t21 * coefm(4) + Nm(3) / t40 * t10 * t21 * coefm(2) + Nm(4) / t39 / Bm * t10 * t21 * coefm(3);

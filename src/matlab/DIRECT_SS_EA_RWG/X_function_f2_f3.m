function X = X_function_f2_f3(theta, Psi, B, Bm, coef, coefm, N, Nm)
%% X_function_f2_f3

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
t8 = sin(Psi);
t9 = t7 * t8;
t12 = 0.1e1 / t3 / B;
t15 = t7 ^ 2;
t17 = sin(theta);
t18 = t15 * t8 * t17;
t22 = cos(theta);
t24 = t15 * t22 * t8;
t40 = 0.1e1 / t2 / B;
t79 = Bm ^ 2;
t80 = t79 ^ 2;
t81 = 0.1e1 / t80;
t86 = 0.1e1 / t80 / Bm;
t107 = 0.1e1 / t79 / Bm;
%
X = N(1) * (t4 * coef(1) * t9 + t12 * coef(10) * t18 + t12 * coef(8) * t24) + N(8) * t4 * t9 * coef(12) + N(6) * t4 * t9 * coef(9) + N(10) * t40 * t9 * coef(3) + N(4) * (t40 * coef(2) * t18 + t40 * coef(4) * t24) + N(3) * (t4 * coef(5) * t24 + t4 * coef(6) * t18 + t40 * coef(11) * t9) + N(2) * (t12 * coef(13) * t24 + t12 * coef(14) * t18 + t4 * coef(7) * t9) + Nm(1) * (t81 * coefm(6) * t9 + t86 * coefm(7) * t18 + t86 * coefm(8) * t24) + Nm(6) * t81 * t9 * coefm(10) + Nm(8) * t81 * t9 * coefm(11) + Nm(4) * (t107 * coefm(1) * t18 + t107 * coefm(3) * t24) + Nm(3) * (t81 * coefm(4) * t24 + t81 * coefm(5) * t18 + t107 * coefm(14) * t9) + Nm(2) * (t86 * coefm(12) * t24 + t86 * coefm(13) * t18 + t81 * coefm(9) * t9) + Nm(10) * t107 * t9 * coefm(2);

function X = Omega_function_g1_f1(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g1_f1

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
% theta_p, theta_q, Psi, GAMMA, coef, K

% OUTPUT DATA
% X   
%

global ko;
%
j = sqrt(-1);
%
t2 = sin(Psi);
t3 = t2 ^ 2;
t4 = cos(Psi);
t5 = t4 ^ 2;
t7 = t3 * t5 * t4;
t8 = GAMMA ^ 2;
t9 = 0.1e1 / t8;
t11 = t7 * t9 * j;
t12 = cos(theta_q);
t13 = ko * t12;
t14 = coef(5);
t15 = cos(theta_p);
t16 = t15 ^ 2;
t20 = t7 * t9;
t21 = j * ko;
t22 = sin(theta_q);
t23 = coef(6);
t24 = t22 * t23;
t27 = ko * t22;
t31 = coef(4);
t35 = sin(theta_p);
t36 = t35 * t15;
t37 = coef(1);
t41 = coef(2);
t45 = coef(3);
t49 = t12 * t14;
t56 = 0.1e1 / t8 / GAMMA;
t57 = t7 * t56;
% Output
X = K(4) * (t11 * t13 * t14 * t16 - t20 * t21 * t24 + t11 * t27 * t23 * t16 - t11 * t27 * t16 * t31 - t11 * t13 * t36 * t37 - t11 * t27 * t36 * t41 - t11 * t13 * t16 * t45 - t20 * t21 * t49) + K(3) * (-t57 * t22 * t16 * t31 + t57 * t49 * t16 - t7 * t56 * t22 * t23 + t57 * t24 * t16 - t7 * t56 * t12 * t14 - t57 * t12 * t35 * t15 * t37 - t57 * t22 * t35 * t15 * t41 - t57 * t12 * t16 * t45);
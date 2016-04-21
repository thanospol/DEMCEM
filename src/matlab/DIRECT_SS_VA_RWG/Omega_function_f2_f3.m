function X = Omega_function_f2_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_f2_f3

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
t6 = t3 * t5;
t7 = GAMMA ^ 2;
t9 = 0.1e1 / t7 / GAMMA;
t10 = t6 * t9;
t11 = cos(theta_q);
t12 = cos(theta_p);
t14 = coef(5);
t17 = sin(theta_q);
t19 = coef(6);
t22 = sin(theta_p);
t24 = coef(7);
t28 = coef(4);
t31 = t2 * t5;
t32 = 0.1e1 / t7;
t34 = j * ko;
t35 = coef(2);
t39 = t3 * t4;
t40 = t39 * t32;
t41 = coef(3);
t45 = coef(1);
t53 = t6 * t32 * j;
t58 = ko * t11;
% Output
X = K(2) * (-t10 * t11 * t12 * t14 - t10 * t17 * t12 * t19 - t10 * t11 * t22 * t24 - t10 * t22 * t17 * t28 - t31 * t32 * t34 * t22 * t35 - t40 * t34 * t11 * t41 - t40 * t34 * t17 * t45) + K(3) * (-t53 * ko * t22 * t17 * t28 - t53 * t58 * t22 * t24 - t53 * t58 * t12 * t14 - t53 * ko * t17 * t12 * t19) + K(1) * (-t31 * t9 * t22 * t35 - t39 * t9 * t11 * t41 - t39 * t9 * t17 * t45);
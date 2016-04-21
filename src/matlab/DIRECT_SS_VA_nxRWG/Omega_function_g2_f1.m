function X = Omega_function_g2_f1(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g2_f1

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
t11 = sin(theta_q);
t12 = cos(theta_p);
t14 = coef(7);
t17 = sin(theta_p);
t18 = t17 * t11;
t19 = coef(8);
t22 = cos(theta_q);
t23 = t22 * t17;
t24 = coef(5);
t28 = coef(6);
t35 = t3 * t5 * t4;
t36 = 0.1e1 / t7;
t37 = t36 * j;
t38 = t35 * t37;
t39 = ko * t22;
t40 = coef(1);
t41 = t12 ^ 2;
t45 = t35 * t36;
t46 = j * ko;
t47 = coef(2);
t48 = t11 * t47;
t51 = ko * t11;
t55 = coef(3);
t59 = coef(4);
t63 = t17 * t12;
t64 = coef(10);
t68 = t22 * t40;
t71 = coef(9);
t78 = t35 * t9;
t98 = t6 * t37;
t115 = -t78 * t11 * t41 * t59 - t78 * t18 * t12 * t64 - t35 * t9 * t22 * t40 - t35 * t9 * t11 * t47 + t78 * t48 * t41 - t78 * t22 * t41 * t55 + t78 * t68 * t41 - t98 * t39 * t12 * t28 - t98 * t51 * t12 * t14 - t98 * ko * t17 * t11 * t19 - t98 * t39 * t17 * t24 - t78 * t23 * t12 * t71;
% Output
X = K(2) * (-t10 * t11 * t12 * t14 - t10 * t18 * t19 - t10 * t23 * t24 - t10 * t22 * t12 * t28) + K(4) * (t38 * t39 * t40 * t41 - t45 * t46 * t48 + t38 * t51 * t47 * t41 - t38 * t39 * t41 * t55 - t38 * t51 * t41 * t59 - t38 * t51 * t63 * t64 - t45 * t46 * t68 - t38 * t39 * t63 * t71) + K(3) * t115;
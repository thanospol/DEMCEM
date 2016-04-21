function X = Omega_function_g3_f1(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g3_f1

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
t11 = sin(theta_p);
t12 = sin(theta_q);
t13 = t11 * t12;
t14 = coef(1);
t17 = cos(theta_p);
t18 = t17 * t12;
t19 = coef(6);
t22 = cos(theta_q);
t23 = t22 * t11;
t24 = coef(7);
t28 = coef(10);
t35 = t3 * t5 * t4;
t36 = t35 * t9;
t37 = coef(2);
t41 = coef(3);
t46 = coef(5);
t49 = t22 * t46;
t50 = t17 ^ 2;
t53 = coef(4);
t54 = t12 * t53;
t58 = coef(8);
t62 = coef(9);
t68 = 0.1e1 / t7;
t69 = t68 * j;
t70 = t6 * t69;
t75 = ko * t17;
t79 = ko * t22;
t86 = -t36 * t18 * t11 * t37 - t36 * t23 * t17 * t41 - t35 * t9 * t22 * t46 + t36 * t49 * t50 + t36 * t54 * t50 - t36 * t22 * t50 * t58 - t36 * t12 * t50 * t62 - t35 * t9 * t12 * t53 - t70 * ko * t11 * t12 * t14 - t70 * t75 * t12 * t19 - t70 * t79 * t11 * t24 - t70 * t79 * t17 * t28;
t89 = t35 * t69;
t97 = t35 * t68;
t98 = j * ko;
t104 = ko * t12;
% Output
X = K(2) * (-t10 * t13 * t14 - t10 * t18 * t19 - t10 * t23 * t24 - t10 * t22 * t17 * t28) + K(3) * t86 + K(4) * (-t89 * t75 * t13 * t37 - t89 * t79 * t11 * t17 * t41 - t97 * t98 * t49 + t89 * t79 * t46 * t50 + t89 * t104 * t53 * t50 - t89 * t79 * t50 * t58 - t89 * t104 * t50 * t62 - t97 * t98 * t54);
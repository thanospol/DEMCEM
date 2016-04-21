function X = Omega_function_g3_f2(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g3_f2

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
t3 = cos(Psi);
t4 = t3 ^ 2;
t5 = t2 * t4;
t6 = GAMMA ^ 2;
t8 = 0.1e1 / t6 / GAMMA;
t9 = cos(theta_p);
t11 = coef(5);
t14 = sin(theta_p);
t16 = coef(4);
t19 = t2 ^ 2;
t20 = t19 * t3;
t21 = sin(theta_q);
t22 = t8 * t21;
t23 = coef(6);
t29 = t19 * t4;
t30 = t29 * t8;
t31 = t14 * t21;
t32 = coef(1);
t35 = cos(theta_q);
t36 = t35 * t14;
t37 = coef(13);
t41 = coef(11);
t44 = t4 * t3;
t45 = t2 * t44;
t47 = t14 * t9;
t48 = coef(14);
t51 = coef(16);
t52 = t8 * t51;
t53 = t9 ^ 2;
t57 = coef(15);
t61 = coef(12);
t64 = 0.1e1 / t6;
t65 = t5 * t64;
t66 = j * ko;
t78 = -t30 * t31 * t32 - t30 * t36 * t37 - t30 * t35 * t9 * t41 - t45 * t8 * t47 * t48 + t45 * t52 * t53 - t45 * t8 * t53 * t57 - t30 * t21 * t9 * t61 - t65 * t66 * t14 * t16 - t45 * t52 - t20 * t64 * t66 * t21 * t23 - t65 * t66 * t9 * t11;
t81 = t19 * t44;
t82 = t81 * t8;
t84 = coef(7);
t88 = coef(8);
t91 = coef(2);
t95 = t64 * j;
t96 = t29 * t95;
t97 = ko * t14;
t101 = coef(10);
t102 = t21 * t101;
t105 = t45 * t64;
t109 = ko * t35;
t121 = coef(3);
t126 = coef(9);
t129 = t35 * t126;
t132 = ko * t21;
t142 = -t82 * t35 * t53 * t84 - t82 * t21 * t53 * t88 - t82 * t36 * t9 * t91 - t96 * t97 * t21 * t32 + t82 * t102 * t53 - t105 * t66 * t53 * t57 - t96 * t109 * t14 * t37 - t81 * t22 * t101 - t96 * t109 * t9 * t41 + t105 * t66 * t51 * t53 - t82 * t31 * t9 * t121 - t81 * t8 * t35 * t126 + t82 * t129 * t53 - t96 * t132 * t9 * t61 - t45 * t95 * t97 * t9 * t48 - t105 * t66 * t51;
t145 = t81 * t95;
t164 = t81 * t64;
% Output
X = K(1) * (-t5 * t8 * t9 * t11 - t5 * t8 * t14 * t16 - t20 * t22 * t23) + K(2) * t78 + K(3) * t142 + K(4) * (-t145 * t109 * t53 * t84 - t145 * t132 * t53 * t88 - t145 * t132 * t47 * t121 + t145 * t109 * t126 * t53 + t145 * t132 * t101 * t53 - t145 * t109 * t47 * t91 - t164 * t66 * t129 - t164 * t66 * t102);
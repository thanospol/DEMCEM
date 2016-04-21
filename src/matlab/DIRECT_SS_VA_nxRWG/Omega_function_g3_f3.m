function X = Omega_function_g3_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g3_f3

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
t5 = t4 * t3;
t6 = t2 * t5;
t7 = GAMMA ^ 2;
t8 = 0.1e1 / t7;
t9 = t6 * t8;
t10 = j * ko;
t11 = coef(6);
t12 = cos(theta_p);
t13 = t12 ^ 2;
t17 = t2 ^ 2;
t18 = t17 * t4;
t19 = t8 * j;
t20 = t18 * t19;
t21 = sin(theta_q);
t22 = ko * t21;
t23 = coef(15);
t27 = sin(theta_p);
t28 = ko * t27;
t29 = coef(7);
t35 = cos(theta_q);
t36 = ko * t35;
t37 = coef(14);
t41 = t17 * t5;
t43 = 0.1e1 / t7 / GAMMA;
t44 = t43 * t21;
t45 = coef(13);
t48 = coef(5);
t52 = t41 * t43;
t53 = t35 * t27;
t54 = coef(8);
t59 = coef(11);
t62 = t43 * t35;
t63 = coef(12);
t66 = t21 * t27;
t67 = coef(9);
t71 = coef(16);
t76 = coef(10);
t80 = coef(17);
t84 = t21 * t45;
t87 = t35 * t63;
t90 = t9 * t10 * t11 * t13 - t20 * t22 * t12 * t23 - t20 * t28 * t21 * t29 - t9 * t10 * t11 - t20 * t36 * t12 * t37 - t41 * t44 * t45 - t9 * t10 * t13 * t48 - t52 * t53 * t12 * t54 - t52 * t21 * t13 * t59 - t41 * t62 * t63 - t52 * t66 * t12 * t67 - t20 * t36 * t27 * t71 - t52 * t35 * t13 * t76 - t6 * t19 * t28 * t12 * t80 + t52 * t84 * t13 + t52 * t87 * t13;
t93 = t41 * t8;
t96 = t41 * t19;
t111 = t27 * t12;
t121 = t2 * t4;
t122 = t121 * t8;
t123 = coef(4);
t127 = t17 * t3;
t128 = t127 * t8;
t129 = coef(2);
t133 = coef(1);
t137 = t18 * t43;
t143 = coef(3);
t155 = t43 * t11;
t162 = -t122 * t10 * t12 * t123 - t128 * t10 * t35 * t129 - t122 * t10 * t27 * t133 - t137 * t21 * t12 * t23 - t137 * t66 * t29 - t128 * t10 * t21 * t143 - t137 * t35 * t12 * t37 - t6 * t43 * t111 * t80 - t137 * t53 * t71 + t6 * t155 * t13 - t6 * t155 - t6 * t43 * t13 * t48;
% Output
X = K(3) * t90 + K(4) * (-t93 * t10 * t84 + t96 * t36 * t63 * t13 + t96 * t22 * t45 * t13 - t96 * t36 * t13 * t76 - t93 * t10 * t87 - t96 * t22 * t13 * t59 - t96 * t22 * t111 * t67 - t96 * t36 * t111 * t54) + K(2) * t162 + K(1) * (-t127 * t44 * t143 - t121 * t43 * t27 * t133 - t127 * t62 * t129 - t121 * t43 * t12 * t123);
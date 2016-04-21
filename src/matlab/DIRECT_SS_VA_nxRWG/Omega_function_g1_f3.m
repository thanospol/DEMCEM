function X = Omega_function_g1_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g1_f3

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
t9 = 0.1e1 / t7 / GAMMA;
t10 = cos(theta_p);
t11 = t10 ^ 2;
t13 = coef(13);
t16 = t2 ^ 2;
t17 = t16 * t4;
t18 = t17 * t9;
t19 = sin(theta_q);
t20 = sin(theta_p);
t21 = t19 * t20;
t22 = coef(7);
t25 = cos(theta_q);
t27 = coef(8);
t31 = coef(9);
t34 = t25 * t20;
t35 = coef(6);
t39 = t20 * t10;
t40 = coef(12);
t43 = coef(1);
t44 = t9 * t43;
t51 = t16 * t5;
t53 = coef(5);
t57 = coef(4);
t60 = t51 * t9;
t61 = coef(11);
t65 = 0.1e1 / t7;
t66 = t65 * j;
t67 = t17 * t66;
t68 = ko * t19;
t72 = ko * t25;
t79 = coef(10);
t84 = coef(3);
t93 = coef(2);
t96 = t25 * t57;
t102 = t6 * t65;
t103 = j * ko;
t107 = t19 * t53;
t115 = -t51 * t9 * t19 * t53 - t51 * t9 * t25 * t57 - t60 * t21 * t10 * t61 - t67 * t68 * t20 * t22 - t67 * t72 * t10 * t27 - t67 * t68 * t10 * t31 - t60 * t34 * t10 * t79 - t60 * t19 * t11 * t84 - t6 * t66 * ko * t20 * t10 * t40 - t60 * t25 * t11 * t93 + t60 * t96 * t11 - t67 * t72 * t20 * t35 - t102 * t103 * t11 * t13 + t60 * t107 * t11 + t102 * t103 * t43 * t11 - t102 * t103 * t43;
t118 = t51 * t66;
t137 = t51 * t65;
% Output
X = K(2) * (-t6 * t9 * t11 * t13 - t18 * t21 * t22 - t18 * t25 * t10 * t27 - t18 * t19 * t10 * t31 - t18 * t34 * t35 - t6 * t9 * t39 * t40 - t6 * t44 + t6 * t44 * t11) + K(3) * t115 + K(4) * (-t118 * t72 * t39 * t79 - t118 * t68 * t39 * t61 - t118 * t68 * t11 * t84 + t118 * t72 * t57 * t11 + t118 * t68 * t53 * t11 - t118 * t72 * t11 * t93 - t137 * t103 * t96 - t137 * t103 * t107);
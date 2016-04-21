function X = Omega_function_g2_f2(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g2_f2

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
t5 = t3 * t4;
t6 = GAMMA ^ 2;
t8 = 0.1e1 / t6 / GAMMA;
t9 = sin(theta_q);
t10 = t8 * t9;
t11 = coef(11);
t14 = t4 ^ 2;
t15 = t2 * t14;
t16 = cos(theta_p);
t18 = coef(12);
t21 = sin(theta_p);
t23 = coef(13);
t29 = t3 * t14;
t30 = 0.1e1 / t6;
t31 = t30 * j;
t32 = t29 * t31;
t33 = ko * t21;
t34 = coef(8);
t38 = t14 * t4;
t39 = t2 * t38;
t41 = coef(16);
t45 = t39 * t30;
t46 = j * ko;
t47 = coef(15);
t50 = t16 ^ 2;
t51 = coef(14);
t55 = t3 * t38;
t56 = t55 * t8;
t57 = cos(theta_q);
t58 = t57 * t21;
t59 = coef(9);
t67 = coef(1);
t71 = coef(2);
t74 = ko * t9;
t75 = coef(7);
t79 = coef(4);
t82 = ko * t57;
t83 = coef(5);
t87 = coef(6);
t91 = t9 * t21;
t92 = coef(10);
t97 = coef(3);
t100 = t9 * t79;
t103 = t57 * t97;
t106 = -t32 * t33 * t9 * t34 - t39 * t31 * t33 * t16 * t41 - t45 * t46 * t47 - t45 * t46 * t50 * t51 - t56 * t58 * t16 * t59 + t45 * t46 * t47 * t50 - t56 * t57 * t50 * t67 - t56 * t9 * t50 * t71 - t32 * t74 * t16 * t75 - t55 * t10 * t79 - t32 * t82 * t21 * t83 - t32 * t82 * t16 * t87 - t56 * t91 * t16 * t92 - t55 * t8 * t57 * t97 + t56 * t100 * t50 + t56 * t103 * t50;
t109 = t55 * t31;
t110 = t21 * t16;
t123 = t55 * t30;
t137 = t15 * t30;
t147 = t29 * t8;
t159 = t8 * t47;
t169 = -t137 * t46 * t16 * t18 - t39 * t8 * t110 * t41 - t39 * t8 * t50 * t51 - t147 * t57 * t16 * t87 - t147 * t58 * t83 - t147 * t9 * t16 * t75 - t137 * t46 * t21 * t23 + t39 * t159 * t50 - t5 * t30 * t46 * t9 * t11 - t147 * t91 * t34 - t39 * t159;
% Output
X = K(1) * (-t5 * t10 * t11 - t15 * t8 * t16 * t18 - t15 * t8 * t21 * t23) + K(3) * t106 + K(4) * (-t109 * t82 * t110 * t59 + t109 * t74 * t79 * t50 + t109 * t82 * t97 * t50 - t109 * t74 * t50 * t71 - t123 * t46 * t100 - t109 * t82 * t50 * t67 - t109 * t74 * t110 * t92 - t123 * t46 * t103) + K(2) * t169;
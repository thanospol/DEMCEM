function X = Omega_function_g2_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g2_f3

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
t7 = 0.1e1 / t6;
t8 = t5 * t7;
t9 = j * ko;
t10 = sin(theta_p);
t11 = coef(7);
t15 = cos(theta_p);
t16 = coef(8);
t20 = t4 * t3;
t21 = t2 * t20;
t23 = 0.1e1 / t6 / GAMMA;
t24 = coef(10);
t25 = t23 * t24;
t27 = t15 ^ 2;
t31 = coef(9);
t34 = t2 ^ 2;
t35 = t34 * t4;
t36 = t35 * t23;
t37 = cos(theta_q);
t38 = t37 * t10;
t39 = coef(15);
t43 = t10 * t15;
t44 = coef(4);
t47 = sin(theta_q);
t49 = coef(16);
t53 = coef(17);
t56 = t10 * t47;
t57 = coef(1);
t60 = t34 * t3;
t61 = t60 * t7;
t62 = coef(5);
t66 = coef(6);
t70 = -t8 * t9 * t10 * t11 - t8 * t9 * t15 * t16 - t21 * t25 + t21 * t25 * t27 - t21 * t23 * t27 * t31 - t36 * t38 * t39 - t21 * t23 * t43 * t44 - t36 * t47 * t15 * t49 - t36 * t37 * t15 * t53 - t36 * t56 * t57 - t61 * t9 * t47 * t62 - t61 * t9 * t37 * t66;
t73 = t23 * t37;
t76 = t23 * t47;
t88 = t34 * t20;
t89 = t88 * t23;
t90 = coef(3);
t94 = coef(2);
t98 = coef(14);
t101 = t21 * t7;
t105 = coef(13);
t106 = t37 * t105;
t111 = t7 * j;
t112 = t35 * t111;
t113 = ko * t47;
t117 = ko * t37;
t122 = coef(11);
t129 = ko * t10;
t134 = coef(12);
t143 = t47 * t98;
t148 = -t89 * t56 * t15 * t90 - t89 * t38 * t15 * t94 - t88 * t76 * t98 + t101 * t9 * t24 * t27 + t89 * t106 * t27 - t88 * t73 * t105 - t112 * t113 * t15 * t49 - t112 * t117 * t15 * t53 - t89 * t37 * t27 * t122 - t112 * t117 * t10 * t39 - t21 * t111 * t129 * t15 * t44 - t89 * t47 * t27 * t134 - t101 * t9 * t27 * t31 - t112 * t129 * t47 * t57 + t89 * t143 * t27 - t101 * t9 * t24;
t151 = t88 * t111;
t170 = t88 * t7;
% Output
X = K(2) * t70 + K(1) * (-t60 * t73 * t66 - t60 * t76 * t62 - t5 * t23 * t10 * t11 - t5 * t23 * t15 * t16) + K(3) * t148 + K(4) * (-t151 * t113 * t43 * t90 - t151 * t117 * t27 * t122 - t151 * t117 * t43 * t94 + t151 * t117 * t105 * t27 + t151 * t113 * t98 * t27 - t151 * t113 * t27 * t134 - t170 * t9 * t106 - t170 * t9 * t143);
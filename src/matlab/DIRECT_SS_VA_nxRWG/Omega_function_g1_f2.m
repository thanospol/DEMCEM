function X = Omega_function_g1_f2(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_g1_f2

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
t6 = t5 * t4;
t7 = t3 * t6;
t8 = GAMMA ^ 2;
t10 = 0.1e1 / t8 / GAMMA;
t11 = t7 * t10;
t12 = cos(theta_q);
t13 = coef(6);
t14 = t12 * t13;
t15 = cos(theta_p);
t16 = t15 ^ 2;
t22 = sin(theta_q);
t24 = coef(7);
t27 = sin(theta_p);
t28 = t22 * t27;
t29 = coef(11);
t33 = t3 * t5;
t34 = 0.1e1 / t8;
t35 = t34 * j;
t36 = t33 * t35;
t37 = ko * t27;
t38 = coef(1);
t42 = t2 * t6;
t43 = t42 * t34;
t44 = j * ko;
t45 = coef(3);
t49 = coef(2);
t56 = coef(8);
t61 = coef(9);
t66 = coef(4);
t70 = coef(5);
t73 = t22 * t24;
t77 = coef(10);
t81 = t11 * t14 * t16 - t7 * t10 * t12 * t13 - t7 * t10 * t22 * t24 - t11 * t28 * t15 * t29 - t36 * t37 * t22 * t38 - t43 * t44 * t16 * t45 + t43 * t44 * t49 * t16 - t43 * t44 * t49 - t36 * ko * t15 * t22 * t56 - t42 * t35 * t37 * t15 * t61 - t11 * t12 * t16 * t66 - t11 * t22 * t16 * t70 + t11 * t73 * t16 - t11 * t12 * t27 * t15 * t77;
t84 = t33 * t10;
t93 = t10 * t49;
t97 = t27 * t15;
t104 = t7 * t35;
t105 = ko * t22;
t112 = ko * t12;
t122 = t7 * t34;
% output
X = K(3) * t81 + K(2) * (-t84 * t28 * t38 - t84 * t15 * t22 * t56 - t42 * t10 * t16 * t45 + t42 * t93 * t16 - t42 * t10 * t97 * t61 - t42 * t93) + K(4) * (t104 * t105 * t24 * t16 - t104 * t105 * t16 * t70 + t104 * t112 * t13 * t16 - t104 * t112 * t16 * t66 - t104 * t105 * t97 * t29 - t122 * t44 * t14 - t122 * t44 * t73 - t104 * t112 * t97 * t77);
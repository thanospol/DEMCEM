function X = Omega_function_f3_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_f3_f3

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
t11 = coef(1);
t14 = t4 ^ 2;
t15 = t2 * t14;
t16 = sin(theta_p);
t18 = coef(2);
t21 = cos(theta_q);
t23 = coef(3);
t26 = cos(theta_p);
t28 = coef(4);
t34 = t3 * t14;
t35 = t34 * t8;
t37 = coef(5);
t41 = coef(6);
t45 = coef(7);
t49 = coef(8);
t52 = 0.1e1 / t6;
t53 = t5 * t52;
t54 = j * ko;
t58 = t15 * t52;
t72 = t34 * t52 * j;
t73 = ko * t9;
t77 = ko * t21;
% Output
X = K(1) * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18 - t5 * t8 * t21 * t23 - t15 * t8 * t26 * t28) + K(2) * (-t35 * t9 * t26 * t37 - t35 * t21 * t16 * t41 - t35 * t21 * t26 * t45 - t35 * t9 * t16 * t49 - t53 * t54 * t9 * t11 - t58 * t54 * t16 * t18 - t53 * t54 * t21 * t23 - t58 * t54 * t26 * t28) + K(3) * (-t72 * t73 * t26 * t37 - t72 * t77 * t16 * t41 - t72 * t77 * t26 * t45 - t72 * t73 * t16 * t49);
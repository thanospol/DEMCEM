function X = Omega_function_f2_f2(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_f2_f2

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
t11 = cos(theta_p);
t12 = sin(theta_q);
t14 = coef(1);
t17 = cos(theta_q);
t18 = sin(theta_p);
t20 = coef(2);
t24 = coef(3);
t27 = t2 * t5;
t28 = 0.1e1 / t7;
t30 = j * ko;
t31 = coef(4);
t35 = t3 * t4;
t37 = coef(5);
t45 = t6 * t28 * j;
% Output
X = K(2) * (-t10 * t11 * t12 * t14 - t10 * t17 * t18 * t20 - t10 * t18 * t12 * t24 - t27 * t28 * t30 * t18 * t31 - t35 * t28 * t30 * t12 * t37) + K(3) * (-t45 * ko * t18 * t12 * t24 - t45 * ko * t17 * t18 * t20 - t45 * ko * t11 * t12 * t14) + K(1) * (-t27 * t9 * t18 * t31 - t35 * t9 * t12 * t37);
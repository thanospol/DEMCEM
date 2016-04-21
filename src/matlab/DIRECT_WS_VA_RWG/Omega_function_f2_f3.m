function X = Omega_function_f2_f3(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_f2_f3

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

t2 = sin(Psi);
t3 = t2 ^ 2;
t4 = cos(Psi);
t5 = t4 ^ 2;
t7 = 0.1e1 / GAMMA;
t8 = t3 * t5 * t7;
t9 = sin(theta_p);
t10 = sin(theta_q);
t15 = cos(theta_q);
t20 = cos(theta_p);
t38 = t3 * t4;
t47 = t2 * t5;
% Output
X = K(4) * (t8 * t9 * t10 * coef(9) + t8 * t15 * t9 * coef(2) + t8 * t20 * t10 * coef(3) + t8 * t15 * t20 * coef(4)) + K(2) * t2 * t4 * t7 * coef(1) + K(3) * (t38 * t7 * t15 * coef(6) + t38 * t7 * t10 * coef(7) + t47 * t7 * t9 * coef(8) + t47 * t7 * t20 * coef(5));
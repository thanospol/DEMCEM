function X = Omega_function_f3_f1(theta_p, theta_q, Psi, GAMMA, coef, K)
%% Omega_function_f3_f1

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  20 September 2011

%  Author:    Athanasios Polimeridis

% References

% A. G. Polimeridis and T. V. Yioultsis, �On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,� IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis and J. R. Mosig, �Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,� Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, �Fast and accurate
% computation of hyper-singular integrals in Galerkin surface integral equation
% formulations via the direct evaluation method,� IEEE Trans.
% Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

% A. G. Polimeridis and J. R. Mosig, �On the direct evaluation of surface integral
% equation impedance matrix elements involving point singularities,� IEEE Antennas
% Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

% INPUT DATA
% theta_p, theta_q, Psi, GAMMA, coef, K

% OUTPUT DATA
% X   
%

t2 = sin(Psi);
t3 = t2 ^ 2;
t4 = cos(Psi);
t5 = t3 * t4;
t6 = 0.1e1 / GAMMA;
t7 = sin(theta_q);
t12 = cos(theta_q);
t20 = t4 ^ 2;
t22 = t3 * t20 * t6;
t23 = cos(theta_p);
t28 = sin(theta_p);
% Output
X = K(3) * (t5 * t6 * t7 * coef(5) + t5 * t6 * t12 * coef(6)) + K(4) * (t22 * t12 * t23 * coef(3) + t22 * t28 * t7 * coef(4) + t22 * t7 * t23 * coef(1) + t22 * t12 * t28 * coef(2));
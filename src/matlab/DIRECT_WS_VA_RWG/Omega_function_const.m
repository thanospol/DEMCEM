function X = Omega_function_const( Psi, GAMMA, coef, K)
%% Omega_function_const

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  08 February 2011

%  Author:    Athanasios Polimeridis

% Reference

% A. G. Polimeridis, and T. V. Yioultsis, “On the direct evaluation of weakly singular
% integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
% Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

% A. G. Polimeridis, and J. R. Mosig, “Complete semi-analytical treatment of weakly
% singular integrals on planar triangles via the direct evaluation method,” Int. J.
% Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

% INPUT DATA
% Psi, GAMMA, coef, K

% OUTPUT DATA
% X   
%

X = K(2) * (sin(Psi) * cos(Psi) * coef(1)) / GAMMA;
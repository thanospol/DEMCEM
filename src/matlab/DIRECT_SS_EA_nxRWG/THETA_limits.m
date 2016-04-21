function [theta_A,theta_B] = THETA_limits(argument)
%% Computing limits of integration along theta 

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
% argument

% OUTPUT DATA
% theta_A, theta_B   
% 

switch argument     
    
    case 1 % I_a
        theta_A = 0;
        theta_B = pi/2;
    case 2 % I_b
        theta_A = pi/2;
        theta_B = pi;
    case 3 % I_c
        theta_A = pi/2;
        theta_B = pi;
    case 4 % I_d
        theta_A = 0;
        theta_B = pi/3;
    case 5 % I_d
        theta_A = pi/3;
        theta_B = pi/2;
    case 6 % I_d
        theta_A = pi/3;
        theta_B = pi/2;
end %switch argument
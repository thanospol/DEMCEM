function [psi_A,psi_B] = PSI_limits(theta,argument)
%% Computing limits of integration along Psi 

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
% psi_A, psi_B   
% 

PsiA=atan(sin(theta)-sqrt(3)*cos(theta));
PsiB=atan(sin(theta)+sqrt(3)*cos(theta));
%
switch argument     
    
    case 1 % I_a
        psi_A = PsiB;
        psi_B = pi/2;
    case 2 % I_b
        psi_A = PsiA;
        psi_B = pi/2;
    case 3 % I_c
        psi_A = 0;
        psi_B = PsiA;
    case 4 % I_d
        psi_A = 0;
        psi_B = PsiB;
    case 5 % I_d
        psi_A = PsiA;
        psi_B = PsiB;
    case 6 % I_d
        psi_A = 0;
        psi_B = PsiA;
end %switch argument
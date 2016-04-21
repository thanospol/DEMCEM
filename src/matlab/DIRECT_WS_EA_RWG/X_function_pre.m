function [N,Nm] = X_function_pre( theta, Psi, tPsiA, tPsiB, PsiA, PsiB, B, Bm, ko)
%% X_function_pre

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
% theta, Psi, tPsiA, tPsiB, PsiA, PsiB, B, Bm, ko

% OUTPUT DATA
% N   
% Nm
%

 if Psi >= PsiB && Psi >= PsiA,
     D = sqrt(3)/sin(Psi);
     N = X2(D,B,ko);
     Nm = X2(D,Bm,ko);
 elseif theta >= pi/2,
     D = sqrt(3)/(cos(Psi)*tPsiA);
     N = X2(D,B,ko);
     Nm = X2(D,Bm,ko);
 elseif Psi >= PsiA,
     D1 = 2*sqrt(3)/(cos(Psi)*(tPsiB+tan(Psi))); 
     D2 = sin(Psi)/sqrt(3); 
     N = X1(Psi,D1,D2,tPsiB,B,ko); 
     Nm = X1(Psi,D1,D2,tPsiB,Bm,ko); 
 else
     D1 = sqrt(3)/(cos(Psi)*sin(theta));
     D2 = (cos(Psi)*tPsiA)/sqrt(3);
     N = X1(Psi,D1,D2,tPsiB,B,ko);
     Nm = X1(Psi,D1,D2,tPsiB,Bm,ko);
 end
function N = X2(D,B,ko)
%% Computing functions N^c and N^d 

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
% psi,D1,D2,tpsiB,B,ko

% OUTPUT DATA
% N(12)   
% 

 j = sqrt(-1);
 a = j*ko*B;
 aD = a*D;
 expaD = exp(-aD);
 T1 = (1-expaD)/aD;
 T2 = (1-T1)/aD;
 T3 = (1-2*T2)/aD;
 T4 = (1-3*T3)/aD;

 N = [1,T1,D*(T1-T2),0,0,1/2,1/3,T2,T3,D*(T2-T3),D*(T3-T4),0];
 N(12) = D*(N(10)-N(11));
 N(4) = D*(N(3)-N(10));
 N(5) = D*(N(4)-N(12));
function N = X1(psi,D1,D2,tpsiB,B,ko)
%% Computing functions N^a and N^b 

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
 D3 = sqrt(3)/(cos(psi)*tpsiB);
 aD3 = a*D3;
 expaD3 = exp(-aD3);
 aD1 = a*D1;
 expaD1 = exp(-aD1);
 D2a = D2/a;

 H=1-D1*D2;

 T11 = (expaD3-expaD1)/aD3;
 T21 = (T11-H*expaD1)/aD3;
 T31 = (2*T21-H^2*expaD1)/aD3;
 T41 = (3*T31-H^3*expaD1)/aD3;

 T12 = D2a*(1-expaD1);
 T22 = D2a*(1-H*expaD1-T12);
 T32 = D2a*(1-H^2*expaD1-2*T22);
 T42 = D2a*(1-H^3*expaD1-3*T32);

 N21 = T11;
 N31 = D3*(T11+T21);
 N81 = T21;
 N91 = T31;
 N101 = D3*(T21+T31);
 N111 = D3*(T31+T41);
 N121 = D3*(N101+N111);
 N41 = D3*(N31+N101);
 N51 = D3*(N41+N121);

 N22 = T12;
 N32 = (T12-T22)/D2;
 N82 = T22;
 N92 = T32;
 N102 = (T22-T32)/D2;
 N112 = (T32-T42)/D2;
 N122 = (N102-N112)/D2;
 N42 = (N32-N102)/D2;
 N52 = (N42-N122)/D2;

 N = [1,N21+N22,N31+N32,N41+N42,N51+N52,1/2,1/3,N81+N82,N91+N92,N101+N102,N111+N112,N121+N122];
 
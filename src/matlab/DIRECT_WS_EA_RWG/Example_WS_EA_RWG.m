%% Example for WS_EA_RWG

%  Licensing: This code is distributed under the GNU LGPL license. 

%  Modified:  19 October 2011

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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)

% global ko = wavenumber
% N_theta,N_psi = order of the Gauss-Legendre quadrature for both dimensions
% of the remaining 2-D smooth integral

% OUTPUT DATA
%             I_DE(1)  = I_f1_f1
%             I_DE(2)  = I_f1_f2
%             I_DE(3)  = I_f1_f3
%             I_DE(4)  = I_f2_f1
%             I_DE(5)  = I_f2_f2
%             I_DE(6)  = I_f2_f3
%             I_DE(7)  = I_f3_f1
%             I_DE(8)  = I_f3_f2
%             I_DE(9)  = I_f3_f3
%

clc
clear all;
format long;

global ko
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1 = [0.0;0.0;0.0];
r2 = [0.5;0.0;0.0];
r3 = [0.0;0.0;0.5];
r4 = [0.0;0.5;0.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ko    = 1;
%
Np_1D = 10;
%
N_theta = Np_1D;
N_psi = Np_1D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic  
for ii=1:100
[I_DE] = DIRECT_WS_EA_RWG(r1,r2,r3,r4,N_theta,N_psi);
end
Time = toc/100*1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PRINT RESULTS                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Np_1D: %d \n', Np_1D)
fprintf('RunTime: %f [msec] \n\n', Time)
fprintf('I_f1_f1: (%4.16e, %4.16e) \n', real(I_DE(1)),imag(I_DE(1)))
fprintf('I_f1_f2: (%4.16e, %4.16e) \n', real(I_DE(2)),imag(I_DE(2)))
fprintf('I_f1_f3: (%4.16e, %4.16e) \n', real(I_DE(3)),imag(I_DE(3)))
fprintf('I_f2_f1: (%4.16e, %4.16e) \n', real(I_DE(4)),imag(I_DE(4)))
fprintf('I_f2_f2: (%4.16e, %4.16e) \n', real(I_DE(5)),imag(I_DE(5)))
fprintf('I_f2_f3: (%4.16e, %4.16e) \n', real(I_DE(6)),imag(I_DE(6)))
fprintf('I_f3_f1: (%4.16e, %4.16e) \n', real(I_DE(7)),imag(I_DE(7)))
fprintf('I_f3_f2: (%4.16e, %4.16e) \n', real(I_DE(8)),imag(I_DE(8)))
fprintf('I_f3_f3: (%4.16e, %4.16e) \n\n', real(I_DE(9)),imag(I_DE(9)))




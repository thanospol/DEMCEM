function [PHI] = Phi_1_1(argument,psi,Apsi)
%% Phi_1_1 function

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
% argument  = 1-8 -> a-h
% psi       = Variable Psi
% Apsi      = alpha(Psi)

% OUTPUT DATA
% Phi_1_1   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_1_1                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
D = sin(psi)/((j*ko)^2*Apsi^3);

switch argument

case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_a1 = A/(2*B);
PSI_a2 = -1;
PSI_a3 = (B/A)*(1-exp(-A/B));
PSI_a  = PSI_a1+PSI_a2+PSI_a3;
PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_b1 = A/(2*C);
PSI_b2 = -1;
PSI_b3 = (C/A)*(1-exp(-A/C));
PSI_b  = PSI_b1+PSI_b2+PSI_b3;
PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_c1 = (A/C)*(1/2-(gamma-gamma^2/2));
PSI_c2 = -(1-gamma);
PSI_c3 = (C/A)*(1-exp(-A*(1-gamma)/C));
PSI_c  = PSI_c1+PSI_c2+PSI_c3;
PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_d1 = -A*(gamma+gamma^2/2)/B;
PSI_d2 = -gamma;
PSI_d3 = (B/A)*(exp(A*(1+gamma)/B)-exp(A/B));
PSI_d  = PSI_d1+PSI_d2+PSI_d3;
PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_e1 = A*(delta^2/2-delta)/B;
PSI_e2 = delta;
PSI_e3 = (B/A)*(exp(-A/B)-exp(-A*(1-delta)/B));
PSI_e  = PSI_e1+PSI_e2+PSI_e3;
PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_f1 = (A/C)*(delta+delta^2/2+1/2);
PSI_f2 = -(1+delta);
PSI_f3 = (C/A)*(1-exp(-A*(1+delta)/C));
PSI_f  = PSI_f1+PSI_f2+PSI_f3;
PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_g1 = A/(2*C);
PSI_g2 = -1;
PSI_g3 = (C/A)*(1-exp(-A/C));
PSI_g  = PSI_g1+PSI_g2+PSI_g3;
PHI    = D*PSI_g;

case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_h1 = -A/(2*B);
PSI_h2 = -1;
PSI_h3 = (B/A)*(exp(A/B)-1);
PSI_h  = PSI_h1+PSI_h2+PSI_h3;
PHI    = D*PSI_h;

end %switch argument
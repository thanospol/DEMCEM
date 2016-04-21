function [PHI] = Phi_3_1(argument,psi,Apsi)
%% Phi_3_1 function

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
% Phi_3_1   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_3_1                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sin(psi)^2/(j*ko*Apsi^2);

switch argument
case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_2_a = 1/(3*B^2);

PSI_5_a1 = 1;
PSI_5_a2 = -B/A*(1-exp(-A/B));
PSI_5_a3 = -1/A*(B-(A+B)*exp(-A/B));
PSI_5_a  = PSI_5_a1+PSI_5_a2+PSI_5_a3;

PSI_a = PSI_2_a/2-PSI_5_a/(j*ko*Apsi)^2;

PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_2_b = 1/(3*C^2);

PSI_5_b1 = 1;
PSI_5_b2 = -C/A*(1-exp(-A/C));
PSI_5_b3 = -1/A*(C-(A+C)*exp(-A/C));
PSI_5_b  = PSI_5_b1+PSI_5_b2+PSI_5_b3;

PSI_b = PSI_2_b/2-PSI_5_b/(j*ko*Apsi)^2;

PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_2_c = (1-gamma)^3/(3*C^2);

PSI_5_c1 = 1-gamma;
PSI_5_c2 = -C/A*(1-exp(-A*(1-gamma)/C));
PSI_5_c3 = -1/A*(C+(A*gamma-C-A)*exp(-A*(1-gamma)/C));
PSI_5_c  = PSI_5_c1+PSI_5_c2+PSI_5_c3;

PSI_c = PSI_2_c/2-PSI_5_c/(j*ko*Apsi)^2;

PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_2_d = ((1+gamma)^3-1)/(3*B^2);

PSI_5_d1 = gamma;
PSI_5_d2 = -B/A*(exp(A*(1+gamma)/B)-exp(A/B));
PSI_5_d3 = -1/A*((A-B)*exp(A/B)+(B-A-A*gamma)*exp(A*(1+gamma)/B));
PSI_5_d  = PSI_5_d1+PSI_5_d2+PSI_5_d3;

PSI_d = PSI_2_d/2-PSI_5_d/(j*ko*Apsi)^2;

PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_2_e = ((1-delta)^3-1)/(3*B^2);

PSI_5_e1 = -delta;
PSI_5_e2 = -B/A*(exp(-A/B)-exp(-A*(1-delta)/B));
PSI_5_e3 = -1/A*((A+B)*exp(-A/B)+(A*delta-A-B)*exp(-A*(1-delta)/B));
PSI_5_e  = PSI_5_e1+PSI_5_e2+PSI_5_e3;

PSI_e = PSI_2_e/2-PSI_5_e/(j*ko*Apsi)^2;

PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_2_f = (1+delta)^3/(3*C^2);

PSI_5_f1 = 1+delta;
PSI_5_f2 = -C/A*(1-exp(-A*(1+delta)/C));
PSI_5_f3 = 1/A*((A+A*delta+C)*exp(-A*(1+delta)/C)-C);
PSI_5_f  = PSI_5_f1+PSI_5_f2+PSI_5_f3;

PSI_f = PSI_2_f/2-PSI_5_f/(j*ko*Apsi)^2;

PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_2_g = 1/(3*C^2);

PSI_5_g1 = 1;
PSI_5_g2 = -C/A*(1-exp(-A/C));
PSI_5_g3 = 1/A*((A+C)*exp(-A/C)-C);
PSI_5_g  = PSI_5_g1+PSI_5_g2+PSI_5_g3;

PSI_g = PSI_2_g/2-PSI_5_g/(j*ko*Apsi)^2;

PHI    = D*PSI_g;
case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_2_h = 1/(3*B^2);

PSI_5_h1 = 1;
PSI_5_h2 = -B/A*(exp(A/B)-1);
PSI_5_h3 = 1/A*((A-B)*exp(A/B)+B);
PSI_5_h  = PSI_5_h1+PSI_5_h2+PSI_5_h3;

PSI_h = PSI_2_h/2-PSI_5_h/(j*ko*Apsi)^2;

PHI    = D*PSI_h;

end %switch argument
function [PHI] = Phi_1_3(argument,psi,Apsi)
%% Phi_1_3 function

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
% Phi_1_3   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_1_3                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sin(psi)/((j*ko)^2*Apsi^2);

switch argument
    
case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_a = 1/(3*B^2);

PSI_1_a = 1/(2*B);

PSI_2_a  = 1-B/A*(1-exp(-A/B));

PSI_a = (j*ko*sin(psi)/2)*PSI_0_a-(sin(psi)/Apsi)*PSI_1_a+(sin(psi)/(j*ko*Apsi^2))*PSI_2_a;

PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_b = 1/(3*C^2);

PSI_1_b = 1/(2*C);

PSI_2_b  = 1-C/A*(1-exp(-A/C));

PSI_b = (j*ko*sin(psi)/2)*PSI_0_b-(sin(psi)/Apsi)*PSI_1_b+(sin(psi)/(j*ko*Apsi^2))*PSI_2_b;

PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_c = (1-gamma)^3/(3*C^2);

PSI_1_c = (1/C)*(1/2-(gamma-gamma^2/2));

PSI_2_c  = (1-gamma)-C/A*(1-exp(-A*(1-gamma)/C));


PSI_c = (j*ko*sin(psi)/2)*PSI_0_c-(sin(psi)/Apsi)*PSI_1_c+(sin(psi)/(j*ko*Apsi^2))*PSI_2_c;

PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_d = ((1+gamma)^3-1)/(3*B^2);

PSI_1_d = -1/B*(gamma+gamma^2/2);

PSI_2_d  = gamma-B/A*(exp(A*(1+gamma)/B)-exp(A/B));

PSI_d = (j*ko*sin(psi)/2)*PSI_0_d-(sin(psi)/Apsi)*PSI_1_d+(sin(psi)/(j*ko*Apsi^2))*PSI_2_d;

PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_e = ((1-delta)^3-1)/(3*B^2);

PSI_1_e = 1/B*(delta^2/2-delta);

PSI_2_e  = -delta-B/A*(exp(-A/B)-exp(-A*(1-delta)/B));

PSI_e = (j*ko*sin(psi)/2)*PSI_0_e-(sin(psi)/Apsi)*PSI_1_e+(sin(psi)/(j*ko*Apsi^2))*PSI_2_e;

PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_f = (1+delta)^3/(3*C^2);

PSI_1_f = 1/C*(delta+delta^2/2+1/2);

PSI_2_f  = (1+delta)-C/A*(1-exp(-A*(1+delta)/C));

PSI_f = (j*ko*sin(psi)/2)*PSI_0_f-(sin(psi)/Apsi)*PSI_1_f+(sin(psi)/(j*ko*Apsi^2))*PSI_2_f;

PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_g = 1/(3*C^2);

PSI_1_g = 1/(2*C);

PSI_2_g  = 1-C/A*(1-exp(-A/C));

PSI_g = (j*ko*sin(psi)/2)*PSI_0_g-(sin(psi)/Apsi)*PSI_1_g+(sin(psi)/(j*ko*Apsi^2))*PSI_2_g;

PHI    = D*PSI_g;
case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_h = 1/(3*B^2);

PSI_1_h = -1/(2*B);

PSI_2_h  = 1-B/A*(exp(A/B)-1);

PSI_h = (j*ko*sin(psi)/2)*PSI_0_h-(sin(psi)/Apsi)*PSI_1_h+(sin(psi)/(j*ko*Apsi^2))*PSI_2_h;

PHI    = D*PSI_h;

end %switch argument
function [PHI] = Phi_2_2(argument,psi,Apsi)
%% Phi_2_2 function

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
% Phi_2_2   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_2_2                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sin(psi)/((j*ko)^2*Apsi^2);

switch argument

    case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_a = 1/(6*B);

PSI_1_a = 1/(12*B);

PSI_2_a  = 1/2-B/A+(B^2/A^2)*(1-exp(-A/B));

PSI_5_a  = B/A-2*B^2/A^2+2*B^3/A^3*(1-exp(-A/B));

PSI_3_a1 = B/A-(B^2/A^2)*(1-exp(-A/B));
PSI_3_a2 = PSI_5_a;
PSI_3_a  = 1/B*(PSI_3_a1-PSI_3_a2);

PSI_4_a  = 1/3;

PSI_a = (cos(psi)/Apsi)*PSI_0_a+j*ko*PSI_1_a-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_a+(cos(psi)/Apsi)*PSI_3_a-(1/Apsi)*(PSI_4_a-PSI_5_a);

PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_b = 1/(6*C);

PSI_1_b = 1/(12*C);

PSI_2_b =  1/2-C/A+(C^2/A^2)*(1-exp(-A/C));

PSI_5_b  = C/A-2*C^2/A^2+2*C^3/A^3*(1-exp(-A/C));

PSI_3_b1 = C/A-(C^2/A^2)*(1-exp(-A/C));
PSI_3_b2 = PSI_5_b;
PSI_3_b  = 1/C*(PSI_3_b1-PSI_3_b2);

PSI_4_b  = 1/3;

PSI_b = (cos(psi)/Apsi)*PSI_0_b+j*ko*PSI_1_b-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_b+(cos(psi)/Apsi)*PSI_3_b-(1/Apsi)*(PSI_4_b-PSI_5_b);

PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_c = (1/C)*((1/2)*(1-gamma^2)-(1/3)*(1-gamma^3));

PSI_1_c = (1/C)*(1/12-(gamma^3/3-gamma^4/4));

PSI_2_c = 1/2*(1-gamma^2)+C/A^2*(C-A+(A*gamma-C)*exp(-A*(1-gamma)/C));

PSI_5_c = C/A*(1-gamma^2*exp(-A*(1-gamma)/C))-2*C^2/A^2*(1-gamma*exp(-A*(1-gamma)/C))+2*C^3/A^3*(1-exp(-A*(1-gamma)/C));

PSI_3_c1 = -C/A^2*(C-A+(A*gamma-C)*exp(-A*(1-gamma)/C));
PSI_3_c2 = PSI_5_c;
PSI_3_c  = 1/C*(PSI_3_c1-PSI_3_c2);

PSI_4_c  = (1-gamma^3)/3;

PSI_c = (cos(psi)/Apsi)*PSI_0_c+j*ko*PSI_1_c-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_c+(cos(psi)/Apsi)*PSI_3_c-(1/Apsi)*(PSI_4_c-PSI_5_c);

PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_d = -1/B*(gamma^2/2+gamma^3/3);

PSI_1_d = -1/B*(gamma^3/3+gamma^4/4);

PSI_2_d = gamma^2/2-B/A^2*(B*exp(A/B)+(A*gamma-B)*exp(A*(1+gamma)/B));

PSI_5_d = B/A*gamma^2*exp(A*(1+gamma)/B)-2*B^2/A^2*gamma*exp(A*(1+gamma)/B)+2*B^3/A^3*(exp(A*(1+gamma)/B)-exp(A/B));

PSI_3_d1 = B/A^2*(B*exp(A/B)+(A*gamma-B)*exp(A*(1+gamma)/B));
PSI_3_d2 = PSI_5_d;
PSI_3_d  = -1/B*(PSI_3_d1+PSI_3_d2);

PSI_4_d  = gamma^3/3;

PSI_d = (cos(psi)/Apsi)*PSI_0_d+j*ko*PSI_1_d-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_d+(cos(psi)/Apsi)*PSI_3_d-(1/Apsi)*(PSI_4_d-PSI_5_d);

PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_e = 1/B*(delta^3/3-delta^2/2);

PSI_1_e = 1/B*(delta^4/4-delta^3/3);

PSI_2_e = -delta^2/2+B/A^2*(B*exp(-A/B)+(A*delta-B)*exp(-A*(1-delta)/B));

PSI_5_e = -B/A*delta^2*exp(-A*(1-delta)/B)+2*B^2/A^2*delta*exp(-A*(1-delta)/B)+2*B^3/A^3*(exp(-A/B)-exp(-A*(1-delta)/B));

PSI_3_e1 = -B/A^2*(B*exp(-A/B)+(A*delta-B)*exp(-A*(1-delta)/B));
PSI_3_e2 = PSI_5_e;
PSI_3_e  = 1/B*(PSI_3_e1-PSI_3_e2);

PSI_4_e  = -delta^3/3;

PSI_e = (cos(psi)/Apsi)*PSI_0_e+j*ko*PSI_1_e-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_e+(cos(psi)/Apsi)*PSI_3_e-(1/Apsi)*(PSI_4_e-PSI_5_e);

PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_f = 1/C*(1/2*(delta^2-1)+1/3*(delta^3+1));

PSI_1_f = 1/C*(delta^3/3+delta^4/4+1/12);

PSI_2_f = 1/2*(delta^2-1)+C/A^2*(A-C+(A*delta+C)*exp(-A*(1+delta)/C));

PSI_5_f = C/A*(1-delta^2*exp(-A*(1+delta)/C))-2*C^2/A^2*(delta*exp(-A*(1+delta)/C)+1)-2*C^3/A^3*(exp(-A*(1+delta)/C)-1);

PSI_3_f1 = -C/A^2*(A-C+(A*delta+C)*exp(-A*(1+delta)/C));
PSI_3_f2 = PSI_5_f;
PSI_3_f  = 1/C*(PSI_3_f1+PSI_3_f2);

PSI_4_f  = (delta^3+1)/3;

PSI_f = (cos(psi)/Apsi)*PSI_0_f+j*ko*PSI_1_f-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_f+(cos(psi)/Apsi)*PSI_3_f-(1/Apsi)*(PSI_4_f-PSI_5_f);

PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_g = -1/(6*C);

PSI_1_g = 1/(12*C);

PSI_2_g = -1/2+C/A^2*(C*exp(-A/C)-C+A);

PSI_5_g = C/A-2*C^2/A^2-2*C^3/A^3*(exp(-A/C)-1);

PSI_3_g1 = -C/A^2*(C*exp(-A/C)-C+A);
PSI_3_g2 = PSI_5_g;
PSI_3_g  = 1/C*(PSI_3_g1+PSI_3_g2);

PSI_4_g  = 1/3;

PSI_g = (cos(psi)/Apsi)*PSI_0_g+j*ko*PSI_1_g-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_g+(cos(psi)/Apsi)*PSI_3_g-(1/Apsi)*(PSI_4_g-PSI_5_g);

PHI    = D*PSI_g;
case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_h = 1/(6*B);

PSI_1_h = -1/(12*B);

PSI_2_h = -1/2-B/A^2*(A+B-B*exp(A/B));

PSI_5_h = -B/A-2*B^2/A^2+2*B^3/A^3*(exp(A/B)-1);

PSI_3_h1 = B/A^2*(A+B-B*exp(A/B));
PSI_3_h2 = PSI_5_h;
PSI_3_h  = -1/B*(PSI_3_h1+PSI_3_h2);

PSI_4_h  = 1/3;

PSI_h = (cos(psi)/Apsi)*PSI_0_h+j*ko*PSI_1_h-(2*cos(psi)/(j*ko*Apsi^2))*PSI_2_h+(cos(psi)/Apsi)*PSI_3_h-(1/Apsi)*(PSI_4_h-PSI_5_h);

PHI    = D*PSI_h;

end %switch argument
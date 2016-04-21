function [PHI] = Phi_1_2(argument,psi,Apsi)
%% Phi_1_2 function

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
% Phi_1_2   
%
global ko
j = sqrt(-1);

A = j*ko*Apsi;
B = cos(psi);
C = sin(psi)/sqrt(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Int_1_2                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
D = sin(psi)/((j*ko)^2*Apsi^2);
        
switch argument

case 1
%%%%%%%%%%%%%%%%%%%%%%    PHI_a   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_a = 1/(6*B);

PSI_1_a = 1/(2*B);

PSI_2_a1 = 1;
PSI_2_a2 = -B/A*(1-exp(-A/B));
PSI_2_a3 = -1/A*(B-(A+B)*exp(-A/B));
PSI_2_a  = PSI_2_a1+PSI_2_a2+PSI_2_a3;

PSI_3_a  = 1-B/A*(1-exp(-A/B));

PSI_4_a  = 1/2-B/A+(B^2/A^2)*(1-exp(-A/B));

PSI_a = j*ko*PSI_0_a+(cos(psi)/Apsi)*PSI_1_a-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_a+PSI_3_a)+(1/Apsi)*PSI_4_a);

PHI    = D*PSI_a;

case 2
%%%%%%%%%%%%%%%%%%%%%%    PHI_b   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_b = 1/(6*C);

PSI_1_b = 1/(2*C);

PSI_2_b1 = 1;
PSI_2_b2 = -C/A*(1-exp(-A/C));
PSI_2_b3 = -1/A*(C-(A+C)*exp(-A/C));
PSI_2_b  = PSI_2_b1+PSI_2_b2+PSI_2_b3;

PSI_3_b  = 1-C/A*(1-exp(-A/C));

PSI_4_b  = 1/2-C/A+(C^2/A^2)*(1-exp(-A/C));

PSI_b = j*ko*PSI_0_b+(cos(psi)/Apsi)*PSI_1_b-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_b+PSI_3_b)+(1/Apsi)*PSI_4_b);

PHI    = D*PSI_b;

case 3
%%%%%%%%%%%%%%%%%%%%%%    PHI_c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_c = (1/C)*((1/2)*(1-gamma^2)-(1/3)*(1-gamma^3));

PSI_1_c = (1/C)*(1/2-(gamma-gamma^2/2));

PSI_2_c1 = 1-gamma;
PSI_2_c2 = -C/A*(1-exp(-A*(1-gamma)/C));
PSI_2_c3 = -1/A*(C+(A*gamma-C-A)*exp(-A*(1-gamma)/C));
PSI_2_c  = PSI_2_c1+PSI_2_c2+PSI_2_c3;

PSI_3_c  = (1-gamma)-C/A*(1-exp(-A*(1-gamma)/C));

PSI_4_c  = 1/2*(1-gamma^2)+C/A^2*(C-A+(A*gamma-C)*exp(-A*(1-gamma)/C));

PSI_c = j*ko*PSI_0_c+(cos(psi)/Apsi)*PSI_1_c-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_c+PSI_3_c)+(1/Apsi)*PSI_4_c);

PHI    = D*PSI_c;

case 4
%%%%%%%%%%%%%%%%%%%%%%    PHI_d   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta    = tan(pi-psi)/sqrt(3);
gamma   = (1-beta)/(1+beta);

PSI_0_d = -1/B*(gamma^2/2+gamma^3/3);

PSI_1_d = -1/B*(gamma+gamma^2/2);

PSI_2_d1 = gamma;
PSI_2_d2 = -B/A*(exp(A*(1+gamma)/B)-exp(A/B));
PSI_2_d3 = -1/A*((A-B)*exp(A/B)+(B-A-A*gamma)*exp(A*(1+gamma)/B));
PSI_2_d  = PSI_2_d1+PSI_2_d2+PSI_2_d3;

PSI_3_d  = gamma-B/A*(exp(A*(1+gamma)/B)-exp(A/B));

PSI_4_d  = gamma^2/2-B/A^2*(B*exp(A/B)+(A*gamma-B)*exp(A*(1+gamma)/B));

PSI_d = j*ko*PSI_0_d+(cos(psi)/Apsi)*PSI_1_d-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_d+PSI_3_d)+(1/Apsi)*PSI_4_d);

PHI    = D*PSI_d;

case 5
%%%%%%%%%%%%%%%%%%%%%%    PHI_e   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_e = 1/B*(delta^3/3-delta^2/2);

PSI_1_e = 1/B*(delta^2/2-delta);

PSI_2_e1 = -delta;
PSI_2_e2 = -B/A*(exp(-A/B)-exp(-A*(1-delta)/B));
PSI_2_e3 = -1/A*((A+B)*exp(-A/B)+(A*delta-A-B)*exp(-A*(1-delta)/B));
PSI_2_e  = PSI_2_e1+PSI_2_e2+PSI_2_e3;

PSI_3_e  = -delta-B/A*(exp(-A/B)-exp(-A*(1-delta)/B));

PSI_4_e  = -delta^2/2+B/A^2*(B*exp(-A/B)+(A*delta-B)*exp(-A*(1-delta)/B));

PSI_e = j*ko*PSI_0_e+(cos(psi)/Apsi)*PSI_1_e-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_e+PSI_3_e)+(1/Apsi)*PSI_4_e);

PHI    = D*PSI_e;

case 6
%%%%%%%%%%%%%%%%%%%%%%    PHI_f   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = tan(psi)/sqrt(3);
delta   = -(1-epsilon)/(1+epsilon);

PSI_0_f = 1/C*(1/2*(delta^2-1)+1/3*(delta^3+1));

PSI_1_f = 1/C*(delta+delta^2/2+1/2);

PSI_2_f1 = 1+delta;
PSI_2_f2 = -C/A*(1-exp(-A*(1+delta)/C));
PSI_2_f3 = 1/A*((A+A*delta+C)*exp(-A*(1+delta)/C)-C);
PSI_2_f  = PSI_2_f1+PSI_2_f2+PSI_2_f3;

PSI_3_f  = (1+delta)-C/A*(1-exp(-A*(1+delta)/C));

PSI_4_f  = 1/2*(delta^2-1)+C/A^2*(A-C+(A*delta+C)*exp(-A*(1+delta)/C));

PSI_f = j*ko*PSI_0_f+(cos(psi)/Apsi)*PSI_1_f-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_f+PSI_3_f)+(1/Apsi)*PSI_4_f);

PHI    = D*PSI_f;

case 7
%%%%%%%%%%%%%%%%%%%%%%    PHI_g   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_g = -1/(6*C);

PSI_1_g = 1/(2*C);

PSI_2_g1 = 1;
PSI_2_g2 = -C/A*(1-exp(-A/C));
PSI_2_g3 = 1/A*((A+C)*exp(-A/C)-C);
PSI_2_g  = PSI_2_g1+PSI_2_g2+PSI_2_g3;

PSI_3_g  = 1-C/A*(1-exp(-A/C));

PSI_4_g  = -1/2+C/A^2*(C*exp(-A/C)-C+A);

PSI_g = j*ko*PSI_0_g+(cos(psi)/Apsi)*PSI_1_g-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_g+PSI_3_g)+(1/Apsi)*PSI_4_g);

PHI    = D*PSI_g;
case 8
%%%%%%%%%%%%%%%%%%%%%%    PHI_h   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI_0_h = 1/(6*B);

PSI_1_h = -1/(2*B);

PSI_2_h1 = 1;
PSI_2_h2 = -B/A*(exp(A/B)-1);
PSI_2_h3 = 1/A*((A-B)*exp(A/B)+B);
PSI_2_h  = PSI_2_h1+PSI_2_h2+PSI_2_h3;

PSI_3_h  = 1-B/A*(exp(A/B)-1);

PSI_4_h  = -1/2-B/A^2*(A+B-B*exp(A/B));

PSI_h = j*ko*PSI_0_h+(cos(psi)/Apsi)*PSI_1_h-((cos(psi)/(j*ko*Apsi^2))*(PSI_2_h+PSI_3_h)+(1/Apsi)*PSI_4_h);

PHI    = D*PSI_h;

end % switch argument
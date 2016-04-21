function [I_DE] = DIRECT_SS_EA_nxRWG(r1,r2,r3,r4,N_theta,N_psi)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the edge adjacent 4-D strongly singular integrals over planar
% triangular elements.

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
% r1,r2,r3, r4 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
% N_theta,N_psi = order of the Gauss-Legendre quadrature for both dimensions
% of the remaining 2-D smooth integral

% OUTPUT DATA
% I_DE(10)  = scalar and potential 4-D strongly singular integrals
%             I_DE(1)  = I_g1_f1
%             I_DE(2)  = I_g1_f2
%             I_DE(3)  = I_g1_f3
%             I_DE(4)  = I_g2_f1
%             I_DE(5)  = I_g2_f2
%             I_DE(6)  = I_g2_f3
%             I_DE(7)  = I_g3_f1
%             I_DE(8)  = I_g3_f2
%             I_DE(9)  = I_g3_f3
%

global ko
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Triangle P                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rp1 = r1;
rp2 = r2;
rp3 = r3;
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ap = (1/2)*norm(cross(rp2-rp1,rp3-rp1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jp = Ap/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Triangle Q                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rq1 = r2;
rq2 = r1;
rq3 = r4;
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aq = (1/2)*norm(cross(rq2-rq1,rq3-rq1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jq = Aq/sqrt(3);
%
alpha = (r2-r1)/2;
beta  = (2*r3-r1-r2)/(2*sqrt(3));
gamma = -(2*r4-r1-r2)/(2*sqrt(3));
% 
[coef_g1_f1,coefm_g1_f1] = coefficients_g1_f1(r1,r2,r3,r4,ko,Ap);
[coef_g1_f2,coefm_g1_f2] = coefficients_g1_f2(r1,r2,r3,r4,ko,Ap);
[coef_g1_f3,coefm_g1_f3] = coefficients_g1_f3(r1,r2,r3,r4,ko,Ap);
[coef_g2_f1,coefm_g2_f1] = coefficients_g2_f1(r1,r2,r3,r4,ko,Ap);
[coef_g2_f2,coefm_g2_f2] = coefficients_g2_f2(r1,r2,r3,r4,ko,Ap);
[coef_g2_f3,coefm_g2_f3] = coefficients_g2_f3(r1,r2,r3,r4,ko,Ap);
[coef_g3_f1,coefm_g3_f1] = coefficients_g3_f1(r1,r2,r3,r4,ko,Ap);
[coef_g3_f2,coefm_g3_f2] = coefficients_g3_f2(r1,r2,r3,r4,ko,Ap);
[coef_g3_f3,coefm_g3_f3] = coefficients_g3_f3(r1,r2,r3,r4,ko,Ap);
%
[ w_theta, z_theta ] = Gauss_1D ( N_theta);
[ w_psi, z_psi ] = Gauss_1D ( N_psi );
%
I_g1_f1 = zeros(6,1);
I_g1_f2 = zeros(6,1);
I_g1_f3 = zeros(6,1);
I_g2_f1 = zeros(6,1);
I_g2_f2 = zeros(6,1);
I_g2_f3 = zeros(6,1);
I_g3_f1 = zeros(6,1);
I_g3_f2 = zeros(6,1);
I_g3_f3 = zeros(6,1);
%  
 for m = 1:6
     I_theta_g1_f1 = 0;
     I_theta_g1_f2 = 0;
     I_theta_g1_f3 = 0;
     I_theta_g2_f1 = 0;
     I_theta_g2_f2 = 0;
     I_theta_g2_f3 = 0;
     I_theta_g3_f1 = 0;
     I_theta_g3_f2 = 0;
     I_theta_g3_f3 = 0;
     for n_theta = 1:N_theta
         [theta_A,theta_B] = THETA_limits(m);
         THETA = ((theta_B-theta_A)/2)*z_theta(n_theta)+(theta_B+theta_A)/2;
         J_theta = (theta_B-theta_A)/2;
         %
         tPsiA=sin(THETA)-sqrt(3)*cos(THETA);
         tPsiB=sin(THETA)+sqrt(3)*cos(THETA);
         %
         PsiA=atan(tPsiA);
         PsiB=atan(tPsiB);
        %
         b0 = beta.'*beta;
         b1 = 2*((alpha.'*beta)*cos(THETA)+(beta.'*gamma)*sin(THETA));
         b2 = norm(alpha)^2*cos(THETA)^2+2*(alpha.'*gamma)*cos(THETA)*sin(THETA)+norm(gamma)^2*sin(THETA)^2;
        %
         B1 = 2*((beta.'*gamma)*sin(THETA)-(alpha.'*beta)*cos(THETA));
         B2 = norm(alpha)^2*cos(THETA)^2-2*(alpha.'*gamma)*cos(THETA)*sin(THETA)+norm(gamma)^2*sin(THETA)^2;
         %
         I_psi_g1_f1 = 0;
         I_psi_g1_f2 = 0;
         I_psi_g1_f3 = 0;
         I_psi_g2_f1 = 0;
         I_psi_g2_f2 = 0;
         I_psi_g2_f3 = 0;
         I_psi_g3_f1 = 0;
         I_psi_g3_f2 = 0;
         I_psi_g3_f3 = 0;
         for n_psi = 1:N_psi
             [psi_A,psi_B] = PSI_limits(THETA,m);
             PSI = ((psi_B-psi_A)/2)*z_psi(n_psi)+(psi_B+psi_A)/2;
             J_psi = (psi_B-psi_A)/2;
             %
             B = sqrt(b0*sin(PSI)^2+b1*cos(PSI)*sin(PSI)+b2*cos(PSI)^2);
             Bm = sqrt(b0*sin(PSI)^2+B1*cos(PSI)*sin(PSI)+B2*cos(PSI)^2);
             %
             [N,Nm] = X_function_pre ( THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, ko);
             %
             X_g1_f1 = X_function_g1_f1(THETA, PSI, B, Bm, coef_g1_f1, coefm_g1_f1, N, Nm);
             X_g1_f2 = X_function_g1_f2(THETA, PSI, B, Bm, coef_g1_f2, coefm_g1_f2, N, Nm);
             X_g1_f3 = X_function_g1_f3(THETA, PSI, B, Bm, coef_g1_f3, coefm_g1_f3, N, Nm);
             X_g2_f1 = X_function_g2_f1(THETA, PSI, B, Bm, coef_g2_f1, coefm_g2_f1, N, Nm);
             X_g2_f2 = X_function_g2_f2(THETA, PSI, B, Bm, coef_g2_f2, coefm_g2_f2, N, Nm);
             X_g2_f3 = X_function_g2_f3(THETA, PSI, B, Bm, coef_g2_f3, coefm_g2_f3, N, Nm);
             X_g3_f1 = X_function_g3_f1(THETA, PSI, B, Bm, coef_g3_f1, coefm_g3_f1, N, Nm);
             X_g3_f2 = X_function_g3_f2(THETA, PSI, B, Bm, coef_g3_f2, coefm_g3_f2, N, Nm);
             X_g3_f3 = X_function_g3_f3(THETA, PSI, B, Bm, coef_g3_f3, coefm_g3_f3, N, Nm);
             %
             I_psi_g1_f1  = I_psi_g1_f1  + w_psi(n_psi)*X_g1_f1;
             I_psi_g1_f2  = I_psi_g1_f2  + w_psi(n_psi)*X_g1_f2;
             I_psi_g1_f3  = I_psi_g1_f3  + w_psi(n_psi)*X_g1_f3;
             I_psi_g2_f1  = I_psi_g2_f1  + w_psi(n_psi)*X_g2_f1;
             I_psi_g2_f2  = I_psi_g2_f2  + w_psi(n_psi)*X_g2_f2;
             I_psi_g2_f3  = I_psi_g2_f3  + w_psi(n_psi)*X_g2_f3;
             I_psi_g3_f1  = I_psi_g3_f1  + w_psi(n_psi)*X_g3_f1;
             I_psi_g3_f2  = I_psi_g3_f2  + w_psi(n_psi)*X_g3_f2;
             I_psi_g3_f3  = I_psi_g3_f3  + w_psi(n_psi)*X_g3_f3;
         end
         I_psi_g1_f1  = J_psi * I_psi_g1_f1;
         I_theta_g1_f1  = I_theta_g1_f1 +w_theta(n_theta)*I_psi_g1_f1;  
         
         I_psi_g1_f2  = J_psi * I_psi_g1_f2;
         I_theta_g1_f2  = I_theta_g1_f2 +w_theta(n_theta)*I_psi_g1_f2; 
         
         I_psi_g1_f3  = J_psi * I_psi_g1_f3;
         I_theta_g1_f3  = I_theta_g1_f3 +w_theta(n_theta)*I_psi_g1_f3;  
         
         I_psi_g2_f1  = J_psi * I_psi_g2_f1;
         I_theta_g2_f1  = I_theta_g2_f1 +w_theta(n_theta)*I_psi_g2_f1;  
         
         I_psi_g2_f2  = J_psi * I_psi_g2_f2;
         I_theta_g2_f2  = I_theta_g2_f2 +w_theta(n_theta)*I_psi_g2_f2;   
         
         I_psi_g2_f3  = J_psi * I_psi_g2_f3;
         I_theta_g2_f3  = I_theta_g2_f3 +w_theta(n_theta)*I_psi_g2_f3;   
         
         I_psi_g3_f1  = J_psi * I_psi_g3_f1;
         I_theta_g3_f1  = I_theta_g3_f1 +w_theta(n_theta)*I_psi_g3_f1;   
         
         I_psi_g3_f2  = J_psi * I_psi_g3_f2;
         I_theta_g3_f2  = I_theta_g3_f2 +w_theta(n_theta)*I_psi_g3_f2;   
         
         I_psi_g3_f3  = J_psi * I_psi_g3_f3;
         I_theta_g3_f3  = I_theta_g3_f3 +w_theta(n_theta)*I_psi_g3_f3;   
     end
     I_g1_f1(m,1)  = J_theta*I_theta_g1_f1;
     I_g1_f2(m,1)  = J_theta*I_theta_g1_f2;
     I_g1_f3(m,1)  = J_theta*I_theta_g1_f3;
     I_g2_f1(m,1)  = J_theta*I_theta_g2_f1;
     I_g2_f2(m,1)  = J_theta*I_theta_g2_f2;
     I_g2_f3(m,1)  = J_theta*I_theta_g2_f3;
     I_g3_f1(m,1)  = J_theta*I_theta_g3_f1;
     I_g3_f2(m,1)  = J_theta*I_theta_g3_f2;
     I_g3_f3(m,1)  = J_theta*I_theta_g3_f3;
 end
 Ig1_f1  = sum(I_g1_f1,1);
 Ising_g1_f1 =  Jp*Jq*Ig1_f1;
 Ising_g1_f1 = Ising_g1_f1 * norm(r3-r2)*norm(r4-r1)/(4*Ap*Aq);
 
 Ig1_f2  = sum(I_g1_f2,1);
 Ising_g1_f2 =  Jp*Jq*Ig1_f2;
 Ising_g1_f2 = Ising_g1_f2 * norm(r3-r2)*norm(r4-r2)/(4*Ap*Aq);
 
 Ig1_f3  = sum(I_g1_f3,1);
 Ising_g1_f3 =  Jp*Jq*Ig1_f3;
 Ising_g1_f3 = Ising_g1_f3 * norm(r3-r2)*norm(r2-r1)/(4*Ap*Aq);
 
 Ig2_f1  = sum(I_g2_f1,1);
 Ising_g2_f1 =  Jp*Jq*Ig2_f1;
 Ising_g2_f1 = Ising_g2_f1 * norm(r3-r1)*norm(r4-r1)/(4*Ap*Aq);
 
 Ig2_f2  = sum(I_g2_f2,1);
 Ising_g2_f2 =  Jp*Jq*Ig2_f2;
 Ising_g2_f2 = Ising_g2_f2 * norm(r3-r1)*norm(r4-r2)/(4*Ap*Aq);
 
 Ig2_f3  = sum(I_g2_f3,1);
 Ising_g2_f3 =  Jp*Jq*Ig2_f3;
 Ising_g2_f3 = Ising_g2_f3 * norm(r3-r1)*norm(r2-r1)/(4*Ap*Aq);
 
 Ig3_f1  = sum(I_g3_f1,1);
 Ising_g3_f1 =  Jp*Jq*Ig3_f1;
 Ising_g3_f1 = Ising_g3_f1 * norm(r2-r1)*norm(r4-r1)/(4*Ap*Aq);
 
 Ig3_f2  = sum(I_g3_f2,1);
 Ising_g3_f2 =  Jp*Jq*Ig3_f2;
 Ising_g3_f2 = Ising_g3_f2 * norm(r2-r1)*norm(r4-r2)/(4*Ap*Aq);
 
 Ig3_f3  = sum(I_g3_f3,1);
 Ising_g3_f3 =  Jp*Jq*Ig3_f3;
 Ising_g3_f3 = Ising_g3_f3 * norm(r2-r1)*norm(r2-r1)/(4*Ap*Aq);
 
 % Final output
I_DE(1) = Ising_g1_f1;
I_DE(2) = Ising_g1_f2;
I_DE(3) = Ising_g1_f3;
I_DE(4) = Ising_g2_f1;
I_DE(5) = Ising_g2_f2;
I_DE(6) = Ising_g2_f3;
I_DE(7) = Ising_g3_f1;
I_DE(8) = Ising_g3_f2;
I_DE(9) = Ising_g3_f3;

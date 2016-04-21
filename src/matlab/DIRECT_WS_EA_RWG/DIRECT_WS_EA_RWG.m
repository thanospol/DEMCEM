function [I_DE] = DIRECT_WS_EA_RWG(r1,r2,r3,r4,N_theta,N_psi)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the edge adjacent 4-D weakly singular integrals over planar
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
% I_DE(10)  = scalar and potential 4-D weakly singular integrals
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

global ko

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = (r2-r1)/2;
beta  = (2*r3-r1-r2)/(2*sqrt(3));
gamma = -(2*r4-r1-r2)/(2*sqrt(3));
% 
[coef_const,coefm_const] = coefficients_const(ko);

[coef_f1_f1,coefm_f1_f1] = coefficients_f1_f1(r1,r2,r3,r4,ko);
[coef_f1_f2,coefm_f1_f2] = coefficients_f1_f2(r1,r2,r3,r4,ko);
[coef_f1_f3,coefm_f1_f3] = coefficients_f1_f3(r1,r2,r3,r4,ko);
[coef_f2_f1,coefm_f2_f1] = coefficients_f2_f1(r1,r2,r3,r4,ko);
[coef_f2_f2,coefm_f2_f2] = coefficients_f2_f2(r1,r2,r3,r4,ko);
[coef_f2_f3,coefm_f2_f3] = coefficients_f2_f3(r1,r2,r3,r4,ko);
[coef_f3_f1,coefm_f3_f1] = coefficients_f3_f1(r1,r2,r3,r4,ko);
[coef_f3_f2,coefm_f3_f2] = coefficients_f3_f2(r1,r2,r3,r4,ko);
[coef_f3_f3,coefm_f3_f3] = coefficients_f3_f3(r1,r2,r3,r4,ko);
%
[ w_theta, z_theta ] = Gauss_1D ( N_theta);
[ w_psi, z_psi ] = Gauss_1D ( N_psi );
%
I_const = zeros(6,1);
I_f1_f1 = zeros(6,1);
I_f1_f2 = zeros(6,1);
I_f1_f3 = zeros(6,1);
I_f2_f1 = zeros(6,1);
I_f2_f2 = zeros(6,1);
I_f2_f3 = zeros(6,1);
I_f3_f1 = zeros(6,1);
I_f3_f2 = zeros(6,1);
I_f3_f3 = zeros(6,1);
%  
 for m = 1:6
     I_theta_const = 0;
     I_theta_f1_f1 = 0;
     I_theta_f1_f2 = 0;
     I_theta_f1_f3 = 0;
     I_theta_f2_f1 = 0;
     I_theta_f2_f2 = 0;
     I_theta_f2_f3 = 0;
     I_theta_f3_f1 = 0;
     I_theta_f3_f2 = 0;
     I_theta_f3_f3 = 0;
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
         I_psi_const = 0;
         I_psi_f1_f1 = 0;
         I_psi_f1_f2 = 0;
         I_psi_f1_f3 = 0;
         I_psi_f2_f1 = 0;
         I_psi_f2_f2 = 0;
         I_psi_f2_f3 = 0;
         I_psi_f3_f1 = 0;
         I_psi_f3_f2 = 0;
         I_psi_f3_f3 = 0;
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
             X_const = X_function_const(PSI, B, Bm, coef_const, coefm_const, N, Nm);
             X_f1_f1 = X_function_f1_f1(THETA, PSI, B, Bm, coef_f1_f1, coefm_f1_f1, N, Nm);
             X_f1_f2 = X_function_f1_f2(THETA, PSI, B, Bm, coef_f1_f2, coefm_f1_f2, N, Nm);
             X_f1_f3 = X_function_f1_f3(THETA, PSI, B, Bm, coef_f1_f3, coefm_f1_f3, N, Nm);
             X_f2_f1 = X_function_f2_f1(THETA, PSI, B, Bm, coef_f2_f1, coefm_f2_f1, N, Nm);
             X_f2_f2 = X_function_f2_f2(THETA, PSI, B, Bm, coef_f2_f2, coefm_f2_f2, N, Nm);
             X_f2_f3 = X_function_f2_f3(THETA, PSI, B, Bm, coef_f2_f3, coefm_f2_f3, N, Nm);
             X_f3_f1 = X_function_f3_f1(THETA, PSI, B, Bm, coef_f3_f1, coefm_f3_f1, N, Nm);
             X_f3_f2 = X_function_f3_f2(THETA, PSI, B, Bm, coef_f3_f2, coefm_f3_f2, N, Nm);
             X_f3_f3 = X_function_f3_f3(THETA, PSI, B, Bm, coef_f3_f3, coefm_f3_f3, N, Nm);
             %
             I_psi_const  = I_psi_const  + w_psi(n_psi)*X_const;
             I_psi_f1_f1  = I_psi_f1_f1  + w_psi(n_psi)*X_f1_f1;
             I_psi_f1_f2  = I_psi_f1_f2  + w_psi(n_psi)*X_f1_f2;
             I_psi_f1_f3  = I_psi_f1_f3  + w_psi(n_psi)*X_f1_f3;
             I_psi_f2_f1  = I_psi_f2_f1  + w_psi(n_psi)*X_f2_f1;
             I_psi_f2_f2  = I_psi_f2_f2  + w_psi(n_psi)*X_f2_f2;
             I_psi_f2_f3  = I_psi_f2_f3  + w_psi(n_psi)*X_f2_f3;
             I_psi_f3_f1  = I_psi_f3_f1  + w_psi(n_psi)*X_f3_f1;
             I_psi_f3_f2  = I_psi_f3_f2  + w_psi(n_psi)*X_f3_f2;
             I_psi_f3_f3  = I_psi_f3_f3  + w_psi(n_psi)*X_f3_f3;
         end
         I_psi_const  = J_psi * I_psi_const;
         I_theta_const  = I_theta_const +w_theta(n_theta)*I_psi_const; 
         
         I_psi_f1_f1  = J_psi * I_psi_f1_f1;
         I_theta_f1_f1  = I_theta_f1_f1 +w_theta(n_theta)*I_psi_f1_f1;   

         I_psi_f1_f2  = J_psi * I_psi_f1_f2;
         I_theta_f1_f2  = I_theta_f1_f2 +w_theta(n_theta)*I_psi_f1_f2; 
         
         I_psi_f1_f3  = J_psi * I_psi_f1_f3;
         I_theta_f1_f3  = I_theta_f1_f3 +w_theta(n_theta)*I_psi_f1_f3;  
 
         I_psi_f2_f1  = J_psi * I_psi_f2_f1;
         I_theta_f2_f1  = I_theta_f2_f1 +w_theta(n_theta)*I_psi_f2_f1; 
         
         I_psi_f2_f2  = J_psi * I_psi_f2_f2;
         I_theta_f2_f2  = I_theta_f2_f2 +w_theta(n_theta)*I_psi_f2_f2;   
         
         I_psi_f2_f3  = J_psi * I_psi_f2_f3;
         I_theta_f2_f3  = I_theta_f2_f3 +w_theta(n_theta)*I_psi_f2_f3;   
         
         I_psi_f3_f1  = J_psi * I_psi_f3_f1;
         I_theta_f3_f1  = I_theta_f3_f1 +w_theta(n_theta)*I_psi_f3_f1;   
         
         I_psi_f3_f2  = J_psi * I_psi_f3_f2;
         I_theta_f3_f2  = I_theta_f3_f2 +w_theta(n_theta)*I_psi_f3_f2;   
         
         I_psi_f3_f3  = J_psi * I_psi_f3_f3;
         I_theta_f3_f3  = I_theta_f3_f3 +w_theta(n_theta)*I_psi_f3_f3;   
     end
     I_const(m,1)  = J_theta*I_theta_const;
     I_f1_f1(m,1)  = J_theta*I_theta_f1_f1;
     I_f1_f2(m,1)  = J_theta*I_theta_f1_f2;
     I_f1_f3(m,1)  = J_theta*I_theta_f1_f3;
     I_f2_f1(m,1)  = J_theta*I_theta_f2_f1;
     I_f2_f2(m,1)  = J_theta*I_theta_f2_f2;
     I_f2_f3(m,1)  = J_theta*I_theta_f2_f3;
     I_f3_f1(m,1)  = J_theta*I_theta_f3_f1;
     I_f3_f2(m,1)  = J_theta*I_theta_f3_f2;
     I_f3_f3(m,1)  = J_theta*I_theta_f3_f3;
 end
 Iconst  = sum(I_const,1);
 
 If1_f1  = sum(I_f1_f1,1);

 If1_f2  = sum(I_f1_f2,1);
 
 If1_f3  = sum(I_f1_f3,1);
 
 If2_f1  = sum(I_f2_f1,1);
 
 If2_f2  = sum(I_f2_f2,1);
 
 If2_f3  = sum(I_f2_f3,1);
 
 If3_f1  = sum(I_f3_f1,1);
 
 If3_f2  = sum(I_f3_f2,1);
 
 If3_f3  = sum(I_f3_f3,1);
 
% Final output
I_DE(1) = norm(r3-r2)*norm(r4-r1)/(12) * (1i * ko * If1_f1 + 4 / (1i *ko) * Iconst);
I_DE(2) = norm(r3-r2)*norm(r4-r2)/(12) * (1i * ko * If1_f2 + 4 / (1i *ko) * Iconst);
I_DE(3) = norm(r3-r2)*norm(r2-r1)/(12) * (1i * ko * If1_f3 + 4 / (1i *ko) * Iconst);
I_DE(4) = norm(r3-r1)*norm(r4-r1)/(12) * (1i * ko * If2_f1 + 4 / (1i *ko) * Iconst);
I_DE(5) = norm(r3-r1)*norm(r4-r2)/(12) * (1i * ko * If2_f2 + 4 / (1i *ko) * Iconst);
I_DE(6) = norm(r3-r1)*norm(r2-r1)/(12) * (1i * ko * If2_f3 + 4 / (1i *ko) * Iconst);
I_DE(7) = norm(r2-r1)*norm(r4-r1)/(12) * (1i * ko * If3_f1 + 4 / (1i *ko) * Iconst);
I_DE(8) = norm(r2-r1)*norm(r4-r2)/(12) * (1i * ko * If3_f2 + 4 / (1i *ko) * Iconst);
I_DE(9) = norm(r2-r1)*norm(r2-r1)/(12) * (1i * ko * If3_f3 + 4 / (1i *ko) * Iconst);
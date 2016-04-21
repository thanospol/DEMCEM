function [I_DE] = DIRECT_SS_VA_nxRWG(r1,r2,r3,r4,r5,N_theta_p,N_theta_q,N_psi)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the vertex adjacent 4-D strongly singular integrals over planar
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
% r1,r2,r3,r4,r5 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r1,r4,r5)
% N_theta_p = order of the GL quadrature for the 1-D smooth integral over theta_p
% N_theta_q = order of the GL quadrature for the 1-D smooth integral over theta_q
% N_psi     = order of the GL quadrature for the 1-D smooth integral over Psi

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
rq1 = r1;
rq2 = r4;
rq3 = r5;
%%%%%%%%%%%%%%%%%%  Triangle Area   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aq = (1/2)*norm(cross(rq2-rq1,rq3-rq1));
%%%%%%%%%%%%%%%%%%  Jacobian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
Jq = Aq/sqrt(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_vo  = (r2-r1)/2;
alpha_v1  = (2*r3-r1-r2)/(2*sqrt(3));
alpha_v2  = (r4-r1)/2;
alpha_v3  = (2*r5-r1-r4)/(2*sqrt(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[coef_g1_f1] = coefficients_g1_f1(r1,r2,r3,r4,r5,Ap);
[coef_g1_f2] = coefficients_g1_f2(r1,r2,r3,r4,r5,Ap);
[coef_g1_f3] = coefficients_g1_f3(r1,r2,r3,r4,r5,Ap);
[coef_g2_f1] = coefficients_g2_f1(r1,r2,r3,r4,r5,Ap);
[coef_g2_f2] = coefficients_g2_f2(r1,r2,r3,r4,r5,Ap);
[coef_g2_f3] = coefficients_g2_f3(r1,r2,r3,r4,r5,Ap);
[coef_g3_f1] = coefficients_g3_f1(r1,r2,r3,r4,r5,Ap);
[coef_g3_f2] = coefficients_g3_f2(r1,r2,r3,r4,r5,Ap);
[coef_g3_f3] = coefficients_g3_f3(r1,r2,r3,r4,r5,Ap);

%-------------------------------------------------------------------------%
%                  3D  Gauss-Legendre Quadrature Rule                     %
%-------------------------------------------------------------------------%
[w_theta_p,z_theta_p] = Gauss_1D(N_theta_p);
[w_theta_q,z_theta_q] = Gauss_1D(N_theta_q);
[w_psi,z_psi]         = Gauss_1D(N_psi);
%---------------------------------------------
%%%%%%%%%%%%%%%%%%%%% Gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    I_theta_p_1_1 = 0;
    I_theta_p_1_2 = 0; 
    I_theta_p_1_3 = 0; 
    I_theta_p_2_1 = 0; 
    I_theta_p_2_2 = 0; 
    I_theta_p_2_3 = 0; 
    I_theta_p_3_1 = 0; 
    I_theta_p_3_2 = 0;
    I_theta_p_3_3 = 0;

        for n_theta_p = 1:N_theta_p
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Int_1, theta_p_A=< theta_p 0 <= theta_p_A  %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta_p_A = 0;
        theta_p_B = pi/3;
        %
        THETA_p = ((theta_p_B-theta_p_A)/2)*z_theta_p(n_theta_p)+(theta_p_B+theta_p_A)/2;
        %
        J_theta_p = (theta_p_B-theta_p_A)/2;
        %
        Lp = (2*sqrt(3))/(sin(THETA_p)+sqrt(3)*cos(THETA_p));
        %
        b_v0 = alpha_vo*cos(THETA_p) + alpha_v1*sin(THETA_p);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    I_theta_q_1_1 = 0;
    I_theta_q_1_2 = 0;
    I_theta_q_1_3 = 0;
    I_theta_q_2_1 = 0;
    I_theta_q_2_2 = 0;
    I_theta_q_2_3 = 0;
    I_theta_q_3_1 = 0;
    I_theta_q_3_2 = 0;
    I_theta_q_3_3 = 0;

        for n_theta_q = 1:N_theta_q
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% Int_1, theta_q_A=< theta_q 0 <= theta_q_A  %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta_q_A = 0;
        theta_q_B = pi/3;
        %
        THETA_q = ((theta_q_B-theta_q_A)/2)*z_theta_q(n_theta_q)+(theta_q_B+theta_q_A)/2;
        %
        J_theta_q = (theta_q_B-theta_q_A)/2;
        %
        Lq = (2*sqrt(3))/(sin(THETA_q)+sqrt(3)*cos(THETA_q));
        %
        b_v1 = alpha_v2*cos(THETA_q) + alpha_v3*sin(THETA_q);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_v = norm(b_v0)^2;
        b_v = -2*dot(b_v0,b_v1);
        c_v = norm(b_v1)^2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        I_psi_1_1 = 0;
        I_psi_1_2 = 0;
        I_psi_1_3 = 0;
        I_psi_2_1 = 0;
        I_psi_2_2 = 0;
        I_psi_2_3 = 0;
        I_psi_3_1 = 0;
        I_psi_3_2 = 0;
        I_psi_3_3 = 0;
   
        for n_psi = 1:N_psi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% psi_A =< PSI <= psi_B %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             psi_1_A = 0;
             psi_1_B = atan(Lq/Lp);
             %
             PSI_1 = ((psi_1_B-psi_1_A)/2)*z_psi(n_psi)+(psi_1_B+psi_1_A)/2;
             %
             J_psi_1 = (psi_1_B-psi_1_A)/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             psi_2_A = atan(Lq/Lp);
             psi_2_B = pi/2;
             %
             PSI_2 = ((psi_2_B-psi_2_A)/2)*z_psi(n_psi)+(psi_2_B+psi_2_A)/2;
             %
             J_psi_2 = (psi_2_B-psi_2_A)/2;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %
             B_1 =sqrt(a_v*cos(PSI_1)^2+b_v*cos(PSI_1)*sin(PSI_1)+c_v*sin(PSI_1)^2);
             %
             %
             B_2 =sqrt(a_v*cos(PSI_2)^2+b_v*cos(PSI_2)*sin(PSI_2)+c_v*sin(PSI_2)^2);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %
             Lp = (2*sqrt(3))/(sin(THETA_p)+sqrt(3)*cos(THETA_p));
             L1 = Lp/cos(PSI_1);

             Lq = (2*sqrt(3))/(sin(THETA_q)+sqrt(3)*cos(THETA_q));
             L2 = Lq/sin(PSI_2);
             %
             K_1 = K_functions(B_1,L1);
             
             K_2 = K_functions(B_2,L2);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             Omega_g1_f1_1 = Omega_function_g1_f1(THETA_p,THETA_q, PSI_1, B_1, coef_g1_f1, K_1);

             Omega_g1_f1_2 = Omega_function_g1_f1(THETA_p,THETA_q, PSI_2, B_2, coef_g1_f1, K_2);
             %
             Omega_g1_f2_1 = Omega_function_g1_f2(THETA_p,THETA_q, PSI_1, B_1, coef_g1_f2, K_1);

             Omega_g1_f2_2 = Omega_function_g1_f2(THETA_p,THETA_q, PSI_2, B_2, coef_g1_f2, K_2);
             %           
             Omega_g1_f3_1 = Omega_function_g1_f3(THETA_p,THETA_q, PSI_1, B_1, coef_g1_f3, K_1);

             Omega_g1_f3_2 = Omega_function_g1_f3(THETA_p,THETA_q, PSI_2, B_2, coef_g1_f3, K_2);
             %           
             Omega_g2_f1_1 = Omega_function_g2_f1(THETA_p,THETA_q, PSI_1, B_1, coef_g2_f1, K_1);

             Omega_g2_f1_2 = Omega_function_g2_f1(THETA_p,THETA_q, PSI_2, B_2, coef_g2_f1, K_2);
             %           
             Omega_g2_f2_1 = Omega_function_g2_f2(THETA_p,THETA_q, PSI_1, B_1, coef_g2_f2, K_1);

             Omega_g2_f2_2 = Omega_function_g2_f2(THETA_p,THETA_q, PSI_2, B_2, coef_g2_f2, K_2);
             %           
             Omega_g2_f3_1 = Omega_function_g2_f3(THETA_p,THETA_q, PSI_1, B_1, coef_g2_f3, K_1);

             Omega_g2_f3_2 = Omega_function_g2_f3(THETA_p,THETA_q, PSI_2, B_2, coef_g2_f3, K_2);
             %           
             Omega_g3_f1_1 = Omega_function_g3_f1(THETA_p,THETA_q, PSI_1, B_1, coef_g3_f1, K_1);

             Omega_g3_f1_2 = Omega_function_g3_f1(THETA_p,THETA_q, PSI_2, B_2, coef_g3_f1, K_2);
             %           
             Omega_g3_f2_1 = Omega_function_g3_f2(THETA_p,THETA_q, PSI_1, B_1, coef_g3_f2, K_1);

             Omega_g3_f2_2 = Omega_function_g3_f2(THETA_p,THETA_q, PSI_2, B_2, coef_g3_f2, K_2);
             %           
             Omega_g3_f3_1 = Omega_function_g3_f3(THETA_p,THETA_q, PSI_1, B_1, coef_g3_f3, K_1);

             Omega_g3_f3_2 = Omega_function_g3_f3(THETA_p,THETA_q, PSI_2, B_2, coef_g3_f3, K_2);
   
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             I_psi_1_1  = I_psi_1_1  + w_psi(n_psi)*(J_psi_1*Omega_g1_f1_1 + J_psi_2*Omega_g1_f1_2);
             %
             I_psi_1_2  = I_psi_1_2  + w_psi(n_psi)*(J_psi_1*Omega_g1_f2_1 + J_psi_2*Omega_g1_f2_2);
             %
             I_psi_1_3  = I_psi_1_3  + w_psi(n_psi)*(J_psi_1*Omega_g1_f3_1 + J_psi_2*Omega_g1_f3_2);
             %
             I_psi_2_1  = I_psi_2_1  + w_psi(n_psi)*(J_psi_1*Omega_g2_f1_1 + J_psi_2*Omega_g2_f1_2);
             %
             I_psi_2_2  = I_psi_2_2  + w_psi(n_psi)*(J_psi_1*Omega_g2_f2_1 + J_psi_2*Omega_g2_f2_2);
             %
             I_psi_2_3  = I_psi_2_3  + w_psi(n_psi)*(J_psi_1*Omega_g2_f3_1 + J_psi_2*Omega_g2_f3_2);
             %
             I_psi_3_1  = I_psi_3_1  + w_psi(n_psi)*(J_psi_1*Omega_g3_f1_1 + J_psi_2*Omega_g3_f1_2);
             %
             I_psi_3_2  = I_psi_3_2  + w_psi(n_psi)*(J_psi_1*Omega_g3_f2_1 + J_psi_2*Omega_g3_f2_2);
             %
             I_psi_3_3  = I_psi_3_3  + w_psi(n_psi)*(J_psi_1*Omega_g3_f3_1 + J_psi_2*Omega_g3_f3_2);

             
        end  %for n_psi = 1:N_psi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_theta_q_1_1  = I_theta_q_1_1 +w_theta_q(n_theta_q)*I_psi_1_1;
    %
    I_theta_q_1_2  = I_theta_q_1_2 +w_theta_q(n_theta_q)*I_psi_1_2;
    %
    I_theta_q_1_3  = I_theta_q_1_3 +w_theta_q(n_theta_q)*I_psi_1_3;
    %
    I_theta_q_2_1  = I_theta_q_2_1 +w_theta_q(n_theta_q)*I_psi_2_1;
    %
    I_theta_q_2_2  = I_theta_q_2_2 +w_theta_q(n_theta_q)*I_psi_2_2;
    %
    I_theta_q_2_3  = I_theta_q_2_3 +w_theta_q(n_theta_q)*I_psi_2_3;
    %
    I_theta_q_3_1  = I_theta_q_3_1 +w_theta_q(n_theta_q)*I_psi_3_1;
    %
    I_theta_q_3_2  = I_theta_q_3_2 +w_theta_q(n_theta_q)*I_psi_3_2;
    %
    I_theta_q_3_3  = I_theta_q_3_3 +w_theta_q(n_theta_q)*I_psi_3_3;

        end  %for n_theta_q = 1:N_theta_q
    
    I_theta_q_1_1 = J_theta_q*I_theta_q_1_1;
    %    
    I_theta_q_1_2 = J_theta_q*I_theta_q_1_2;
    %
    I_theta_q_1_3 = J_theta_q*I_theta_q_1_3;
    %
    I_theta_q_2_1 = J_theta_q*I_theta_q_2_1;
    %
    I_theta_q_2_2 = J_theta_q*I_theta_q_2_2;
    %
    I_theta_q_2_3 = J_theta_q*I_theta_q_2_3;
    %
    I_theta_q_3_1 = J_theta_q*I_theta_q_3_1;
    %
    I_theta_q_3_2 = J_theta_q*I_theta_q_3_2;
    %
    I_theta_q_3_3 = J_theta_q*I_theta_q_3_3;

    %
    I_theta_p_1_1  = I_theta_p_1_1 +w_theta_p(n_theta_p)*I_theta_q_1_1;
    %
    I_theta_p_1_2  = I_theta_p_1_2 +w_theta_p(n_theta_p)*I_theta_q_1_2;
    %
    I_theta_p_1_3  = I_theta_p_1_3 +w_theta_p(n_theta_p)*I_theta_q_1_3;
    %
    I_theta_p_2_1  = I_theta_p_2_1 +w_theta_p(n_theta_p)*I_theta_q_2_1;
    %
    I_theta_p_2_2  = I_theta_p_2_2 +w_theta_p(n_theta_p)*I_theta_q_2_2;
    %
    I_theta_p_2_3  = I_theta_p_2_3 +w_theta_p(n_theta_p)*I_theta_q_2_3;
    %
    I_theta_p_3_1  = I_theta_p_3_1 +w_theta_p(n_theta_p)*I_theta_q_3_1;
    %
    I_theta_p_3_2  = I_theta_p_3_2 +w_theta_p(n_theta_p)*I_theta_q_3_2;
    %
    I_theta_p_3_3  = I_theta_p_3_3 +w_theta_p(n_theta_p)*I_theta_q_3_3;

        end  %for n_theta_p = 1:N_theta_p
    
    I_theta_p_1_1 = J_theta_p*I_theta_p_1_1;
    %
    I_theta_p_1_2 = J_theta_p*I_theta_p_1_2;
    %
    I_theta_p_1_3 = J_theta_p*I_theta_p_1_3;
    %
    I_theta_p_2_1 = J_theta_p*I_theta_p_2_1;
    %
    I_theta_p_2_2 = J_theta_p*I_theta_p_2_2;
    %
    I_theta_p_2_3 = J_theta_p*I_theta_p_2_3;
    %
    I_theta_p_3_1 = J_theta_p*I_theta_p_3_1;
    %
    I_theta_p_3_2 = J_theta_p*I_theta_p_3_2;
    %
    I_theta_p_3_3 = J_theta_p*I_theta_p_3_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Output 
I_DE(1)  = Jp * Jq * norm(r3-r2)*norm(r5-r4) / (4*Ap*Aq) * I_theta_p_1_1;
I_DE(2)  = Jp * Jq * norm(r3-r2)*norm(r5-r1) / (4*Ap*Aq) * I_theta_p_1_2;
I_DE(3)  = Jp * Jq * norm(r3-r2)*norm(r4-r1) / (4*Ap*Aq) * I_theta_p_1_3;
I_DE(4)  = Jp * Jq * norm(r3-r1)*norm(r5-r4) / (4*Ap*Aq) * I_theta_p_2_1;
I_DE(5)  = Jp * Jq * norm(r3-r1)*norm(r5-r1) / (4*Ap*Aq) * I_theta_p_2_2;
I_DE(6)  = Jp * Jq * norm(r3-r1)*norm(r4-r1) / (4*Ap*Aq) * I_theta_p_2_3;
I_DE(7)  = Jp * Jq * norm(r2-r1)*norm(r5-r4) / (4*Ap*Aq) * I_theta_p_3_1;
I_DE(8)  = Jp * Jq * norm(r2-r1)*norm(r5-r1) / (4*Ap*Aq) * I_theta_p_3_2;
I_DE(9)  = Jp * Jq * norm(r2-r1)*norm(r4-r1) / (4*Ap*Aq) * I_theta_p_3_3;
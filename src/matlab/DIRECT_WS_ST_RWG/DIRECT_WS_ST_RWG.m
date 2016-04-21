function [I_DE] = DIRECT_WS_ST_RWG(r1,r2,r3,Np_1D)
%% Main body of the DIRECT EVALUATION method for the
% evaluation of the coincident 4-D weakly singular integrals over planar
% triangular elements.

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
% r1,r2,r3 = point vectors of the triangular element's vertices
% Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
% Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r3)
% Np_1D = order of the Gauss-Legendre quadrature rule

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

global ko
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Asing0 = norm(r1-r2)^2/4;
Asing1 = dot(r1-r2,r1+r2-2*r3)/(2*sqrt(3));
Asing2 = norm(r1+r2-2*r3)^2/12; 
%
Asing3 = (1/4)        *Asing0   -(sqrt(3)/4)  *Asing1      +(3/4)        *Asing2;
Asing4 = (sqrt(3)/2)  *Asing0   -(1/2)        *Asing1      -(sqrt(3)/2)  *Asing2;
Asing5 = (3/4)        *Asing0   +(sqrt(3)/4)  *Asing1      +(1/4)        *Asing2;
%
Asing6 = (1/4)        *Asing0   +(sqrt(3)/4)  *Asing1       +(3/4)       *Asing2;
Asing7 = -(sqrt(3)/2) *Asing0   -(1/2)        *Asing1       +(sqrt(3)/2) *Asing2;
Asing8 = (3/4)        *Asing0   -(sqrt(3)/4)  *Asing1       +(1/4)       *Asing2;
%
Asing = [Asing0 Asing3 Asing6
         Asing1 Asing4 Asing7
         Asing2 Asing5 Asing8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  1D  Gauss-Legendre Quadrature Rule                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w,z] = Gauss_1D(Np_1D);

Isub = zeros(3,3,3);
%---------------------------------------------
for m = 1:3
    %%%%%%%%% Define the coefs of the appropriate subtriangle
    acc = Asing(1,m);
    acs = Asing(2,m);
    ass = Asing(3,m);
    %%%%%%%%%%%%%%%%%%%%% Gauss quadrature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Int1_1_1=0;Int2_1_1=0;Int3_1_1=0;
    Int1_1_2=0;Int2_1_2=0;Int3_1_2=0;
    Int1_1_3=0;Int2_1_3=0;Int3_1_3=0;
    Int1_2_1=0;Int2_2_1=0;Int3_2_1=0;
    Int1_2_2=0;Int2_2_2=0;Int3_2_2=0;
    Int1_2_3=0;Int2_2_3=0;Int3_2_3=0;
    Int1_3_1=0;Int2_3_1=0;Int3_3_1=0;
    Int1_3_2=0;Int2_3_2=0;Int3_3_2=0;
    Int1_3_3=0;Int2_3_3=0;Int3_3_3=0;
    for kk = 1:Np_1D
        %%%%%%%%%%%%%%% Int1,  0 =< PSI <= pi/3 %%%%%%%%%%%%%%%%%%%%%%%%%%%
        PSI_1a = 0;
        PSI_1b = pi/3;
        PSI_1k = ((PSI_1b-PSI_1a)*z(kk)+(PSI_1b+PSI_1a))/2;
        a_PSI_1k = sqrt(acc*cos(PSI_1k)^2-acs*cos(PSI_1k)*sin(PSI_1k)+ass*sin(PSI_1k)^2);
        % Int1_1_1 
        PHI_a_1_1 = Phi_1_1(1,PSI_1k,a_PSI_1k);
        PHI_e_1_1 = Phi_1_1(5,PSI_1k,a_PSI_1k);
        PHI_f_1_1 = Phi_1_1(6,PSI_1k,a_PSI_1k);
        
        F1_1_1 = PHI_a_1_1+PHI_e_1_1+PHI_f_1_1;       
        Int1_1_1  = Int1_1_1+w(kk)*F1_1_1;
        % Int1_1_2 
        PHI_a_1_2 = Phi_1_2(1,PSI_1k,a_PSI_1k);
        PHI_e_1_2 = Phi_1_2(5,PSI_1k,a_PSI_1k);
        PHI_f_1_2 = Phi_1_2(6,PSI_1k,a_PSI_1k);
        
        F1_1_2 = PHI_a_1_2+PHI_e_1_2+PHI_f_1_2;       
        Int1_1_2  = Int1_1_2+w(kk)*F1_1_2;
        % Int1_1_3 
        PHI_a_1_3 = Phi_1_3(1,PSI_1k,a_PSI_1k);
        PHI_e_1_3 = Phi_1_3(5,PSI_1k,a_PSI_1k);
        PHI_f_1_3 = Phi_1_3(6,PSI_1k,a_PSI_1k);
        
        F1_1_3 = PHI_a_1_3+PHI_e_1_3+PHI_f_1_3;       
        Int1_1_3  = Int1_1_3+w(kk)*F1_1_3;
        % Int1_2_1 
        PHI_a_2_1 = Phi_2_1(1,PSI_1k,a_PSI_1k);
        PHI_e_2_1 = Phi_2_1(5,PSI_1k,a_PSI_1k);
        PHI_f_2_1 = Phi_2_1(6,PSI_1k,a_PSI_1k);
        
        F1_2_1 = PHI_a_2_1+PHI_e_2_1+PHI_f_2_1;       
        Int1_2_1  = Int1_2_1+w(kk)*F1_2_1;
        % Int1_2_2 
        PHI_a_2_2 = Phi_2_2(1,PSI_1k,a_PSI_1k);
        PHI_e_2_2 = Phi_2_2(5,PSI_1k,a_PSI_1k);
        PHI_f_2_2 = Phi_2_2(6,PSI_1k,a_PSI_1k);
        
        F1_2_2 = PHI_a_2_2+PHI_e_2_2+PHI_f_2_2;       
        Int1_2_2  = Int1_2_2+w(kk)*F1_2_2;
        % Int1_2_3 
        PHI_a_2_3 = Phi_2_3(1,PSI_1k,a_PSI_1k);
        PHI_e_2_3 = Phi_2_3(5,PSI_1k,a_PSI_1k);
        PHI_f_2_3 = Phi_2_3(6,PSI_1k,a_PSI_1k);
        
        F1_2_3 = PHI_a_2_3+PHI_e_2_3+PHI_f_2_3;       
        Int1_2_3  = Int1_2_3+w(kk)*F1_2_3;
        % Int1_3_1 
        PHI_a_3_1 = Phi_3_1(1,PSI_1k,a_PSI_1k);
        PHI_e_3_1 = Phi_3_1(5,PSI_1k,a_PSI_1k);
        PHI_f_3_1 = Phi_3_1(6,PSI_1k,a_PSI_1k);
        
        F1_3_1 = PHI_a_3_1+PHI_e_3_1+PHI_f_3_1;       
        Int1_3_1  = Int1_3_1+w(kk)*F1_3_1;
        % Int1_3_2 
        PHI_a_3_2 = Phi_3_2(1,PSI_1k,a_PSI_1k);
        PHI_e_3_2 = Phi_3_2(5,PSI_1k,a_PSI_1k);
        PHI_f_3_2 = Phi_3_2(6,PSI_1k,a_PSI_1k);
        
        F1_3_2 = PHI_a_3_2+PHI_e_3_2+PHI_f_3_2;       
        Int1_3_2  = Int1_3_2+w(kk)*F1_3_2;
        % Int1_3_3 
        PHI_a_3_3 = Phi_3_3(1,PSI_1k,a_PSI_1k);
        PHI_e_3_3 = Phi_3_3(5,PSI_1k,a_PSI_1k);
        PHI_f_3_3 = Phi_3_3(6,PSI_1k,a_PSI_1k);
        
        F1_3_3 = PHI_a_3_3+PHI_e_3_3+PHI_f_3_3;       
        Int1_3_3  = Int1_3_3+w(kk)*F1_3_3;
        %%%%%%%%%%%%%%% Int2,  pi/3 =< PSI <= 2pi/3 %%%%%%%%%%%%%%%%%%%%%%%
        PSI_2a = pi/3;
        PSI_2b = 2*pi/3;
        PSI_2k = ((PSI_2b-PSI_2a)*z(kk)+(PSI_2b+PSI_2a))/2;
        a_PSI_2k = sqrt(acc*cos(PSI_2k)^2-acs*cos(PSI_2k)*sin(PSI_2k)+ass*sin(PSI_2k)^2);
        % Int2_1_1 
        PHI_b_1_1 = Phi_1_1(2,PSI_2k,a_PSI_2k);
        PHI_g_1_1 = Phi_1_1(7,PSI_2k,a_PSI_2k);
                
        F2_1_1 = PHI_b_1_1+PHI_g_1_1;       
        Int2_1_1  = Int2_1_1+w(kk)*F2_1_1;
        % Int2_1_2 
        PHI_b_1_2 = Phi_1_2(2,PSI_2k,a_PSI_2k);
        PHI_g_1_2 = Phi_1_2(7,PSI_2k,a_PSI_2k);
                
        F2_1_2 = PHI_b_1_2+PHI_g_1_2;       
        Int2_1_2  = Int2_1_2+w(kk)*F2_1_2;
        % Int2_1_3 
        PHI_b_1_3 = Phi_1_3(2,PSI_2k,a_PSI_2k);
        PHI_g_1_3 = Phi_1_3(7,PSI_2k,a_PSI_2k);
                
        F2_1_3 = PHI_b_1_3+PHI_g_1_3;       
        Int2_1_3  = Int2_1_3+w(kk)*F2_1_3;
        % Int2_2_1 
        PHI_b_2_1 = Phi_2_1(2,PSI_2k,a_PSI_2k);
        PHI_g_2_1 = Phi_2_1(7,PSI_2k,a_PSI_2k);
                
        F2_2_1 = PHI_b_2_1+PHI_g_2_1;       
        Int2_2_1  = Int2_2_1+w(kk)*F2_2_1;
        % Int2_2_2 
        PHI_b_2_2 = Phi_2_2(2,PSI_2k,a_PSI_2k);
        PHI_g_2_2 = Phi_2_2(7,PSI_2k,a_PSI_2k);
                
        F2_2_2 = PHI_b_2_2+PHI_g_2_2;       
        Int2_2_2  = Int2_2_2+w(kk)*F2_2_2;
        % Int2_2_3 
        PHI_b_2_3 = Phi_2_3(2,PSI_2k,a_PSI_2k);
        PHI_g_2_3 = Phi_2_3(7,PSI_2k,a_PSI_2k);
                
        F2_2_3 = PHI_b_2_3+PHI_g_2_3;       
        Int2_2_3  = Int2_2_3+w(kk)*F2_2_3;
        % Int2_3_1 
        PHI_b_3_1 = Phi_3_1(2,PSI_2k,a_PSI_2k);
        PHI_g_3_1 = Phi_3_1(7,PSI_2k,a_PSI_2k);
                
        F2_3_1 = PHI_b_3_1+PHI_g_3_1;       
        Int2_3_1  = Int2_3_1+w(kk)*F2_3_1;
        % Int2_3_2 
        PHI_b_3_2 = Phi_3_2(2,PSI_2k,a_PSI_2k);
        PHI_g_3_2 = Phi_3_2(7,PSI_2k,a_PSI_2k);
                
        F2_3_2 = PHI_b_3_2+PHI_g_3_2;       
        Int2_3_2  = Int2_3_2+w(kk)*F2_3_2;
        % Int2_3_3 
        PHI_b_3_3 = Phi_3_3(2,PSI_2k,a_PSI_2k);
        PHI_g_3_3 = Phi_3_3(7,PSI_2k,a_PSI_2k);
                
        F2_3_3 = PHI_b_3_3+PHI_g_3_3;       
        Int2_3_3  = Int2_3_3+w(kk)*F2_3_3;
        %%%%%%%%%%%%%%% Int3,  2pi/3 =< PSI <= pi %%%%%%%%%%%%%%%%%%%%%%%%%
        PSI_3a = 2*pi/3;
        PSI_3b = pi;
        PSI_3k = ((PSI_3b-PSI_3a)*z(kk)+(PSI_3b+PSI_3a))/2;
        a_PSI_3k = sqrt(acc*cos(PSI_3k)^2-acs*cos(PSI_3k)*sin(PSI_3k)+ass*sin(PSI_3k)^2);
        % Int3_1_1 
        PHI_c_1_1 = Phi_1_1(3,PSI_3k,a_PSI_3k);
        PHI_d_1_1 = Phi_1_1(4,PSI_3k,a_PSI_3k);
        PHI_h_1_1 = Phi_1_1(8,PSI_3k,a_PSI_3k);
                
        F3_1_1 = PHI_c_1_1+PHI_d_1_1+PHI_h_1_1;       
        Int3_1_1  = Int3_1_1+w(kk)*F3_1_1;
        % Int3_1_2 
        PHI_c_1_2 = Phi_1_2(3,PSI_3k,a_PSI_3k);
        PHI_d_1_2 = Phi_1_2(4,PSI_3k,a_PSI_3k);
        PHI_h_1_2 = Phi_1_2(8,PSI_3k,a_PSI_3k);
                
        F3_1_2 = PHI_c_1_2+PHI_d_1_2+PHI_h_1_2;       
        Int3_1_2  = Int3_1_2+w(kk)*F3_1_2;
        % Int3_1_3 
        PHI_c_1_3 = Phi_1_3(3,PSI_3k,a_PSI_3k);
        PHI_d_1_3 = Phi_1_3(4,PSI_3k,a_PSI_3k);
        PHI_h_1_3 = Phi_1_3(8,PSI_3k,a_PSI_3k);
                
        F3_1_3 = PHI_c_1_3+PHI_d_1_3+PHI_h_1_3;       
        Int3_1_3  = Int3_1_3+w(kk)*F3_1_3;
        % Int3_2_1 
        PHI_c_2_1 = Phi_2_1(3,PSI_3k,a_PSI_3k);
        PHI_d_2_1 = Phi_2_1(4,PSI_3k,a_PSI_3k);
        PHI_h_2_1 = Phi_2_1(8,PSI_3k,a_PSI_3k);
                
        F3_2_1 = PHI_c_2_1+PHI_d_2_1+PHI_h_2_1;       
        Int3_2_1  = Int3_2_1+w(kk)*F3_2_1;
        % Int3_2_2 
        PHI_c_2_2 = Phi_2_2(3,PSI_3k,a_PSI_3k);
        PHI_d_2_2 = Phi_2_2(4,PSI_3k,a_PSI_3k);
        PHI_h_2_2 = Phi_2_2(8,PSI_3k,a_PSI_3k);
                
        F3_2_2 = PHI_c_2_2+PHI_d_2_2+PHI_h_2_2;       
        Int3_2_2  = Int3_2_2+w(kk)*F3_2_2;
        % Int3_2_3 
        PHI_c_2_3 = Phi_2_3(3,PSI_3k,a_PSI_3k);
        PHI_d_2_3 = Phi_2_3(4,PSI_3k,a_PSI_3k);
        PHI_h_2_3 = Phi_2_3(8,PSI_3k,a_PSI_3k);
                
        F3_2_3 = PHI_c_2_3+PHI_d_2_3+PHI_h_2_3;       
        Int3_2_3  = Int3_2_3+w(kk)*F3_2_3;
        % Int3_3_1 
        PHI_c_3_1 = Phi_3_1(3,PSI_3k,a_PSI_3k);
        PHI_d_3_1 = Phi_3_1(4,PSI_3k,a_PSI_3k);
        PHI_h_3_1 = Phi_3_1(8,PSI_3k,a_PSI_3k);
                
        F3_3_1 = PHI_c_3_1+PHI_d_3_1+PHI_h_3_1;       
        Int3_3_1  = Int3_3_1+w(kk)*F3_3_1;
        % Int3_3_2 
        PHI_c_3_2 = Phi_3_2(3,PSI_3k,a_PSI_3k);
        PHI_d_3_2 = Phi_3_2(4,PSI_3k,a_PSI_3k);
        PHI_h_3_2 = Phi_3_2(8,PSI_3k,a_PSI_3k);
                
        F3_3_2 = PHI_c_3_2+PHI_d_3_2+PHI_h_3_2;       
        Int3_3_2  = Int3_3_2+w(kk)*F3_3_2;
        % Int3_3_3 
        PHI_c_3_3 = Phi_3_3(3,PSI_3k,a_PSI_3k);
        PHI_d_3_3 = Phi_3_3(4,PSI_3k,a_PSI_3k);
        PHI_h_3_3 = Phi_3_3(8,PSI_3k,a_PSI_3k);
                
        F3_3_3 = PHI_c_3_3+PHI_d_3_3+PHI_h_3_3;       
        Int3_3_3  = Int3_3_3+w(kk)*F3_3_3;
    end % for kk=1:Npp
    % Isub_1_1
     Int1_1_1 = ((PSI_1b-PSI_1a)/2)*Int1_1_1;
     Int2_1_1 = ((PSI_2b-PSI_2a)/2)*Int2_1_1;
     Int3_1_1 = ((PSI_3b-PSI_3a)/2)*Int3_1_1;
          
     Isub(1,1,m) = Int1_1_1+Int2_1_1+Int3_1_1;
     % Isub_1_2
     Int1_1_2 = ((PSI_1b-PSI_1a)/2)*Int1_1_2;
     Int2_1_2 = ((PSI_2b-PSI_2a)/2)*Int2_1_2;
     Int3_1_2 = ((PSI_3b-PSI_3a)/2)*Int3_1_2;
          
     Isub(1,2,m) = Int1_1_2+Int2_1_2+Int3_1_2;
     % Isub_1_3
     Int1_1_3 = ((PSI_1b-PSI_1a)/2)*Int1_1_3;
     Int2_1_3 = ((PSI_2b-PSI_2a)/2)*Int2_1_3;
     Int3_1_3 = ((PSI_3b-PSI_3a)/2)*Int3_1_3;
          
     Isub(1,3,m) = Int1_1_3+Int2_1_3+Int3_1_3;
     % Isub_2_1
     Int1_2_1 = ((PSI_1b-PSI_1a)/2)*Int1_2_1;
     Int2_2_1 = ((PSI_2b-PSI_2a)/2)*Int2_2_1;
     Int3_2_1 = ((PSI_3b-PSI_3a)/2)*Int3_2_1;
          
     Isub(2,1,m) = Int1_2_1+Int2_2_1+Int3_2_1;
     % Isub_2_2
     Int1_2_2 = ((PSI_1b-PSI_1a)/2)*Int1_2_2;
     Int2_2_2 = ((PSI_2b-PSI_2a)/2)*Int2_2_2;
     Int3_2_2 = ((PSI_3b-PSI_3a)/2)*Int3_2_2;
          
     Isub(2,2,m) = Int1_2_2+Int2_2_2+Int3_2_2;
     % Isub_2_3
     Int1_2_3 = ((PSI_1b-PSI_1a)/2)*Int1_2_3;
     Int2_2_3 = ((PSI_2b-PSI_2a)/2)*Int2_2_3;
     Int3_2_3 = ((PSI_3b-PSI_3a)/2)*Int3_2_3;
          
     Isub(2,3,m) = Int1_2_3+Int2_2_3+Int3_2_3;
     % Isub_3_1
     Int1_3_1 = ((PSI_1b-PSI_1a)/2)*Int1_3_1;
     Int2_3_1 = ((PSI_2b-PSI_2a)/2)*Int2_3_1;
     Int3_3_1 = ((PSI_3b-PSI_3a)/2)*Int3_3_1;
          
     Isub(3,1,m) = Int1_3_1+Int2_3_1+Int3_3_1;
     % Isub_3_2
     Int1_3_2 = ((PSI_1b-PSI_1a)/2)*Int1_3_2;
     Int2_3_2 = ((PSI_2b-PSI_2a)/2)*Int2_3_2;
     Int3_3_2 = ((PSI_3b-PSI_3a)/2)*Int3_3_2;
          
     Isub(3,2,m) = Int1_3_2+Int2_3_2+Int3_3_2;
     % Isub_3_3
     Int1_3_3 = ((PSI_1b-PSI_1a)/2)*Int1_3_3;
     Int2_3_3 = ((PSI_2b-PSI_2a)/2)*Int2_3_3;
     Int3_3_3 = ((PSI_3b-PSI_3a)/2)*Int3_3_3;
          
     Isub(3,3,m) = Int1_3_3+Int2_3_3+Int3_3_3;
     
end % for m = 1:3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I_simplex] = DIRECT_post(Isub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_1 = r2-r3;
l_2 = r3-r1;
l_3 = r1-r2;
%
L1 = norm(l_1);
L2 = norm(l_2);
L3 = norm(l_3);
%
L11 = dot(l_1,l_1);
L22 = dot(l_2,l_2);
L33 = dot(l_3,l_3);
%
L12 = dot(l_1,l_2);
L13 = dot(l_1,l_3);
L23 = dot(l_2,l_3);
%
I_vector(1) = L22 * I_simplex(9) - 2 * L23 * I_simplex(6)     + L33 * I_simplex(5);
I_vector(2) = L23 * I_simplex(7) - L12     * I_simplex(9)     - L33 * I_simplex(4) + L13 * I_simplex(6);
I_vector(3) = L12 * I_simplex(8) - L22     * I_simplex(7)     - L13 * I_simplex(5) + L23 * I_simplex(4);
I_vector(4) = L23 * I_simplex(3) - L33     * I_simplex(2)     - L12 * I_simplex(9) + L13 * I_simplex(8);
I_vector(5) = L33 * I_simplex(1) - 2 * L13 * I_simplex(3)     + L11 * I_simplex(9);
I_vector(6) = L13 * I_simplex(2) - L23     * I_simplex(1)     - L11 * I_simplex(8)  + L12 * I_simplex(7);
I_vector(7) = L12 * I_simplex(6) - L13     * I_simplex(5)     - L22 * I_simplex(3)  + L23 * I_simplex(2);
I_vector(8) = L13 * I_simplex(4) - L11     * I_simplex(6)     - L23 * I_simplex(1)  + L12 * I_simplex(3);
I_vector(9) = L11 * I_simplex(5) - 2 * L12 * I_simplex(4)     + L22 * I_simplex(1);
%
I_scalar    = I_simplex(10);
% Final output
I_DE(1) = (L1*L1) / 12 * (1i * ko * I_vector(1) + 4 / (1i *ko) * I_scalar);
I_DE(2) = (L1*L2) / 12 * (1i * ko * I_vector(2) + 4 / (1i *ko) * I_scalar);
I_DE(3) = (L1*L3) / 12 * (1i * ko * I_vector(3) + 4 / (1i *ko) * I_scalar);
I_DE(4) = (L2*L1) / 12 * (1i * ko * I_vector(4) + 4 / (1i *ko) * I_scalar);
I_DE(5) = (L2*L2) / 12 * (1i * ko * I_vector(5) + 4 / (1i *ko) * I_scalar);
I_DE(6) = (L2*L3) / 12 * (1i * ko * I_vector(6) + 4 / (1i *ko) * I_scalar);
I_DE(7) = (L3*L1) / 12 * (1i * ko * I_vector(7) + 4 / (1i *ko) * I_scalar);
I_DE(8) = (L3*L2) / 12 * (1i * ko * I_vector(8) + 4 / (1i *ko) * I_scalar);
I_DE(9) = (L3*L3) / 12 * (1i * ko * I_vector(9) + 4 / (1i *ko) * I_scalar);
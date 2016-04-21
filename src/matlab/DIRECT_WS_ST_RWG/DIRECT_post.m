function [I_simplex] = DIRECT_post(Isub)
%% Post-processing

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
% Isub = Integrals for the subtriangles

% OUTPUT DATA
%             I_simplex(1)  = I_z1_z1'
%             I_simplex(2)  = I_z1_z2'
%             I_simplex(3)  = I_z1_z3'
%             I_simplex(4)  = I_z2_z1'
%             I_simplex(5)  = I_z2_z2'
%             I_simplex(6)  = I_z2_z3'
%             I_simplex(7)  = I_z3_z1'
%             I_simplex(8)  = I_z3_z2'
%             I_simplex(9)  = I_z3_z3'
%             I_simplex(10) = I_const

Ieta_xi = zeros(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(1,1)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi(1,1) = Isub(1,1,1)+Isub(1,1,2)+Isub(1,1,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(1,2)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_1_2_a = Isub(1,2,1);

Ieta_xi_1_2_b = 1/2*Isub(1,1,2)-1/2*Isub(1,2,2)-sqrt(3)/2*Isub(1,3,2);

Ieta_xi_1_2_c = -1/2*Isub(1,1,3)-1/2*Isub(1,2,3)+sqrt(3)/2*Isub(1,3,3);

Ieta_xi(1,2)  = Ieta_xi_1_2_a+Ieta_xi_1_2_b+Ieta_xi_1_2_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(1,3)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_1_3_a = Isub(1,3,1);

Ieta_xi_1_3_b = sqrt(3)/2*Isub(1,1,2)+sqrt(3)/2*Isub(1,2,2)-1/2*Isub(1,3,2);

Ieta_xi_1_3_c = sqrt(3)/2*Isub(1,1,3)-sqrt(3)/2*Isub(1,2,3)-1/2*Isub(1,3,3);

Ieta_xi(1,3)  = Ieta_xi_1_3_a+Ieta_xi_1_3_b+Ieta_xi_1_3_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(2,1)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_2_1_a = Isub(2,1,1);

Ieta_xi_2_1_b = 1/2*Isub(1,1,2)-1/2*Isub(2,1,2)-sqrt(3)/2*Isub(3,1,2);

Ieta_xi_2_1_c = -1/2*Isub(1,1,3)-1/2*Isub(2,1,3)+sqrt(3)/2*Isub(3,1,3);

Ieta_xi(2,1)  = Ieta_xi_2_1_a+Ieta_xi_2_1_b+Ieta_xi_2_1_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(2,2)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_2_2_a = Isub(2,2,1);

Ieta_xi_2_2_b = 1/4*        Isub(1,1,2)   -1/4*         Isub(1,2,2)   -sqrt(3)/4*     Isub(1,3,2)...
               -1/4*        Isub(2,1,2)   +1/4*         Isub(2,2,2)   +sqrt(3)/4*     Isub(2,3,2)...
               -sqrt(3)/4*  Isub(3,1,2)   +sqrt(3)/4*   Isub(3,2,2)   +3/4*           Isub(3,3,2);

Ieta_xi_2_2_c = 1/4*       Isub(1,1,3)    +1/4*        Isub(1,2,3)    -sqrt(3)/4*      Isub(1,3,3)...
               +1/4*       Isub(2,1,3)    +1/4*        Isub(2,2,3)    -sqrt(3)/4*      Isub(2,3,3)...
               -sqrt(3)/4* Isub(3,1,3)    -sqrt(3)/4*  Isub(3,2,3)    +3/4*            Isub(3,3,3);
           
Ieta_xi(2,2)  = Ieta_xi_2_2_a+Ieta_xi_2_2_b+Ieta_xi_2_2_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(2,3)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_2_3_a = Isub(2,3,1);

Ieta_xi_2_3_b = sqrt(3)/4*  Isub(1,1,2)   +sqrt(3)/4*   Isub(1,2,2)   -1/4*           Isub(1,3,2)...
               -sqrt(3)/4*  Isub(2,1,2)   -sqrt(3)/4*   Isub(2,2,2)   +1/4*           Isub(2,3,2)...
               -3/4*        Isub(3,1,2)   -3/4*         Isub(3,2,2)   +sqrt(3)/4*     Isub(3,3,2);

Ieta_xi_2_3_c =-sqrt(3)/4* Isub(1,1,3)    +sqrt(3)/4*  Isub(1,2,3)    +1/4*            Isub(1,3,3)...
               -sqrt(3)/4* Isub(2,1,3)    +sqrt(3)/4*  Isub(2,2,3)    +1/4*            Isub(2,3,3)...
               +3/4*       Isub(3,1,3)    -3/4*        Isub(3,2,3)    -sqrt(3)/4*      Isub(3,3,3);
           
Ieta_xi(2,3)  = Ieta_xi_2_3_a+Ieta_xi_2_3_b+Ieta_xi_2_3_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(3,1)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_3_1_a = Isub(3,1,1);

Ieta_xi_3_1_b = sqrt(3)/2*Isub(1,1,2)+sqrt(3)/2*Isub(2,1,2)-1/2*Isub(3,1,2);

Ieta_xi_3_1_c = sqrt(3)/2*Isub(1,1,3)-sqrt(3)/2*Isub(2,1,3)-1/2*Isub(3,1,3);

Ieta_xi(3,1)  = Ieta_xi_3_1_a+Ieta_xi_3_1_b+Ieta_xi_3_1_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(3,2)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_3_2_a = Isub(3,2,1);

Ieta_xi_3_2_b = sqrt(3)/4*  Isub(1,1,2)   -sqrt(3)/4*   Isub(1,2,2)   -3/4*           Isub(1,3,2)...
               +sqrt(3)/4*  Isub(2,1,2)   -sqrt(3)/4*   Isub(2,2,2)   -3/4*           Isub(2,3,2)...
               -1/4*        Isub(3,1,2)   +1/4*         Isub(3,2,2)   +sqrt(3)/4*     Isub(3,3,2);

Ieta_xi_3_2_c =-sqrt(3)/4* Isub(1,1,3)    -sqrt(3)/4*  Isub(1,2,3)    +3/4*            Isub(1,3,3)...
               +sqrt(3)/4* Isub(2,1,3)    +sqrt(3)/4*  Isub(2,2,3)    -3/4*            Isub(2,3,3)...
               +1/4*       Isub(3,1,3)    +1/4*        Isub(3,2,3)    -sqrt(3)/4*      Isub(3,3,3);
           
Ieta_xi(3,2)  = Ieta_xi_3_2_a+Ieta_xi_3_2_b+Ieta_xi_3_2_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ieta_xi(3,3)                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ieta_xi_3_3_a = Isub(3,3,1);

Ieta_xi_3_3_b = 3/4*        Isub(1,1,2)   +3/4*         Isub(1,2,2)   -sqrt(3)/4*     Isub(1,3,2)...
               +3/4*        Isub(2,1,2)   +3/4*         Isub(2,2,2)   -sqrt(3)/4*     Isub(2,3,2)...
               -sqrt(3)/4*  Isub(3,1,2)   -sqrt(3)/4*   Isub(3,2,2)   +1/4*           Isub(3,3,2);

Ieta_xi_3_3_c = 3/4*       Isub(1,1,3)    -3/4*        Isub(1,2,3)    -sqrt(3)/4*      Isub(1,3,3)...
               -3/4*       Isub(2,1,3)    +3/4*        Isub(2,2,3)    +sqrt(3)/4*      Isub(2,3,3)...
               -sqrt(3)/4* Isub(3,1,3)    +sqrt(3)/4*  Isub(3,2,3)    +1/4*            Isub(3,3,3);
           
Ieta_xi(3,3)  = Ieta_xi_3_3_a+Ieta_xi_3_3_b+Ieta_xi_3_3_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising0                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_simplex(10) =  Ieta_xi(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(1,1)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_1_1 =  3*        Ieta_xi(1,1)   -3*         Ieta_xi(1,2)   -sqrt(3)*     Ieta_xi(1,3)...
            -3*        Ieta_xi(2,1)   +3*         Ieta_xi(2,2)   +sqrt(3)*     Ieta_xi(2,3)...
            -sqrt(3)*  Ieta_xi(3,1)   +sqrt(3)*   Ieta_xi(3,2)   +             Ieta_xi(3,3);

Ising_1_1    = 1/12*Ising_1_1;
           
I_simplex(1)   = Ising_1_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(1,2)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_1_2 =  3*        Ieta_xi(1,1)   +3*         Ieta_xi(1,2)   -sqrt(3)*     Ieta_xi(1,3)...
            -3*        Ieta_xi(2,1)   -3*         Ieta_xi(2,2)   +sqrt(3)*     Ieta_xi(2,3)...
            -sqrt(3)*  Ieta_xi(3,1)   -sqrt(3)*   Ieta_xi(3,2)   +             Ieta_xi(3,3);

Ising_1_2    = 1/12*Ising_1_2;
           
I_simplex(2)   = Ising_1_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(1,3)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_1_3 =  sqrt(3)*Ieta_xi(1,3) - sqrt(3)*Ieta_xi(2,3) - Ieta_xi(3,3);

Ising_1_3    = 1/6*Ising_1_3;
           
I_simplex(3)   = Ising_1_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(2,1)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_2_1 =  3*        Ieta_xi(1,1)   -3*         Ieta_xi(1,2)   -sqrt(3)*     Ieta_xi(1,3)...
            +3*        Ieta_xi(2,1)   -3*         Ieta_xi(2,2)   -sqrt(3)*     Ieta_xi(2,3)...
            -sqrt(3)*  Ieta_xi(3,1)   +sqrt(3)*   Ieta_xi(3,2)   +             Ieta_xi(3,3);

Ising_2_1    = 1/12*Ising_2_1;
           
I_simplex(4)   = Ising_2_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(2,2)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_2_2 =  3*        Ieta_xi(1,1)   +3*         Ieta_xi(1,2)   -sqrt(3)*     Ieta_xi(1,3)...
            +3*        Ieta_xi(2,1)   +3*         Ieta_xi(2,2)   -sqrt(3)*     Ieta_xi(2,3)...
            -sqrt(3)*  Ieta_xi(3,1)   -sqrt(3)*   Ieta_xi(3,2)   +             Ieta_xi(3,3);

Ising_2_2    = 1/12*Ising_2_2;
           
I_simplex(5)   = Ising_2_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(2,3)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_2_3 =  sqrt(3)*Ieta_xi(1,3) + sqrt(3)*Ieta_xi(2,3) - Ieta_xi(3,3);

Ising_2_3    = 1/6*Ising_2_3;
           
I_simplex(6)   = Ising_2_3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(3,1)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_3_1 =  sqrt(3)*Ieta_xi(3,1) - sqrt(3)*Ieta_xi(3,2) - Ieta_xi(3,3);

Ising_3_1    = 1/6*Ising_3_1;
           
I_simplex(7)   = Ising_3_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(3,2)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_3_2 =  sqrt(3)*Ieta_xi(3,1) + sqrt(3)*Ieta_xi(3,2) - Ieta_xi(3,3);

Ising_3_2    = 1/6*Ising_3_2;
           
I_simplex(8)   = Ising_3_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Ising(3,3)                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ising_3_3    = Ieta_xi(3,3);            

Ising_3_3    = 1/3*Ising_3_3;
           
I_simplex(9)   = Ising_3_3;
/*
 * Copyright (c) 2010 Athanasios Polimeridis
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU LGPL for more details.
 */


#include "demcem_ss_va_nxrwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"
#include <cmath>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void demcem_ss_va_nxrwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const complex<double> ko, const int N_theta_p, const int N_theta_q, const int N_psi, complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double rp1[3], rp2[3], rp3[3], rp21[3], rp31[3], Ap_cross[3];
	double Ap, Jp;
	double rq1[3], rq2[3], rq3[3], rq21[3], rq31[3], Aq_cross[3];
	double Aq, Jq;
	double r21[3], r31[3], r32[3], r41[3], r54[3], r51[3];
	//
	complex<double> Ising_g1_f1, Ising_g1_f2, Ising_g1_f3, Ising_g2_f1, Ising_g2_f2, Ising_g2_g3, Ising_g3_f1, Ising_g3_f2, Ising_g3_f3;
	//
	double theta_p_A, theta_p_B, THETA_p, J_theta_p, Lp;
	double theta_q_A, theta_q_B, THETA_q, J_theta_q, Lq;

	double b_v0[3], b_v1[3];

	double alpha_vo[3], alpha_v1[3], alpha_v2[3], alpha_v3[3];

	double a_v;
	double b_v;
	double c_v;

	double psi_1_A, psi_1_B, PSI_1, J_psi_1, psi_2_A, psi_2_B, PSI_2, J_psi_2;

	double B_1;
	double B_2;

	double L1, L2;

	complex<double> K_1[4], K_2[4];

	complex<double> Omega_g1_f1_1, Omega_g1_f2_1, Omega_g1_f3_1, Omega_g2_f1_1, Omega_g2_f2_1, Omega_g2_f3_1, Omega_g3_f1_1, Omega_g3_f2_1, Omega_g3_f3_1;
	complex<double> Omega_g1_f1_2, Omega_g1_f2_2, Omega_g1_f3_2, Omega_g2_f1_2, Omega_g2_f2_2, Omega_g2_f3_2, Omega_g3_f1_2, Omega_g3_f2_2, Omega_g3_f3_2;
	//
	double WPSI, WTHETA_p, WTHETA_q;
	//
	complex<double> I_theta_p_1_1, I_theta_p_1_2, I_theta_p_1_3, I_theta_p_2_1, I_theta_p_2_2, I_theta_p_2_3, I_theta_p_3_1, I_theta_p_3_2, I_theta_p_3_3; 
	complex<double> I_theta_q_1_1, I_theta_q_1_2, I_theta_q_1_3, I_theta_q_2_1, I_theta_q_2_2, I_theta_q_2_3, I_theta_q_3_1, I_theta_q_3_2, I_theta_q_3_3; 
	complex<double> I_psi_1_1, I_psi_1_2, I_psi_1_3, I_psi_2_1, I_psi_2_2, I_psi_2_3, I_psi_3_1, I_psi_3_2, I_psi_3_3; 

	// 2. Coefficients' parameters

	complex<double> coef_g1_f1[6];
	complex<double> coef_g1_f2[11];
	complex<double> coef_g1_f3[13];
	complex<double> coef_g2_f1[10];
	complex<double> coef_g2_f2[16];
	complex<double> coef_g2_f3[17];
	complex<double> coef_g3_f1[10];
	complex<double> coef_g3_f2[16];
	complex<double> coef_g3_f3[17];

	// 3. Quadrature parameters

	double* w_theta_p;                            // Pointers to arrays of integration points and weights for theta_p-interval
	double* z_theta_p;
	double* w_theta_q;                            // Pointers to arrays of integration points and weights for theta_q-interval
	double* z_theta_q;
	double* w_psi;                                // Pointers to arrays of integration points and weights for psi-interval
	double* z_psi;
	
	// ************************************************
	//			MAIN CODE
	// ************************************************
	// Evaluate Jacobians for equilateral parameter space
	for (int i = 0; i < 3; i++)
	{
			rp1[i] =  r1[i];
			rp2[i] =  r2[i];
			rp3[i] =  r3[i];
			rp21[i]   = rp2[i] - rp1[i];
	        rp31[i]   = rp3[i] - rp1[i];
	}
	//
    vector_cross(rp21,rp31,Ap_cross);
	Ap = 0.5 * sqrt(vector_dot(Ap_cross,Ap_cross));
	Jp = Ap / sqrt(3.0);
	//
	for (int i = 0; i < 3; i++)
	{
			rq1[i] =  r1[i];
			rq2[i] =  r4[i];
			rq3[i] =  r5[i];
			rq21[i]   = rq2[i] - rq1[i];
	        rq31[i]   = rq3[i] - rq1[i];
	}
	//
    vector_cross(rq21,rq31,Aq_cross);
	Aq = 0.5 * sqrt(vector_dot(Aq_cross,Aq_cross));
	Jq = Aq / sqrt(3.0);
	// Evaluate alpha, beta, gamma r21, r31, r32, r41, r42 parameters
	for (int i = 0; i < 3; i++)
	{
			alpha_vo[i] = (r2[i] - r1[i]) / 2;
			alpha_v1[i] = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			alpha_v2[i] = (r4[i] - r1[i]) / 2;
			alpha_v3[i] = (2 * r5[i] - r1[i] - r4[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r51[i]   = r5[i] - r1[i];
			r54[i]   = r5[i] - r4[i];
	}
	// Get the coefficients
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g1_f1, 11 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g1_f2, 12 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g1_f3, 13 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g2_f1, 21 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g2_f2, 22 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g2_f3, 23 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g3_f1, 31 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g3_f2, 32 );
	coefficients_ss_va_nxrwg ( r1,  r2,  r3,  r4, r5,  ko, Ap, coef_g3_f3, 33 );
	// Get the weights and abscissas for the 2 1-D quadratures
		// 1. Allocate space for the arrays
	 w_theta_p = new double [N_theta_p];
	 z_theta_p = new double [N_theta_p];

	 w_theta_q = new double [N_theta_q];
	 z_theta_q = new double [N_theta_q];

	 w_psi   = new double [N_psi];
	 z_psi   = new double [N_psi];
		// 2. Get the weights and abscissas
	 gl_quad ( N_theta_p, z_theta_p, w_theta_p );
	 gl_quad ( N_theta_q, z_theta_q, w_theta_q );
	 gl_quad ( N_psi, z_psi, w_psi );
     // Initialization of I_theta_p
	    I_theta_p_1_1 = 0.0;
		I_theta_p_1_2 = 0.0;
		I_theta_p_1_3 = 0.0;
		I_theta_p_2_1 = 0.0;
		I_theta_p_2_2 = 0.0;
		I_theta_p_2_3 = 0.0;
		I_theta_p_3_1 = 0.0;
		I_theta_p_3_2 = 0.0;
		I_theta_p_3_3 = 0.0;
	 //
     theta_p_A = 0.0;
     theta_p_B = M_PI / 3.0;
     //
     J_theta_p = (theta_p_B - theta_p_A) / 2.0;
	 for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
	 {

        THETA_p = ((theta_p_B - theta_p_A) / 2.0) * z_theta_p[n_theta_p] + (theta_p_B + theta_p_A) / 2.0;
        
        //
        Lp = (2.0 * sqrt(3.0) ) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p) );
        //
			for (int i = 0; i < 3; i++)
				{
					b_v0[i]   = alpha_vo[i] * cos(THETA_p) + alpha_v1[i] * sin(THETA_p);
				}
		// Initialization of I_theta_q
		I_theta_q_1_1 = 0.0;
		I_theta_q_1_2 = 0.0;
		I_theta_q_1_3 = 0.0;
		I_theta_q_2_1 = 0.0;
		I_theta_q_2_2 = 0.0;
		I_theta_q_2_3 = 0.0;
		I_theta_q_3_1 = 0.0;
		I_theta_q_3_2 = 0.0;
		I_theta_q_3_3 = 0.0;
		//
         theta_q_A = 0.0;
         theta_q_B = M_PI / 3.0;
         //
         J_theta_q = (theta_q_B - theta_q_A) / 2.0;
		 for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 {

			THETA_q = ((theta_q_B - theta_q_A) / 2.0) * z_theta_q[n_theta_q] + (theta_q_B + theta_q_A) / 2.0;
			
			//
			Lq = (2.0 * sqrt(3.0) ) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q) );
			//
			for (int i = 0; i < 3; i++)
				{
					b_v1[i]   = alpha_v2[i] * cos(THETA_q) + alpha_v3[i] * sin(THETA_q);
				}
			//
			a_v = vector_dot(b_v0,b_v0);
			b_v = -2.0 * vector_dot(b_v0,b_v1);
			c_v = vector_dot(b_v1,b_v1);
			// Initialization of I_psi
			I_psi_1_1 = 0.0;
			I_psi_1_2 = 0.0;
			I_psi_1_3 = 0.0;
			I_psi_2_1 = 0.0;
			I_psi_2_2 = 0.0;
			I_psi_2_3 = 0.0;
			I_psi_3_1 = 0.0;
			I_psi_3_2 = 0.0;
			I_psi_3_3 = 0.0;
			// psi_A =< PSI <= psi_B
             psi_1_A = 0.0;
             psi_1_B = atan(Lq / Lp);
             //
             J_psi_1 = (psi_1_B - psi_1_A) / 2.0;
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {

				 PSI_1 = ((psi_1_B-psi_1_A) / 2.0) * z_psi[n_psi] + (psi_1_B + psi_1_A) / 2.0;
				 
                 //
				 psi_2_A = atan(Lq/Lp);
				 psi_2_B = M_PI / 2.0;
				 //
				 PSI_2 = ((psi_2_B - psi_2_A) / 2.0) * z_psi[n_psi] + (psi_2_B + psi_2_A) / 2.0;
				 //
				 J_psi_2 = (psi_2_B - psi_2_A) / 2.0;
				 //
				 B_1 = sqrt(a_v * pow(cos(PSI_1),2.0) + b_v * cos(PSI_1) * sin(PSI_1) + c_v * pow(sin(PSI_1),2.0));
				 //
				 B_2 = sqrt(a_v * pow(cos(PSI_2),2.0) + b_v * cos(PSI_2) * sin(PSI_2) + c_v * pow(sin(PSI_2),2.0));
				 //
				 Lp = (2.0 * sqrt(3.0)) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p));
				 L1 = Lp / cos(PSI_1);

				 Lq = (2.0 * sqrt(3.0)) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q));
				 L2 = Lq / sin(PSI_2);
				 //

				 k_functions_ss_va_nxrwg(B_1, L1, K_1, ko);
				 k_functions_ss_va_nxrwg(B_2, L2, K_2, ko);
				 //
				 Omega_g1_f1_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g1_f1, ko, K_1, 11);
				 Omega_g1_f1_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g1_f1, ko, K_2, 11);
				 //
				 Omega_g1_f2_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g1_f2, ko, K_1, 12);
				 Omega_g1_f2_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g1_f2, ko, K_2, 12);
				 //
				 Omega_g1_f3_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g1_f3, ko, K_1, 13);
				 Omega_g1_f3_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g1_f3, ko, K_2, 13);
				 //
				 Omega_g2_f1_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g2_f1, ko, K_1, 21);
				 Omega_g2_f1_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g2_f1, ko, K_2, 21);
				 //
				 Omega_g2_f2_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g2_f2, ko, K_1, 22);
				 Omega_g2_f2_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g2_f2, ko, K_2, 22);
				 //
				 Omega_g2_f3_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g2_f3, ko, K_1, 23);
				 Omega_g2_f3_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g2_f3, ko, K_2, 23);
				 //
				 Omega_g3_f1_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g3_f1, ko, K_1, 31);
				 Omega_g3_f1_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g3_f1, ko, K_2, 31);
				 //
				 Omega_g3_f2_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g3_f2, ko, K_1, 32);
				 Omega_g3_f2_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g3_f2, ko, K_2, 32);
				 //
				 Omega_g3_f3_1 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_1, B_1, coef_g3_f3, ko, K_1, 33);
				 Omega_g3_f3_2 = omega_functions_ss_va_nxrwg(THETA_p, THETA_q, PSI_2, B_2, coef_g3_f3, ko, K_2, 33);
				 //
         		 WPSI = w_psi[n_psi];
				 //
				 I_psi_1_1   +=  WPSI * (J_psi_1 * Omega_g1_f1_1 + J_psi_2 * Omega_g1_f1_2);
				 I_psi_1_2   +=  WPSI * (J_psi_1 * Omega_g1_f2_1 + J_psi_2 * Omega_g1_f2_2);
				 I_psi_1_3   +=  WPSI * (J_psi_1 * Omega_g1_f3_1 + J_psi_2 * Omega_g1_f3_2);
				 I_psi_2_1   +=  WPSI * (J_psi_1 * Omega_g2_f1_1 + J_psi_2 * Omega_g2_f1_2);
				 I_psi_2_2   +=  WPSI * (J_psi_1 * Omega_g2_f2_1 + J_psi_2 * Omega_g2_f2_2);
				 I_psi_2_3   +=  WPSI * (J_psi_1 * Omega_g2_f3_1 + J_psi_2 * Omega_g2_f3_2);
				 I_psi_3_1   +=  WPSI * (J_psi_1 * Omega_g3_f1_1 + J_psi_2 * Omega_g3_f1_2);
				 I_psi_3_2   +=  WPSI * (J_psi_1 * Omega_g3_f2_1 + J_psi_2 * Omega_g3_f2_2);
				 I_psi_3_3   +=  WPSI * (J_psi_1 * Omega_g3_f3_1 + J_psi_2 * Omega_g3_f3_2);
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 WTHETA_q = w_theta_q[n_theta_q];
			 //
			 I_theta_q_1_1     +=  WTHETA_q * I_psi_1_1;
			 I_theta_q_1_2     +=  WTHETA_q * I_psi_1_2;
			 I_theta_q_1_3     +=  WTHETA_q * I_psi_1_3;
			 I_theta_q_2_1     +=  WTHETA_q * I_psi_2_1;
			 I_theta_q_2_2     +=  WTHETA_q * I_psi_2_2;
			 I_theta_q_2_3     +=  WTHETA_q * I_psi_2_3;
			 I_theta_q_3_1     +=  WTHETA_q * I_psi_3_1;
			 I_theta_q_3_2     +=  WTHETA_q * I_psi_3_2;
			 I_theta_q_3_3     +=  WTHETA_q * I_psi_3_3;
			 
		 } //end for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 //
		 I_theta_q_1_1     *=  J_theta_q ;
		 I_theta_q_1_2     *=  J_theta_q ;
		 I_theta_q_1_3     *=  J_theta_q ;
		 I_theta_q_2_1     *=  J_theta_q ;
		 I_theta_q_2_2     *=  J_theta_q ;
		 I_theta_q_2_3     *=  J_theta_q ;
		 I_theta_q_3_1     *=  J_theta_q ;
		 I_theta_q_3_2     *=  J_theta_q ;
		 I_theta_q_3_3     *=  J_theta_q ;
		 //
		 WTHETA_p = w_theta_p[n_theta_p];
		 //
		 I_theta_p_1_1    += WTHETA_p * I_theta_q_1_1;
		 I_theta_p_1_2    += WTHETA_p * I_theta_q_1_2;
		 I_theta_p_1_3    += WTHETA_p * I_theta_q_1_3;
		 I_theta_p_2_1    += WTHETA_p * I_theta_q_2_1;
		 I_theta_p_2_2    += WTHETA_p * I_theta_q_2_2;
		 I_theta_p_2_3    += WTHETA_p * I_theta_q_2_3;
		 I_theta_p_3_1    += WTHETA_p * I_theta_q_3_1;
		 I_theta_p_3_2    += WTHETA_p * I_theta_q_3_2;
		 I_theta_p_3_3    += WTHETA_p * I_theta_q_3_3;

	 } //end for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )

	 delete[] w_theta_p;
	 delete[] w_theta_q;
	 delete[] w_psi;
	 delete[] z_theta_p;
	 delete[] z_theta_q;
	 delete[] z_psi;
	 	 
	 //
	 I_theta_p_1_1     *=  J_theta_p ;
	 I_theta_p_1_2     *=  J_theta_p ;
	 I_theta_p_1_3     *=  J_theta_p ;
	 I_theta_p_2_1     *=  J_theta_p ;
	 I_theta_p_2_2     *=  J_theta_p ;
	 I_theta_p_2_3     *=  J_theta_p ;
	 I_theta_p_3_1     *=  J_theta_p ;
	 I_theta_p_3_2     *=  J_theta_p ;
	 I_theta_p_3_3     *=  J_theta_p ;
	 //
	 double PRECOMP = Jp * Jq / (4 * Ap * Aq);
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r54,r54)) * I_theta_p_1_1;
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r51,r51)) * I_theta_p_1_2;
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * I_theta_p_1_3;
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r54,r54)) * I_theta_p_2_1;
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r51,r51)) * I_theta_p_2_2;
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * I_theta_p_2_3;
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r54,r54)) * I_theta_p_3_1;
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r51,r51)) * I_theta_p_3_2;
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * I_theta_p_3_3;
}

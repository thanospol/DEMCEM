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


#include "demcem_ss_ea_nxrwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"
#include <cmath>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void demcem_ss_ea_nxrwg (const double r1[],const double r2[],const double r3[],const double r4[], const double ko, const int N_theta, const int N_psi, complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double rp1[3], rp2[3], rp3[3], rp21[3], rp31[3], Ap_cross[3];
	double Ap, Jp;
	double rq1[3], rq2[3], rq3[3], rq21[3], rq31[3], Aq_cross[3];
	double Aq, Jq;
	double alpha[3], beta[3], gamma[3], r21[3], r31[3], r32[3], r41[3], r42[3];
	//
	complex<double> Ising_g1_f1, Ising_g1_f2, Ising_g1_f3, Ising_g2_f1, Ising_g2_f2, Ising_g2_f3, Ising_g3_f1, Ising_g3_f2, Ising_g3_f3;
	//
	double THETA;
	double J_theta;
	double tPsiA;
	double tPsiB;
	double PsiA;
	double PsiB;
	double psi_A;
	double psi_B;
	//
	complex<double> N[12], Nm[12];
	//
	double PSI;
	double J_psi;
	double B;
	double Bm;
	//
	double b0;
	double b1;
	double b2;
	double B1;
	double B2;
	//
	double WPSI, WTHETA;
	//
	complex<double> I_theta_g1_f1, I_theta_g1_f2, I_theta_g1_f3, I_theta_g2_f1, I_theta_g2_f2, I_theta_g2_f3, I_theta_g3_f1, I_theta_g3_f2, I_theta_g3_f3;
	complex<double> I_psi_g1_f1, I_psi_g1_f2, I_psi_g1_f3, I_psi_g2_f1, I_psi_g2_f2, I_psi_g2_f3, I_psi_g3_f1, I_psi_g3_f2, I_psi_g3_f3;
	//
	complex<double> X_g1_f1, X_g1_f2, X_g1_f3, X_g2_f1, X_g2_f2, X_g2_f3, X_g3_f1, X_g3_f2, X_g3_f3;
	//
	complex<double> I_g1_f1[6], I_g1_f2[6], I_g1_f3[6], I_g2_f1[6], I_g2_f2[6], I_g2_f3[6], I_g3_f1[6], I_g3_f2[6], I_g3_f3[6]; 
	//
	complex<double> Ig1_f1, Ig1_f2, Ig1_f3, Ig2_f1, Ig2_f2, Ig2_f3, Ig3_f1, Ig3_f2, Ig3_f3; 

	// 2. Coefficients' parameters

	complex<double> coef_g1_f1[20], coefm_g1_f1[20];
	complex<double> coef_g1_f2[27], coefm_g1_f2[27];
	complex<double> coef_g1_f3[48], coefm_g1_f3[48];
	complex<double> coef_g2_f1[27], coefm_g2_f1[27];
	complex<double> coef_g2_f2[20], coefm_g2_f2[20];
	complex<double> coef_g2_f3[48], coefm_g2_f3[48];
	complex<double> coef_g3_f1[27], coefm_g3_f1[27];
	complex<double> coef_g3_f2[27], coefm_g3_f2[27];
	complex<double> coef_g3_f3[48], coefm_g3_f3[48];

	// 3. Quadrature parameters

	double* w_theta;                            // Pointers to arrays of integration points and weights for theta-interval
	double* z_theta;
	double* w_psi;                              // Pointers to arrays of integration points and weights for psi-interval
	double* z_psi;

	// 4.

	double theta_A;
	double theta_B;

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
			rq2[i] =  r2[i];
			rq3[i] =  r4[i];
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
			alpha[i] = (r2[i] - r1[i]) / 2;
			beta[i]  = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			gamma[i] = -(2 * r4[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r42[i]   = r4[i] - r2[i];
	}
	// Get the coefficients
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f1,  coefm_g1_f1, 11 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f2,  coefm_g1_f2, 12 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f3,  coefm_g1_f3, 13 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f1,  coefm_g2_f1, 21 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f2,  coefm_g2_f2, 22 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f3,  coefm_g2_f3, 23 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f1,  coefm_g3_f1, 31 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f2,  coefm_g3_f2, 32 );
	coefficients_ss_ea_nxrwg ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f3,  coefm_g3_f3, 33 );
	// Get the weights and abscissas for the 2 1-D quadratures
		// 1. Allocate space for the arrays
	 w_theta = new double [N_theta];
	 z_theta = new double [N_theta];

	 w_psi   = new double [N_psi];
	 z_psi   = new double [N_psi];
		// 2. Get the weights and abscissas
	 gl_quad ( N_theta, z_theta, w_theta );
	 gl_quad ( N_psi, z_psi, w_psi );
     // Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_g1_f1[im] = 0.0;
		 I_g1_f2[im] = 0.0;
		 I_g1_f3[im] = 0.0;
		 I_g2_f1[im] = 0.0;
		 I_g2_f2[im] = 0.0;
		 I_g2_f3[im] = 0.0;
		 I_g3_f1[im] = 0.0;
		 I_g3_f2[im] = 0.0;
		 I_g3_f3[im] = 0.0;
	 }
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_g1_f1 = 0.0;
		 I_theta_g1_f2 = 0.0;
		 I_theta_g1_f3 = 0.0;
		 I_theta_g2_f1 = 0.0;
		 I_theta_g2_f2 = 0.0;
		 I_theta_g2_f3 = 0.0;
		 I_theta_g3_f1 = 0.0;
		 I_theta_g3_f2 = 0.0;
		 I_theta_g3_f3 = 0.0;
		 //
		 theta_limits_ss_ea_nxrwg ( m, &theta_A, &theta_B );
		 //
         J_theta = (theta_B - theta_A) * 0.5;
		 for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 {
			 THETA   = ( (theta_B - theta_A) * 0.5) * z_theta[n_theta] + (theta_B + theta_A) * 0.5;
			 //
			 tPsiA   = sin(THETA) - sqrt(3.0) * cos(THETA);
			 tPsiB   = sin(THETA) + sqrt(3.0) * cos(THETA);
			 //
			 PsiA = atan(tPsiA);
			 PsiB = atan(tPsiB);
			 // Evaluate b0, b1, b2
			 b0 = vector_dot(beta,beta);
			 b1 = 2 * ( vector_dot(alpha,beta) * cos(THETA) + vector_dot(beta,gamma) * sin(THETA));
			 b2 = vector_dot(alpha,alpha)  * pow(cos(THETA),2) + 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     // Evaluate B1, B2 (B_plus, B_minus)
			 B1 = 2 * ( vector_dot(beta,gamma) * sin(THETA) - vector_dot(alpha,beta) * cos(THETA));
			 B2 = vector_dot(alpha,alpha) * pow(cos(THETA),2) - 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     //
			 psi_limits_ss_ea_nxrwg ( m, PsiA, PsiB, &psi_A, &psi_B );
			 //
             J_psi = (psi_B - psi_A) * 0.5;
             //
			 I_psi_g1_f1 = 0.0;
			 I_psi_g1_f2 = 0.0;
			 I_psi_g1_f3 = 0.0;
			 I_psi_g2_f1 = 0.0;
			 I_psi_g2_f2 = 0.0;
			 I_psi_g2_f3 = 0.0;
			 I_psi_g3_f1 = 0.0;
			 I_psi_g3_f2 = 0.0;
			 I_psi_g3_f3 = 0.0;
			 //
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 PSI   = ( (psi_B - psi_A) * 0.5) * z_psi[n_psi] + (psi_B + psi_A) * 0.5;				 
				 //				 
				 B     = sqrt(b0 * pow(sin(PSI),2) + b1 * cos(PSI) * sin(PSI) + b2 * pow(cos(PSI),2) );
				 Bm    = sqrt(b0 * pow(sin(PSI),2) + B1 * cos(PSI) * sin(PSI) + B2 * pow(cos(PSI),2) );
				 // Pre-processing: Evaluate N[12], Nm[12]
				 x_function_pre_ss_ea_nxrwg ( THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, N, Nm, ko);
				 // Evaluate the integrands
				 X_g1_f1 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g1_f1, coefm_g1_f1, N, Nm, 11);
				 X_g1_f2 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g1_f2, coefm_g1_f2, N, Nm, 12);
				 X_g1_f3 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g1_f3, coefm_g1_f3, N, Nm, 13);
				 X_g2_f1 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g2_f1, coefm_g2_f1, N, Nm, 21);
				 X_g2_f2 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g2_f2, coefm_g2_f2, N, Nm, 22);
				 X_g2_f3 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g2_f3, coefm_g2_f3, N, Nm, 23);
				 X_g3_f1 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g3_f1, coefm_g3_f1, N, Nm, 31);
				 X_g3_f2 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g3_f2, coefm_g3_f2, N, Nm, 32);
				 X_g3_f3 = x_functions_ss_ea_nxrwg (THETA, PSI, B, Bm, coef_g3_f3, coefm_g3_f3, N, Nm, 33);
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_g1_f1 +=  WPSI * X_g1_f1;
				 I_psi_g1_f2 +=  WPSI * X_g1_f2;
				 I_psi_g1_f3 +=  WPSI * X_g1_f3;
				 I_psi_g2_f1 +=  WPSI * X_g2_f1;
				 I_psi_g2_f2 +=  WPSI * X_g2_f2;
				 I_psi_g2_f3 +=  WPSI * X_g2_f3;
				 I_psi_g3_f1 +=  WPSI * X_g3_f1;
				 I_psi_g3_f2 +=  WPSI * X_g3_f2;
				 I_psi_g3_f3 +=  WPSI * X_g3_f3;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_g1_f1 *=  J_psi;
			 I_psi_g1_f2 *=  J_psi;
			 I_psi_g1_f3 *=  J_psi;
			 I_psi_g2_f1 *=  J_psi;
			 I_psi_g2_f2 *=  J_psi;
			 I_psi_g2_f3 *=  J_psi;
			 I_psi_g3_f1 *=  J_psi;
			 I_psi_g3_f2 *=  J_psi;
			 I_psi_g3_f3 *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_g1_f1 +=  WTHETA * I_psi_g1_f1;
			 I_theta_g1_f2 +=  WTHETA * I_psi_g1_f2;
			 I_theta_g1_f3 +=  WTHETA * I_psi_g1_f3;
			 I_theta_g2_f1 +=  WTHETA * I_psi_g2_f1;
			 I_theta_g2_f2 +=  WTHETA * I_psi_g2_f2;
			 I_theta_g2_f3 +=  WTHETA * I_psi_g2_f3;
			 I_theta_g3_f1 +=  WTHETA * I_psi_g3_f1;
			 I_theta_g3_f2 +=  WTHETA * I_psi_g3_f2;
			 I_theta_g3_f3 +=  WTHETA * I_psi_g3_f3;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_g1_f1[m-1] = J_theta * I_theta_g1_f1;
		 I_g1_f2[m-1] = J_theta * I_theta_g1_f2;
		 I_g1_f3[m-1] = J_theta * I_theta_g1_f3;
		 I_g2_f1[m-1] = J_theta * I_theta_g2_f1;
		 I_g2_f2[m-1] = J_theta * I_theta_g2_f2;
		 I_g2_f3[m-1] = J_theta * I_theta_g2_f3;
		 I_g3_f1[m-1] = J_theta * I_theta_g3_f1;
		 I_g3_f2[m-1] = J_theta * I_theta_g3_f2;
		 I_g3_f3[m-1] = J_theta * I_theta_g3_f3;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 Ig1_f1 = I_g1_f1[0] + I_g1_f1[1] + I_g1_f1[2] + I_g1_f1[3] + I_g1_f1[4] + I_g1_f1[5];
	 Ig1_f2 = I_g1_f2[0] + I_g1_f2[1] + I_g1_f2[2] + I_g1_f2[3] + I_g1_f2[4] + I_g1_f2[5];
	 Ig1_f3 = I_g1_f3[0] + I_g1_f3[1] + I_g1_f3[2] + I_g1_f3[3] + I_g1_f3[4] + I_g1_f3[5];
	 Ig2_f1 = I_g2_f1[0] + I_g2_f1[1] + I_g2_f1[2] + I_g2_f1[3] + I_g2_f1[4] + I_g2_f1[5];
	 Ig2_f2 = I_g2_f2[0] + I_g2_f2[1] + I_g2_f2[2] + I_g2_f2[3] + I_g2_f2[4] + I_g2_f2[5];
	 Ig2_f3 = I_g2_f3[0] + I_g2_f3[1] + I_g2_f3[2] + I_g2_f3[3] + I_g2_f3[4] + I_g2_f3[5];
	 Ig3_f1 = I_g3_f1[0] + I_g3_f1[1] + I_g3_f1[2] + I_g3_f1[3] + I_g3_f1[4] + I_g3_f1[5];
	 Ig3_f2 = I_g3_f2[0] + I_g3_f2[1] + I_g3_f2[2] + I_g3_f2[3] + I_g3_f2[4] + I_g3_f2[5];
	 Ig3_f3 = I_g3_f3[0] + I_g3_f3[1] + I_g3_f3[2] + I_g3_f3[3] + I_g3_f3[4] + I_g3_f3[5];
	 //
	 double PRECOMP = Jp * Jq / (4 * Ap * Aq);
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * Ig1_f1;
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r42,r42)) * Ig1_f2;
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r21,r21)) * Ig1_f3;
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * Ig2_f1;
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r42,r42)) * Ig2_f2;
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r21,r21)) * Ig2_f3;
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * Ig3_f1;
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r42,r42)) * Ig3_f2;
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r21,r21)) * Ig3_f3;
}

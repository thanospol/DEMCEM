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

#include <cmath>
#include "demcem_ws_ea_rwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"


using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void demcem_ws_ea_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double ko, const int N_theta, const int N_psi, complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double alpha[3], beta[3], gamma[3], r21[3], r31[3], r32[3], r41[3], r42[3];
	//
	complex<double> Ising_const, Ising_f1_f1, Ising_f1_f2, Ising_f1_f3, Ising_f2_f1, Ising_f2_f2, Ising_f2_f3, Ising_f3_f1, Ising_f3_f2, Ising_f3_f3;
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
	complex<double> I_theta_const, I_theta_f1_f1, I_theta_f1_f2, I_theta_f1_f3, I_theta_f2_f1, I_theta_f2_f2, I_theta_f2_f3, I_theta_f3_f1, I_theta_f3_f2, I_theta_f3_f3;
	complex<double> I_psi_const, I_psi_f1_f1, I_psi_f1_f2, I_psi_f1_f3, I_psi_f2_f1, I_psi_f2_f2, I_psi_f2_f3, I_psi_f3_f1, I_psi_f3_f2, I_psi_f3_f3;
	//
	complex<double> X_const, X_f1_f1, X_f1_f2, X_f1_f3, X_f2_f1, X_f2_f2, X_f2_f3, X_f3_f1, X_f3_f2, X_f3_f3;
	//
	complex<double> I_const[6], I_f1_f1[6], I_f1_f2[6], I_f1_f3[6], I_f2_f1[6], I_f2_f2[6], I_f2_f3[6], I_f3_f1[6], I_f3_f2[6], I_f3_f3[6]; 
	//
	complex<double> Iconst, If1_f1, If1_f2, If1_f3, If2_f1, If2_f2, If2_f3, If3_f1, If3_f2, If3_f3; 

	// 2. Coefficients' parameters

	complex<double> coef_const[3],  coefm_const[3];
	complex<double> coef_f1_f1[40], coefm_f1_f1[40];
	complex<double> coef_f1_f2[43], coefm_f1_f2[43];
	complex<double> coef_f1_f3[43], coefm_f1_f3[43];
	complex<double> coef_f2_f1[43], coefm_f2_f1[43];
	complex<double> coef_f2_f2[40], coefm_f2_f2[40];
	complex<double> coef_f2_f3[43], coefm_f2_f3[43];
	complex<double> coef_f3_f1[43], coefm_f3_f1[43];
	complex<double> coef_f3_f2[43], coefm_f3_f2[43];
	complex<double> coef_f3_f3[43], coefm_f3_f3[43];

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
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_const,  coefm_const, 0 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f1_f1,  coefm_f1_f1, 11 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f1_f2,  coefm_f1_f2, 12 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f1_f3,  coefm_f1_f3, 13 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f2_f1,  coefm_f2_f1, 21 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f2_f2,  coefm_f2_f2, 22 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f2_f3,  coefm_f2_f3, 23 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f3_f1,  coefm_f3_f1, 31 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f3_f2,  coefm_f3_f2, 32 );
	coefficients_ws_ea_rwg ( r1,  r2,  r3,  r4,  ko,  coef_f3_f3,  coefm_f3_f3, 33 );
	// Get the weights and abscissas for the 2 1-D quadratures
		// 1. Allocate space for the arrays
	 w_theta = new double [N_theta];
	 z_theta = new double [N_theta];

	 w_psi   = new double [N_psi];
	 z_psi   = new double [N_psi];
		// 2. Get the weights and abscissas
	 gl_quad ( N_theta, z_theta, w_theta );
	 gl_quad ( N_psi, z_psi, w_psi );
     /*// Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_const[im] = 0.0;
		 I_f1_f1[im] = 0.0;
		 I_f1_f2[im] = 0.0;
		 I_f1_f3[im] = 0.0;
		 I_f2_f1[im] = 0.0;
		 I_f2_f2[im] = 0.0;
		 I_f2_f3[im] = 0.0;
		 I_f3_f1[im] = 0.0;
		 I_f3_f2[im] = 0.0;
		 I_f3_f3[im] = 0.0;
	 }*/
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_const = 0.0;
		 I_theta_f1_f1 = 0.0;
		 I_theta_f1_f2 = 0.0;
		 I_theta_f1_f3 = 0.0;
		 I_theta_f2_f1 = 0.0;
		 I_theta_f2_f2 = 0.0;
		 I_theta_f2_f3 = 0.0;
		 I_theta_f3_f1 = 0.0;
		 I_theta_f3_f2 = 0.0;
		 I_theta_f3_f3 = 0.0;
		 //
		 theta_limits_ws_ea_rwg ( m, &theta_A, &theta_B );
         J_theta = (theta_B - theta_A) * 0.5;
		 //
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
			 I_psi_const = 0.0;
			 I_psi_f1_f1 = 0.0;
			 I_psi_f1_f2 = 0.0;
			 I_psi_f1_f3 = 0.0;
			 I_psi_f2_f1 = 0.0;
			 I_psi_f2_f2 = 0.0;
			 I_psi_f2_f3 = 0.0;
			 I_psi_f3_f1 = 0.0;
			 I_psi_f3_f2 = 0.0;
			 I_psi_f3_f3 = 0.0;
			 //
             //
             psi_limits_ws_ea_rwg ( m, PsiA, PsiB, &psi_A, &psi_B );
             J_psi = (psi_B - psi_A) * 0.5;
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 PSI   = ( (psi_B - psi_A) * 0.5) * z_psi[n_psi] + (psi_B + psi_A) * 0.5;
				 //				 
				 B     = sqrt(b0 * pow(sin(PSI),2) + b1 * cos(PSI) * sin(PSI) + b2 * pow(cos(PSI),2) );
				 Bm    = sqrt(b0 * pow(sin(PSI),2) + B1 * cos(PSI) * sin(PSI) + B2 * pow(cos(PSI),2) );
				 // Pre-processing: Evaluate N[12], Nm[12]
				 x_functions_pre_ws_ea_rwg ( THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, N, Nm, ko);
				 // Evaluate the integrands
				 X_const = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_const, coefm_const, N, Nm, 0);
				 X_f1_f1 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f1_f1, coefm_f1_f1, N, Nm, 11);
				 X_f1_f2 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f1_f2, coefm_f1_f2, N, Nm, 12);
				 X_f1_f3 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f1_f3, coefm_f1_f3, N, Nm, 13);
				 X_f2_f1 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f2_f1, coefm_f2_f1, N, Nm, 21);
				 X_f2_f2 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f2_f2, coefm_f2_f2, N, Nm, 22);
				 X_f2_f3 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f2_f3, coefm_f2_f3, N, Nm, 23);
				 X_f3_f1 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f3_f1, coefm_f3_f1, N, Nm, 31);
				 X_f3_f2 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f3_f2, coefm_f3_f2, N, Nm, 32);
				 X_f3_f3 = x_functions_ws_ea_rwg (THETA, PSI, B, Bm, coef_f3_f3, coefm_f3_f3, N, Nm, 33);
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_const +=  WPSI * X_const;
				 I_psi_f1_f1 +=  WPSI * X_f1_f1;
				 I_psi_f1_f2 +=  WPSI * X_f1_f2;
				 I_psi_f1_f3 +=  WPSI * X_f1_f3;
				 I_psi_f2_f1 +=  WPSI * X_f2_f1;
				 I_psi_f2_f2 +=  WPSI * X_f2_f2;
				 I_psi_f2_f3 +=  WPSI * X_f2_f3;
				 I_psi_f3_f1 +=  WPSI * X_f3_f1;
				 I_psi_f3_f2 +=  WPSI * X_f3_f2;
				 I_psi_f3_f3 +=  WPSI * X_f3_f3;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_const *=  J_psi;
			 I_psi_f1_f1 *=  J_psi;
			 I_psi_f1_f2 *=  J_psi;
			 I_psi_f1_f3 *=  J_psi;
			 I_psi_f2_f1 *=  J_psi;
			 I_psi_f2_f2 *=  J_psi;
			 I_psi_f2_f3 *=  J_psi;
			 I_psi_f3_f1 *=  J_psi;
			 I_psi_f3_f2 *=  J_psi;
			 I_psi_f3_f3 *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_const +=  WTHETA * I_psi_const;
			 I_theta_f1_f1 +=  WTHETA * I_psi_f1_f1;
			 I_theta_f1_f2 +=  WTHETA * I_psi_f1_f2;
			 I_theta_f1_f3 +=  WTHETA * I_psi_f1_f3;
			 I_theta_f2_f1 +=  WTHETA * I_psi_f2_f1;
			 I_theta_f2_f2 +=  WTHETA * I_psi_f2_f2;
			 I_theta_f2_f3 +=  WTHETA * I_psi_f2_f3;
			 I_theta_f3_f1 +=  WTHETA * I_psi_f3_f1;
			 I_theta_f3_f2 +=  WTHETA * I_psi_f3_f2;
			 I_theta_f3_f3 +=  WTHETA * I_psi_f3_f3;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_const[m-1] = J_theta * I_theta_const;
		 I_f1_f1[m-1] = J_theta * I_theta_f1_f1;
		 I_f1_f2[m-1] = J_theta * I_theta_f1_f2;
		 I_f1_f3[m-1] = J_theta * I_theta_f1_f3;
		 I_f2_f1[m-1] = J_theta * I_theta_f2_f1;
		 I_f2_f2[m-1] = J_theta * I_theta_f2_f2;
		 I_f2_f3[m-1] = J_theta * I_theta_f2_f3;
		 I_f3_f1[m-1] = J_theta * I_theta_f3_f1;
		 I_f3_f2[m-1] = J_theta * I_theta_f3_f2;
		 I_f3_f3[m-1] = J_theta * I_theta_f3_f3;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 Iconst = I_const[0] + I_const[1] + I_const[2] + I_const[3] + I_const[4] + I_const[5];
	 If1_f1 = I_f1_f1[0] + I_f1_f1[1] + I_f1_f1[2] + I_f1_f1[3] + I_f1_f1[4] + I_f1_f1[5];
	 If1_f2 = I_f1_f2[0] + I_f1_f2[1] + I_f1_f2[2] + I_f1_f2[3] + I_f1_f2[4] + I_f1_f2[5];
	 If1_f3 = I_f1_f3[0] + I_f1_f3[1] + I_f1_f3[2] + I_f1_f3[3] + I_f1_f3[4] + I_f1_f3[5];
	 If2_f1 = I_f2_f1[0] + I_f2_f1[1] + I_f2_f1[2] + I_f2_f1[3] + I_f2_f1[4] + I_f2_f1[5];
	 If2_f2 = I_f2_f2[0] + I_f2_f2[1] + I_f2_f2[2] + I_f2_f2[3] + I_f2_f2[4] + I_f2_f2[5];
	 If2_f3 = I_f2_f3[0] + I_f2_f3[1] + I_f2_f3[2] + I_f2_f3[3] + I_f2_f3[4] + I_f2_f3[5];
	 If3_f1 = I_f3_f1[0] + I_f3_f1[1] + I_f3_f1[2] + I_f3_f1[3] + I_f3_f1[4] + I_f3_f1[5];
	 If3_f2 = I_f3_f2[0] + I_f3_f2[1] + I_f3_f2[2] + I_f3_f2[3] + I_f3_f2[4] + I_f3_f2[5];
	 If3_f3 = I_f3_f3[0] + I_f3_f3[1] + I_f3_f3[2] + I_f3_f3[3] + I_f3_f3[4] + I_f3_f3[5];
	 //
	 double PRECOMP = 1.0 / 12.0;
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If1_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If1_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If1_f3 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If2_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If2_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If2_f3 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * (Iunit * ko * If3_f1 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r42,r42)) * (Iunit * ko * If3_f2 + 4.0 / (Iunit * ko) * Iconst);
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r21,r21)) * (Iunit * ko * If3_f3 + 4.0 / (Iunit * ko) * Iconst);
}

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
#include "demcem_ws_st_rwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT_ST_RWG
// ***********************************************************************

void demcem_ws_st_rwg (const double r1[],const double r2[],const double r3[], const double ko, const int Np_1D, complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double l_1[3], l_2[3], l_3[3], rp12[3], rp123[3];
	double Asing0, Asing1, Asing2, Asing3, Asing4, Asing5, Asing6, Asing7, Asing8;
	double Asing[3][3];
	double acc, acs, ass;

	complex<double> Int1, Int2, Int3;

	double  PSI_1k, a_PSI_1k;
	double  PSI_2k, a_PSI_2k;
	double  PSI_3k, a_PSI_3k;

	complex<double> PHI_a, PHI_b, PHI_c, PHI_d, PHI_e, PHI_f, PHI_g, PHI_h, F1, F2, F3;

	complex<double> Isub[3][3][3];

	complex<double> I_simplex[10], I_vector[9], I_scalar;
	//

	// 2. Quadrature parameters

	double* w;                            // Pointers to arrays of integration points and weights 
	double* z;

	// ************************************************
	//			MAIN CODE
	// ************************************************
	for (int i = 0; i < 3; i++)
	{
			rp12[i]    = r1[i] - r2[i];
			rp123[i]   = r1[i] + r2[i] - 2.0 * r3[i];
			//
			l_1[i]     = r2[i] - r3[i];
			l_2[i]     = r3[i] - r1[i];
			l_3[i]     = r1[i] - r2[i];
	}
	// Evaluate A parameters
	Asing0 = vector_dot(rp12,rp12) / 4.0;
	Asing1 = vector_dot(rp12,rp123)  / (2.0 * sqrt(3.0));
	Asing2 = vector_dot(rp123,rp123) / 12.0; 
	//
	Asing3 = (1.0/4.0)        *Asing0   -(sqrt(3.0)/4.0)  *Asing1      +(3.0/4.0)        *Asing2;
	Asing4 = (sqrt(3.0)/2.0)  *Asing0   -(1.0/2.0)        *Asing1      -(sqrt(3.0)/2.0)  *Asing2;
	Asing5 = (3.0/4.0)        *Asing0   +(sqrt(3.0)/4.0)  *Asing1      +(1.0/4.0)        *Asing2;
	//
	Asing6 = (1.0/4.0)        *Asing0   +(sqrt(3.0)/4.0)  *Asing1       +(3.0/4.0)       *Asing2;
	Asing7 = -(sqrt(3.0)/2.0) *Asing0   -(1.0/2.0)        *Asing1       +(sqrt(3.0)/2.0) *Asing2;
	Asing8 = (3.0/4.0)        *Asing0   -(sqrt(3.0)/4.0)  *Asing1       +(1.0/4.0)       *Asing2;
	//
	Asing[0][0] = Asing0;
	Asing[1][0] = Asing1;
	Asing[2][0] = Asing2;
	Asing[0][1] = Asing3;
	Asing[1][1] = Asing4;
	Asing[2][1] = Asing5;
	Asing[0][2] = Asing6;
	Asing[1][2] = Asing7;
	Asing[2][2] = Asing8;
	// Get the weights and abscissas for the 1-D quadrature
	// 1. Allocate space for the arrays
	 w = new double [Np_1D];
	 z = new double [Np_1D];
	// 2. Get the weights and abscissas
	 gl_quad ( Np_1D, z, w );
    
    // assign values to the constant angles
    double PSI_1a = 0.0;
    double PSI_1b = M_PI / 3.0;
    double PSI_2a = M_PI / 3.0;
    double PSI_2b = 2.0 * M_PI / 3.0;
    double PSI_3a = 2.0 * M_PI / 3.0;
    double PSI_3b = M_PI;

	 // for loops

	 for ( int column = 0 ; column <  3 ; column++ )
	 {
		 for ( int row = 0 ; row <  3 ; row++ )
		 {
			 for ( int m = 0 ; m <  3 ; m++ )
			 {
				 // Define the coefs of the appropriate subtriangle
				 acc = Asing[0][m];
				 acs = Asing[1][m];
				 ass = Asing[2][m];
				 // Gauss quadrature
				 Int1 = 0.0;
				 Int2 = 0.0;
				 Int3 = 0.0;
				 //
				 for ( int kk = 0 ; kk <  Np_1D ; kk++ )
				 {
					 // Int1,  0 =< PSI <= pi/3
					 
					 PSI_1k = ((PSI_1b - PSI_1a) * z[kk] + (PSI_1b + PSI_1a)) / 2.0;
					 a_PSI_1k = sqrt(acc * pow(cos(PSI_1k),2.0) - acs * cos(PSI_1k) * sin(PSI_1k) + ass * pow(sin(PSI_1k),2.0));

					 PHI_a = phi_functions_ws_st_rwg(row, column, 1, PSI_1k, a_PSI_1k, ko);
					 PHI_e = phi_functions_ws_st_rwg(row, column, 5, PSI_1k, a_PSI_1k, ko);
					 PHI_f = phi_functions_ws_st_rwg(row, column, 6, PSI_1k, a_PSI_1k, ko);

					 F1 = PHI_a + PHI_e + PHI_f;       
					 Int1  = Int1 + w[kk] * F1;
					 // Int2,  pi/3 =< PSI <= 2pi/3

					 PSI_2k = ((PSI_2b - PSI_2a) * z[kk] + (PSI_2b + PSI_2a)) / 2.0;
					 a_PSI_2k = sqrt(acc * pow(cos(PSI_2k),2.0) - acs * cos(PSI_2k) * sin(PSI_2k) + ass * pow(sin(PSI_2k),2.0));

					 PHI_b = phi_functions_ws_st_rwg(row, column, 2, PSI_2k, a_PSI_2k, ko);
					 PHI_g = phi_functions_ws_st_rwg(row, column, 7, PSI_2k, a_PSI_2k, ko);

					 F2 = PHI_b + PHI_g ;       
					 Int2  = Int2 + w[kk] * F2;
					 // Int3,  2pi/3 =< PSI <= pi
					 PSI_3k = ((PSI_3b - PSI_3a) * z[kk] + (PSI_3b + PSI_3a)) / 2.0;
					 a_PSI_3k = sqrt(acc * pow(cos(PSI_3k),2.0) - acs * cos(PSI_3k) * sin(PSI_3k) + ass * pow(sin(PSI_3k),2.0));

					 PHI_c = phi_functions_ws_st_rwg(row, column, 3, PSI_3k, a_PSI_3k, ko);
					 PHI_d = phi_functions_ws_st_rwg(row, column, 4, PSI_3k, a_PSI_3k, ko);
					 PHI_h = phi_functions_ws_st_rwg(row, column, 8, PSI_3k, a_PSI_3k, ko);

					 F3 = PHI_c + PHI_d + PHI_h ;       
					 Int3  = Int3 + w[kk] * F3;

				 }// for ( int kk = 0 ; kk <  Np_1D ; kk++ )

				 Int1 = ((PSI_1b - PSI_1a) / 2.0) * Int1;
				 Int2 = ((PSI_2b - PSI_2a) / 2.0) * Int2;
				 Int3 = ((PSI_3b - PSI_3a) / 2.0) * Int3;
				 //

				 Isub[row][column][m] = Int1 + Int2 + Int3;

			 }// for ( int m = 0 ; m <  3 ; m++ )

		 }// for ( int row = 0 ; row <  3 ; row++ )
	 }// for ( int column = 0 ; column <  3 ; column++ )

	// Post-processing

	post_ws_st_rwg ( Isub, I_simplex );

	// 
	double L11  = vector_dot(l_1,l_1);
	double L22  = vector_dot(l_2,l_2);
	double L33  = vector_dot(l_3,l_3);
	//
	double L12  = vector_dot(l_1,l_2);
	double L13  = vector_dot(l_1,l_3);
	double L23  = vector_dot(l_2,l_3);
	//
	double L1   = sqrt(L11);
	double L2   = sqrt(L22);
	double L3   = sqrt(L33);
	//
	I_vector[0] = L22 * I_simplex[8] - 2 * L23 * I_simplex[5]     + L33 * I_simplex[4];
	I_vector[1] = L23 * I_simplex[6] - L12     * I_simplex[8]     - L33 * I_simplex[3] + L13 * I_simplex[5];
	I_vector[2] = L12 * I_simplex[7] - L22     * I_simplex[6]     - L13 * I_simplex[4] + L23 * I_simplex[3];
	I_vector[3] = L23 * I_simplex[2] - L33     * I_simplex[1]     - L12 * I_simplex[8] + L13 * I_simplex[7];
	I_vector[4] = L33 * I_simplex[0] - 2 * L13 * I_simplex[2]     + L11 * I_simplex[8];
	I_vector[5] = L13 * I_simplex[1] - L23     * I_simplex[0]     - L11 * I_simplex[7]  + L12 * I_simplex[6];
	I_vector[6] = L12 * I_simplex[5] - L13     * I_simplex[4]     - L22 * I_simplex[2]  + L23 * I_simplex[1];
	I_vector[7] = L13 * I_simplex[3] - L11     * I_simplex[5]     - L23 * I_simplex[0]  + L12 * I_simplex[2];
	I_vector[8] = L11 * I_simplex[4] - 2 * L12 * I_simplex[3]     + L22 * I_simplex[0];
	//
	I_scalar    = I_simplex[9];
	// Final output
	I_DE[0] = (L1*L1) / 12.0 * (Iunit * ko * I_vector[0] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[1] = (L1*L2) / 12.0 * (Iunit * ko * I_vector[1] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[2] = (L1*L3) / 12.0 * (Iunit * ko * I_vector[2] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[3] = (L2*L1) / 12.0 * (Iunit * ko * I_vector[3] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[4] = (L2*L2) / 12.0 * (Iunit * ko * I_vector[4] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[5] = (L2*L3) / 12.0 * (Iunit * ko * I_vector[5] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[6] = (L3*L1) / 12.0 * (Iunit * ko * I_vector[6] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[7] = (L3*L2) / 12.0 * (Iunit * ko * I_vector[7] + 4.0 / (Iunit *ko) * I_scalar);
	I_DE[8] = (L3*L3) / 12.0 * (Iunit * ko * I_vector[8] + 4.0 / (Iunit *ko) * I_scalar);
}
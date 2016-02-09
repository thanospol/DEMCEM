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

#include "demcem_ws_st_rwg.h"
#include <cmath>

using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT_post
// ***********************************************************************

void post_ws_st_rwg ( complex<double> Isub[3][3][3], complex<double> I_simplex[] )
{
	complex<double> Ieta_xi[3][3];

	// Ieta_xi(1,1)
	Ieta_xi[0][0] = Isub[0][0][0] + Isub[0][0][1] + Isub[0][0][2];
	// Ieta_xi(1,2)
	complex<double> Ieta_xi_1_2_a =  Isub[0][1][0];
	complex<double> Ieta_xi_1_2_b =  0.5 * Isub[0][0][1] - 0.5 * Isub[0][1][1] - 0.5 * sqrt(3.0) * Isub[0][2][1];
	complex<double> Ieta_xi_1_2_c = -0.5 * Isub[0][0][2] - 0.5 * Isub[0][1][2] + 0.5 * sqrt(3.0) * Isub[0][2][2];

	Ieta_xi[0][1] = Ieta_xi_1_2_a + Ieta_xi_1_2_b + Ieta_xi_1_2_c;
	// Ieta_xi(1,3)
	complex<double> Ieta_xi_1_3_a =  Isub[0][2][0];
	complex<double> Ieta_xi_1_3_b =  0.5 * sqrt(3.0) * Isub[0][0][1] + 0.5 * sqrt(3.0) * Isub[0][1][1] - 0.5 * Isub[0][2][1];
	complex<double> Ieta_xi_1_3_c =  0.5 * sqrt(3.0) * Isub[0][0][2] - 0.5 * sqrt(3.0) * Isub[0][1][2] - 0.5 * Isub[0][2][2];

	Ieta_xi[0][2] = Ieta_xi_1_3_a + Ieta_xi_1_3_b + Ieta_xi_1_3_c;
	// Ieta_xi(2,1)
	complex<double> Ieta_xi_2_1_a =  Isub[1][0][0];
	complex<double> Ieta_xi_2_1_b =  0.5 * Isub[0][0][1] - 0.5 * Isub[1][0][1] - 0.5 * sqrt(3.0) * Isub[2][0][1];
	complex<double> Ieta_xi_2_1_c = -0.5 * Isub[0][0][2] - 0.5 * Isub[1][0][2] + 0.5 * sqrt(3.0) * Isub[2][0][2];
	
	Ieta_xi[1][0] = Ieta_xi_2_1_a + Ieta_xi_2_1_b + Ieta_xi_2_1_c;
	// Ieta_xi(2,2)
	complex<double> Ieta_xi_2_2_a =  Isub[1][1][0];
	complex<double> Ieta_xi_2_2_b =  0.25 *             Isub[0][0][1] - 0.25 *             Isub[0][1][1] - 0.25 * sqrt(3.0) * Isub[0][2][1]
	                               - 0.25 *             Isub[1][0][1] + 0.25 *             Isub[1][1][1] + 0.25 * sqrt(3.0) * Isub[1][2][1]
								   - 0.25 * sqrt(3.0) * Isub[2][0][1] + 0.25 * sqrt(3.0) * Isub[2][1][1] + 0.75 *             Isub[2][2][1];

	complex<double> Ieta_xi_2_2_c =  0.25 *             Isub[0][0][2] + 0.25 *             Isub[0][1][2] - 0.25 * sqrt(3.0) * Isub[0][2][2]
	                               + 0.25 *             Isub[1][0][2] + 0.25 *             Isub[1][1][2] - 0.25 * sqrt(3.0) * Isub[1][2][2]
								   - 0.25 * sqrt(3.0) * Isub[2][0][2] - 0.25 * sqrt(3.0) * Isub[2][1][2] + 0.75 *             Isub[2][2][2];

	Ieta_xi[1][1] = Ieta_xi_2_2_a + Ieta_xi_2_2_b + Ieta_xi_2_2_c;
	// Ieta_xi(2,3)
	complex<double> Ieta_xi_2_3_a =  Isub[1][2][0];
	complex<double> Ieta_xi_2_3_b =  0.25 * sqrt(3.0) * Isub[0][0][1] + 0.25 * sqrt(3.0) * Isub[0][1][1] - 0.25             * Isub[0][2][1]
	                               - 0.25 * sqrt(3.0) * Isub[1][0][1] - 0.25 * sqrt(3.0) * Isub[1][1][1] + 0.25 *             Isub[1][2][1]
								   - 0.75 *             Isub[2][0][1] - 0.75 *             Isub[2][1][1] + 0.25 * sqrt(3.0) * Isub[2][2][1];

	complex<double> Ieta_xi_2_3_c =- 0.25 * sqrt(3.0) * Isub[0][0][2] + 0.25 * sqrt(3.0) * Isub[0][1][2] + 0.25 *             Isub[0][2][2]
	                               - 0.25 * sqrt(3.0) * Isub[1][0][2] + 0.25 * sqrt(3.0) * Isub[1][1][2] + 0.25 *             Isub[1][2][2]
								   + 0.75 *             Isub[2][0][2] - 0.75 *             Isub[2][1][2] - 0.25 * sqrt(3.0) * Isub[2][2][2];

	Ieta_xi[1][2] = Ieta_xi_2_3_a + Ieta_xi_2_3_b + Ieta_xi_2_3_c;
	// Ieta_xi(3,1)
	complex<double> Ieta_xi_3_1_a =  Isub[2][0][0];
	complex<double> Ieta_xi_3_1_b =  0.5 * sqrt(3.0) * Isub[0][0][1] + 0.5 * sqrt(3.0) * Isub[1][0][1] - 0.5 * Isub[2][0][1];
	complex<double> Ieta_xi_3_1_c =  0.5 * sqrt(3.0) * Isub[0][0][2] - 0.5 * sqrt(3.0) * Isub[1][0][2] - 0.5 * Isub[2][0][2];
	
	Ieta_xi[2][0] = Ieta_xi_3_1_a + Ieta_xi_3_1_b + Ieta_xi_3_1_c;
	// Ieta_xi(3,2)
	complex<double> Ieta_xi_3_2_a =  Isub[2][1][0];
	complex<double> Ieta_xi_3_2_b =  0.25 * sqrt(3.0) * Isub[0][0][1] - 0.25 * sqrt(3.0) * Isub[0][1][1] - 0.75 *             Isub[0][2][1]
	                               + 0.25 * sqrt(3.0) * Isub[1][0][1] - 0.25 * sqrt(3.0) * Isub[1][1][1] - 0.75 *             Isub[1][2][1]
								   - 0.25 *             Isub[2][0][1] + 0.25 *             Isub[2][1][1] + 0.25 * sqrt(3.0) * Isub[2][2][1];

	complex<double> Ieta_xi_3_2_c =- 0.25 * sqrt(3.0) * Isub[0][0][2] - 0.25 * sqrt(3.0) * Isub[0][1][2] + 0.75 *             Isub[0][2][2]
	                               + 0.25 * sqrt(3.0) * Isub[1][0][2] + 0.25 * sqrt(3.0) * Isub[1][1][2] - 0.75 *             Isub[1][2][2]
								   + 0.25 *             Isub[2][0][2] + 0.25 *             Isub[2][1][2] - 0.25 * sqrt(3.0) * Isub[2][2][2];

	Ieta_xi[2][1] = Ieta_xi_3_2_a + Ieta_xi_3_2_b + Ieta_xi_3_2_c;
	// Ieta_xi(3,3)
	complex<double> Ieta_xi_3_3_a =  Isub[2][2][0];
	complex<double> Ieta_xi_3_3_b =  0.75 *             Isub[0][0][1] + 0.75 *             Isub[0][1][1] - 0.25 * sqrt(3.0) * Isub[0][2][1]
	                               + 0.75 *             Isub[1][0][1] + 0.75 *             Isub[1][1][1] - 0.25 * sqrt(3.0) * Isub[1][2][1]
								   - 0.25 * sqrt(3.0) * Isub[2][0][1] - 0.25 * sqrt(3.0) * Isub[2][1][1] + 0.25 *             Isub[2][2][1];

	complex<double> Ieta_xi_3_3_c =  0.75 *             Isub[0][0][2] - 0.75 *             Isub[0][1][2] - 0.25 * sqrt(3.0) * Isub[0][2][2]
	                               - 0.75 *             Isub[1][0][2] + 0.75 *             Isub[1][1][2] + 0.25 * sqrt(3.0) * Isub[1][2][2]
								   - 0.25 * sqrt(3.0) * Isub[2][0][2] + 0.25 * sqrt(3.0) * Isub[2][1][2] + 0.25 *             Isub[2][2][2];

	Ieta_xi[2][2] = Ieta_xi_3_3_a + Ieta_xi_3_3_b + Ieta_xi_3_3_c;
	// ***********************************************************************

	// 
	I_simplex[9] = Ieta_xi[0][0];

	// 
	complex<double> Ising_1_1 =  3.0 *       Ieta_xi[0][0] - 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           - 3.0 *       Ieta_xi[1][0] + 3.0 *       Ieta_xi[1][1] + sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[0] = 1.0 / 12.0 * Ising_1_1;
	// 
	complex<double> Ising_1_2 =  3.0 *       Ieta_xi[0][0] + 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           - 3.0 *       Ieta_xi[1][0] - 3.0 *       Ieta_xi[1][1] + sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[1] = 1.0 / 12.0 * Ising_1_2;
	// 
	complex<double> Ising_1_3 =  sqrt(3.0) * Ieta_xi[0][2] - sqrt(3.0) * Ieta_xi[1][2] -             Ieta_xi[2][2];

	I_simplex[2] = 1.0 / 6.0 * Ising_1_3;
	// 
	complex<double> Ising_2_1 =  3.0 *       Ieta_xi[0][0] - 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           + 3.0 *       Ieta_xi[1][0] - 3.0 *       Ieta_xi[1][1] - sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[3] = 1.0 / 12.0 * Ising_2_1;
	// 
	complex<double> Ising_2_2 =  3.0 *       Ieta_xi[0][0] + 3.0 *       Ieta_xi[0][1] - sqrt(3.0) * Ieta_xi[0][2]
	                           + 3.0 *       Ieta_xi[1][0] + 3.0 *       Ieta_xi[1][1] - sqrt(3.0) * Ieta_xi[1][2]
							   - sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] +             Ieta_xi[2][2];

	I_simplex[4] = 1.0 / 12.0 * Ising_2_2;
	// 
	complex<double> Ising_2_3 =  sqrt(3.0) * Ieta_xi[0][2] + sqrt(3.0) * Ieta_xi[1][2] -             Ieta_xi[2][2];

	I_simplex[5] = 1.0 / 6.0 * Ising_2_3;
	// 
	complex<double> Ising_3_1 =  sqrt(3.0) * Ieta_xi[2][0] - sqrt(3.0) * Ieta_xi[2][1] -             Ieta_xi[2][2];

	I_simplex[6] = 1.0 / 6.0 * Ising_3_1;
	// 
	complex<double> Ising_3_2 =  sqrt(3.0) * Ieta_xi[2][0] + sqrt(3.0) * Ieta_xi[2][1] -             Ieta_xi[2][2];

	I_simplex[7] = 1.0 / 6.0 * Ising_3_2;
	// 
	I_simplex[8] = 1.0 / 3.0 * Ieta_xi[2][2];
}
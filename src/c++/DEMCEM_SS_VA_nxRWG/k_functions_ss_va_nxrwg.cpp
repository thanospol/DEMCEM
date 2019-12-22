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
#include "demcem_constants.h"

#include <cmath>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void X_function_pre
// ***********************************************************************

void k_functions_ss_va_nxrwg (double GAMMA, double L, complex<double> K[], const complex<double> ko)
{
	complex<double> a  = Iunit * ko * GAMMA;
	//
	K[0] = (1.0 / pow(a,2.0) ) * (1.0 - exp(-a * L) - a * L * exp(-a * L) );
	K[1] = (1.0 / pow(a,3.0) ) * (2.0 - 2.0 * exp(-a * L) - 2.0 * a * L * exp(-a * L) - pow(a,2.0) * pow(L,2.0) * exp(-a * L) );
	K[2] = -(pow(L,3.0) * exp(-a * L) ) / a + (3.0 / a) * K[1];
	K[3] = -(pow(L,4.0) * exp(-a * L) ) / a + (4.0 / a) * K[2];
}
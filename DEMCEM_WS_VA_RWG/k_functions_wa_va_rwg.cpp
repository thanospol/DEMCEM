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
#include "demcem_ws_va_rwg.h"
#include "demcem_constants.h"


using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void k_functions
// ***********************************************************************

void k_functions_wa_va_rwg (double GAMMA, double L, complex<double> K[], const double ko)
{
	complex<double> a  = Iunit * ko * GAMMA;
	//
	K[0] = (1.0 / pow(a,2.0) ) * (1.0 - exp(-a * L) - a * L * exp(-a * L) );
	K[1] = (1.0 / pow(a,3.0) ) * (2.0 - 2.0 * exp(-a * L) - 2.0 * a * L * exp(-a * L) - pow(a,2.0) * pow(L,2.0) * exp(-a * L) );
	K[2] = -(pow(L,3.0) * exp(-a * L) ) / a + (3.0 / a) * K[1];
	K[3] = -(pow(L,4.0) * exp(-a * L) ) / a + (4.0 / a) * K[2];
}
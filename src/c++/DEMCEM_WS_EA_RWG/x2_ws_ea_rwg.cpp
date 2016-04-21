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


// ***********************************************************************
//			IMPLEMENTATION OF void void x2
// ***********************************************************************

void x2_ws_ea_rwg ( complex<double> D, complex<double> a, complex<double> N[] )
{
	complex<double> aD, expaD, T1, T2, T3, T4;
	//
	aD    = a * D;
	expaD = exp(-aD);
	//
	T1    = (1.0 - expaD) / aD;
	T2    = (1.0 - T1) / aD;
	T3    = (1.0 - 2.0 * T2) / aD;
	T4    = (1.0 - 3.0 * T3) / aD;
	// Final Output
	N[0]   = 1.0;
	N[1]   = T1;
	N[2]   = D * (T1 - T2);
	N[3]   = 0.0;
	N[4]   = 0.0;
	N[5]   = 0.5;
	N[6]   = 1.0 / 3.0;
	N[7]   = T2;
	N[8]   = T3;
	N[9]   = D * (T2 - T3);
	N[10]  = D * (T3 - T4);
	N[11]  = 0.0;
	//
	N[11]  = D * (N[9] - N[10]);
	N[3]   = D * (N[2] - N[9]);
	N[4]   = D * (N[3] - N[11]);
}
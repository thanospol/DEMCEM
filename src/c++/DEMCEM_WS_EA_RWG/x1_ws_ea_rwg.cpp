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
//			IMPLEMENTATION OF void void x1
// ***********************************************************************

void x1_ws_ea_rwg ( double psi, complex<double> D1, complex<double> D2, double tpsiB, complex<double> a, complex<double> N[] )
{
	complex<double> D3, aD3, expaD3, aD1, expaD1, D2a, H, T11, T21, T31, T41, T12, T22, T32, T42;
	complex<double>	N21, N31, N81, N91, N101, N111, N121, N41, N51, N22, N32, N82, N92, N102, N112, N122, N42, N52;
	//
	D3     = sqrt(3.0) / (cos(psi) * tpsiB);
	aD3    = a * D3;
	expaD3 = exp(-aD3);
	aD1    = a * D1;
	expaD1 = exp(-aD1);
	D2a    = D2 / a;
	//
	H      = 1.0 - D1 * D2;
	//
	T11    = (expaD3 - expaD1) / aD3;
	T21    = (T11 - H * expaD1) / aD3;
	T31    = (2.0 * T21 - pow(H,2) * expaD1) / aD3;
	T41    = (3.0 * T31 - pow(H,3) * expaD1) / aD3;
	//
	T12    = D2a * (1.0 - expaD1);
	T22    = D2a * (1.0 - H * expaD1 - T12);
	T32    = D2a * (1.0 - pow(H,2) * expaD1 - 2.0 * T22);
	T42    = D2a * (1.0 - pow(H,3) * expaD1 - 3.0 * T32);
	//
	N21    = T11;
	N31    = D3 * (T11 + T21);
	N81    = T21;
	N91    = T31;
	N101   = D3 * (T21 + T31);
	N111   = D3 * (T31 + T41);
	N121   = D3 * (N101 + N111);
	N41    = D3 * (N31 + N101);
	N51    = D3 * (N41 + N121);
	//
	N22    = T12;
	N32    = (T12 - T22) / D2;
	N82    = T22;
	N92    = T32;
	N102   = (T22 - T32) / D2;
	N112   = (T32 - T42) / D2;
	N122   = (N102 - N112) / D2;
	N42    = (N32 - N102) / D2;
	N52    = (N42 - N122) / D2;
	// Final Output
	N[0]   = 1.0;
	N[1]   = N21+N22;
	N[2]   = N31+N32;
	N[3]   = N41+N42;
	N[4]   = N51+N52;
	N[5]   = 0.5;
	N[6]   = 1.0 / 3.0;
	N[7]   = N81+N82;
	N[8]   = N91+N92;
	N[9]   = N101+N102;
	N[10]  = N111+N112;
	N[11]  = N121+N122;
}
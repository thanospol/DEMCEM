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
#include "demcem_constants.h"
#include <cmath>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void X_function_pre
// ***********************************************************************

void x_function_pre_ss_ea_nxrwg (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const double ko)
{
	double  D, D1, D2;
	//
	if (Psi >= PsiB && Psi >= PsiA)
	{
		D  = sqrt(3.0) / sin(Psi);
		x2_ss_ea_nxrwg(D, ko, B, N);
		x2_ss_ea_nxrwg (D, ko, Bm, Nm);
	}
	else if (theta >= M_PI_2)
	{
		D  = sqrt(3.0) / (cos(Psi) * tPsiA);
		x2_ss_ea_nxrwg (D, ko, B, N);
		x2_ss_ea_nxrwg (D, ko, Bm, Nm);
	}
	else if (Psi >= PsiA)
	{
		D1 = 2.0 * sqrt(3.0) / (cos(Psi) * ( tPsiB + tan(Psi) )  );
		D2 = sin(Psi) / sqrt(3.0);
		x1_ss_ea_nxrwg( Psi, D1, D2, tPsiB, ko, B,  N);
		x1_ss_ea_nxrwg( Psi, D1, D2, tPsiB, ko, Bm, Nm);
	}
	else
	{
		D1 = sqrt(3.0) / (cos(Psi) * sin(theta) );
		D2 = ( cos(Psi) * tPsiA ) / sqrt(3.0);
		x1_ss_ea_nxrwg( Psi, D1, D2, tPsiB, ko, B, N);
		x1_ss_ea_nxrwg( Psi, D1, D2, tPsiB, ko, Bm, Nm);
	}

}
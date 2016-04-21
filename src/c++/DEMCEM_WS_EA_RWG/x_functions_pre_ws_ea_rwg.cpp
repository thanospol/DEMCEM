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
#include "demcem_constants.h"


// ***********************************************************************
//			IMPLEMENTATION OF void x_function_pre
// ***********************************************************************

void x_functions_pre_ws_ea_rwg (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const double ko)
{
	double  D, D1, D2;
	//
	complex<double> a  = Iunit * ko * B;
	complex<double> am = Iunit * ko * Bm;
	//
	if (Psi >= PsiB && Psi >= PsiA)
	{
		D  = sqrt(3.0) / sin(Psi);
		x2_ws_ea_rwg(D, a, N);
		x2_ws_ea_rwg (D, am, Nm);
	}
	else if (theta >= M_PI_2)
	{
		D  = sqrt(3.0) / (cos(Psi) * tPsiA);
		x2_ws_ea_rwg (D, a, N);
		x2_ws_ea_rwg (D, am, Nm);
	}
	else if (Psi >= PsiA)
	{
		D1 = 2.0 * sqrt(3.0) / (cos(Psi) * ( tPsiB + tan(Psi) )  );
		D2 = sin(Psi) / sqrt(3.0);
		x1_ws_ea_rwg(Psi, D1, D2, tPsiB, a, N);
		x1_ws_ea_rwg(Psi, D1, D2, tPsiB, am, Nm);
	}
	else
	{
		D1 = sqrt(3.0) / (cos(Psi) * sin(theta) );
		D2 = ( cos(Psi) * tPsiA ) / sqrt(3.0);
		x1_ws_ea_rwg(Psi, D1, D2, tPsiB, a, N);
		x1_ws_ea_rwg(Psi, D1, D2, tPsiB, am, Nm);
	}

}
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


#include "demcem_ss_ea_rwg.h"

// ***********************************************************************
//			IMPLEMENTATION OF void psi_limits
// ***********************************************************************

void psi_limits_ss_ea_rwg ( int argument, double PsiA, double PsiB, double *psi_A, double *psi_B )
{
	switch(argument)
	{
		case 1:
			 *psi_A = PsiB;
			 *psi_B = M_PI_2;
			 break;
		case 2:
			 *psi_A = PsiA;
			 *psi_B = M_PI_2;
			 break;
		case 3:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
		case 4:
			 *psi_A = 0;
			 *psi_B = PsiB;
			 break;
		case 5:
			 *psi_A = PsiA;
			 *psi_B = PsiB;
			 break;
		case 6:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
	}
}
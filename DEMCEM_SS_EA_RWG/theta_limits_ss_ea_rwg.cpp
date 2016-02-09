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

void theta_limits_ss_ea_rwg ( int argument, double *theta_A, double *theta_B )
{
	switch(argument)
	{
		case 1:
			 *theta_A = 0.0;
			 *theta_B = M_PI_2 ;
			 break;
		case 2:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 3:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 4:
			 *theta_A = 0.0;
			 *theta_B = M_PI / 3;
			 break;
		case 5:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
		case 6:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
	}
}
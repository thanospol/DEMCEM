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
#include "demcem_constants.h"

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_functions
// ***********************************************************************

complex<double> phi_functions_ws_st_rwg ( int ROW, int COLUMN, int argument, double psi, double Apsi, const complex<double> ko )
{
	complex<double> PHI_out;
	complex<double> j   = Iunit;

	complex<double> PSI_1 ;
	complex<double> PSI_2 ;
	complex<double> PSI_3 ;
	complex<double> PSI_ ;

	complex<double> A = j * ko * Apsi;
	complex<double> B = cos(psi);
	complex<double> C = sin(psi) / sqrt(3.0);

	// Int_1_1 
	if (ROW == 0  && COLUMN == 0)
	{
		                                     
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 11);
		
	}// end if (ROW == 0  && COLUMN == 0)
	// Int_1_2
	else if (ROW == 0  && COLUMN == 1)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 12 );
		
	}// end if (ROW == 0  && COLUMN == 1)
	// Int_1_3
	else if (ROW == 0  && COLUMN == 2)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 13 );
		
	}// end if (ROW == 0  && COLUMN == 2)
	// Int_2_1
	else if (ROW == 1  && COLUMN == 0)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 21 );
		
	}// end if (ROW == 1  && COLUMN == 0)
	// Int_2_2
	else if (ROW == 1  && COLUMN == 1)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 22 );
		
	}// end if (ROW == 1  && COLUMN == 1)
	// Int_2_3
	else if (ROW == 1  && COLUMN == 2)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 23 );
		
	}// end if (ROW == 1  && COLUMN == 2)
	// Int_3_1
	else if (ROW == 2  && COLUMN == 0)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 31 );
		
	}// end if (ROW == 2  && COLUMN == 0)
	// Int_3_2
	else if (ROW == 2  && COLUMN == 1)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 32 );
		
	}// end if (ROW == 2  && COLUMN == 1)
	// Int_3_3
	else if (ROW == 2  && COLUMN == 2)
	{
		                                      
		PHI_out = phi_ws_st_rwg (  A,  B,  C, argument,  psi,  Apsi,  ko, 33 );
		
	}// end if (ROW == 2  && COLUMN == 2)
	// Final Output
	return PHI_out;
}
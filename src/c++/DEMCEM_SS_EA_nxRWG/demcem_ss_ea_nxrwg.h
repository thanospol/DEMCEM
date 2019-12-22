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


#if !defined _DEMCEM_SS_EA_NXRGW_H_
#define _DEMCEM_SS_EA_NXRGW_H_

#include <cmath>
#include <complex>

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************

// Common functions for all cases

void x2_ss_ea_nxrwg  ( complex<double> D, const complex<double> ko, complex<double> a, complex<double> N[] );
void x1_ss_ea_nxrwg  ( double psi, complex<double> D1, complex<double> D2, double tpsiB, const complex<double> ko, complex<double> B, complex<double> N[] );
void gl_quad ( int n, double x[], double w[] );
void theta_limits_ss_ea_nxrwg ( int argument, double *theta_A, double *theta_B );
void psi_limits_ss_ea_nxrwg ( int argument, double PsiA, double PsiB, double *psi_A, double *psi_B );
void x_function_pre_ss_ea_nxrwg (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const complex<double> ko);

// Case-dependent functions
void coefficients_ss_ea_nxrwg (const double r1[],const double r2[],const double r3[],const double r4[],const complex<double> ko, double Ap, complex<double> coef[], complex<double> coefm[], int flag );

complex<double> x_functions_ss_ea_nxrwg (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[], int flag);



#endif
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


#if !defined _DEMCEM_WS_EA_RWG_H_
#define _DEMCEM_WS_EA_RWG_H_

#include <cmath>
#include <complex>

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************

void x2_ws_ea_rwg  ( complex<double> D, complex<double> a, complex<double> N[] );

void x1_ws_ea_rwg  ( double psi, complex<double> D1, complex<double> D2, double tpsiB, complex<double> a, complex<double> N[] );

void gl_quad ( int n, double x[], double w[] );

void theta_limits_ws_ea_rwg ( int argument, double *theta_A, double *theta_B );

void psi_limits_ws_ea_rwg ( int argument, double PsiA, double PsiB, double *psi_A, double *psi_B );

void x_functions_pre_ws_ea_rwg (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const double ko);

void coefficients_ws_ea_rwg (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[], int flag );

complex<double> x_functions_ws_ea_rwg (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[], int flag);


#endif 
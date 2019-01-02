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

#if !defined _DEMCEM_SS_VA_RWG_H_
#define _DEMCEM_SS_VA_RWG_H_

#include <cmath>
#include <complex>

using namespace std;

/* DECLARATION OF FUNCTIONS */

void gl_quad ( int n, double x[], double w[] );

void k_functions_ss_va_rwg (double GAMMA, double L, complex<double> K[], const complex<double> ko);

void coefficients_ss_va_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const complex<double> ko, complex<double> coef[], int flag );

complex<double> omega_functions_ss_va_rwg (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const complex<double> ko, complex<double> K[], int flag);

#endif 
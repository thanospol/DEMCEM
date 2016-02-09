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

#if !defined _DEMCEM_WS_ST_RWG_H_
#define _DEMCEM_WS_ST_RWG_H_

#include <cmath>
#include <complex>


using namespace std;

/* DECLARATION OF FUNCTIONS */

void post_ws_st_rwg ( complex<double> Isub[3][3][3], complex<double> I_simplex[]);

complex<double> phi_functions_ws_st_rwg ( int ROW, int COLUMN, int argument, double psi, double Apsi, const double ko);

complex<double> phi_ws_st_rwg ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko , int flag);

void gl_quad ( int n, double x[], double w[] );

#endif 
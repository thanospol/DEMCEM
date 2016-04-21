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
#include "demcem_ws_va_rwg.h"

using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> omega_functions
// ***********************************************************************


complex<double> omega_functions_wa_va_rwg (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], complex<double> K[], int flag)
{                                 
	complex<double> X;
	//
    if (flag == 0)
    {
        X = K[1] * (sin(Psi) * cos(Psi) * coef[0]) / GAMMA;
    }// if (flag == 0)
    if (flag == 11)
    {
        double t2 = sin(Psi);
        double t3 = t2 * t2;
        double t4 = cos(Psi);
        double t5 = t4 * t4;
        double t8 = t3 * t5 / GAMMA;
        double t9 = sin(theta_p);
        double t10 = sin(theta_q);
        double t15 = cos(theta_q);
        double t16 = cos(theta_p);
        //
        X = K[3] * (t8 * t9 * t10 * coef[0] + t8 * t15 * t16 * coef[1] + t8 * t10 * t16 * coef[2] + t8 * t15 * t9 * coef[3]);
    }// if (flag == 11)
    
    if (flag == 12)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_p);
        complex<double> t12 = sin(theta_p);
        complex<double> t20 = t2 * t2;
        complex<double> t22 = t20 * t4 * t6;
        complex<double> t23 = sin(theta_q);
        complex<double> t28 = cos(theta_q);
        //
        X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t12 * t23 * coef[0] + t22 * t28 * t7 * coef[3] + t22 * t28 * t12 * coef[4] + t22 * t23 * t7 * coef[5]);
        
    }// if (flag == 12)
    
    if (flag == 13)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_p);
        complex<double> t12 = sin(theta_p);
        complex<double> t20 = t2 * t2;
        complex<double> t22 = t20 * t4 * t6;
        complex<double> t23 = sin(theta_q);
        complex<double> t28 = cos(theta_q);
        //
        X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t12 * t23 * coef[0] + t22 * t28 * t7 * coef[3] + t22 * t28 * t12 * coef[4] + t22 * t23 * t7 * coef[5]);
        
    }// if (flag == 13)
    
    if (flag == 21)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t3 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_q);
        complex<double> t12 = sin(theta_q);
        complex<double> t20 = t4 * t4;
        complex<double> t22 = t3 * t20 * t6;
        complex<double> t23 = sin(theta_p);
        complex<double> t28 = cos(theta_p);
        //
        X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[5]) + K[3] * (t22 * t23 * t12 * coef[0] + t22 * t7 * t28 * coef[2] + t22 * t12 * t28 * coef[3] + t22 * t7 * t23 * coef[4]);
        
    }// if (flag == 21)
    
    if (flag == 22)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t3 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_q);
        complex<double> t12 = sin(theta_q);
        complex<double> t17 = t4 * t4;
        complex<double> t18 = t2 * t17;
        complex<double> t19 = cos(theta_p);
        complex<double> t24 = sin(theta_p);
        complex<double> t39 = t3 * t17 * t6;
        //
        X = K[2] * (t5 * t6 * t7 * coef[4] + t5 * t6 * t12 * coef[5] + t18 * t6 * t19 * coef[2] + t18 * t6 * t24 * coef[3]) + K[1] * t2 * t4 * t6 * coef[0] + K[3] * (t39 * t24 * t12 * coef[1] + t39 * t19 * t12 * coef[6] + t39 * t7 * t24 * coef[7] + t39 * t7 * t19 * coef[8]);
        
    }// if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_p);
        complex<double> t12 = sin(theta_p);
        complex<double> t17 = t2 * t2;
        complex<double> t18 = t17 * t3;
        complex<double> t19 = cos(theta_q);
        complex<double> t24 = sin(theta_q);
        complex<double> t39 = t17 * t4 * t6;
        //
        X = K[2] * (t5 * t6 * t7 * coef[8] + t5 * t6 * t12 * coef[5] + t18 * t6 * t19 * coef[6] + t18 * t6 * t24 * coef[7]) + K[1] * t2 * t3 * t6 * coef[0] + K[3] * (t39 * t12 * t24 * coef[1] + t39 * t19 * t7 * coef[2] + t39 * t19 * t12 * coef[3] + t39 * t24 * t7 * coef[4]);
        
    }// if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t3 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = sin(theta_q);
        complex<double> t12 = cos(theta_q);
        complex<double> t20 = t4 * t4;
        complex<double> t22 = t3 * t20 * t6;
        complex<double> t23 = sin(theta_p);
        complex<double> t28 = cos(theta_p);
        //
        X = K[2] * (t5 * t6 * t7 * coef[1] + t5 * t6 * t12 * coef[2]) + K[3] * (t22 * t23 * t7 * coef[0] + t22 * t12 * t28 * coef[3] + t22 * t7 * t28 * coef[4] + t22 * t12 * t23 * coef[5]);
        
    }// if (flag == 31)
    
    if (flag == 32)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = sin(theta_p);
        complex<double> t12 = t2 * t2;
        complex<double> t13 = t12 * t3;
        complex<double> t14 = cos(theta_q);
        complex<double> t19 = cos(theta_p);
        complex<double> t24 = sin(theta_q);
        complex<double> t39 = t12 * t4 * t6;
        //
        X = K[2] * (t5 * t6 * t7 * coef[8] + t13 * t6 * t14 * coef[2] + t5 * t6 * t19 * coef[7] + t13 * t6 * t24 * coef[6]) + K[1] * t2 * t3 * t6 * coef[0] + K[3] * (t39 * t7 * t24 * coef[1] + t39 * t14 * t7 * coef[5] + t39 * t14 * t19 * coef[3] + t39 * t24 * t19 * coef[4]);
        
    }// if (flag == 32)
    
    if (flag == 33)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t3 * t4;
        complex<double> t6 = 0.1e1 / GAMMA;
        complex<double> t7 = cos(theta_q);
        complex<double> t12 = sin(theta_q);
        complex<double> t17 = t4 * t4;
        complex<double> t18 = t2 * t17;
        complex<double> t19 = sin(theta_p);
        complex<double> t24 = cos(theta_p);
        complex<double> t39 = t3 * t17 * t6;
        //
        X = K[2] * (t5 * t6 * t7 * coef[8] + t5 * t6 * t12 * coef[2] + t18 * t6 * t19 * coef[3] + t18 * t6 * t24 * coef[7]) + K[1] * t2 * t4 * t6 * coef[0] + K[3] * (t39 * t19 * t12 * coef[1] + t39 * t12 * t24 * coef[5] + t39 * t7 * t19 * coef[6] + t39 * t7 * t24 * coef[4]);
        
    }// if (flag == 33)
    
	// Final Output
	return X;
}



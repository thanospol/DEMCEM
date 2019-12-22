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
#include "demcem_ss_va_rwg.h"
#include "demcem_constants.h"


using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> omega_functions
// ***********************************************************************


complex<double> omega_functions_ss_va_rwg (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const complex<double> ko, complex<double> K[], int flag)

{                                 
	complex<double> X;
	//
    complex<double> j = Iunit;
    
    if (flag == 12)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t10 = t6 / t7 / GAMMA;
        complex<double> t11 = sin(theta_p);
        complex<double> t12 = sin(theta_q);
        complex<double> t14 = coef[0];
        complex<double> t17 = cos(theta_p);
        complex<double> t19 = coef[1];
        complex<double> t27 = t6 / t7 * j;
        //
        X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t17 * t12 * t19) + K[2] * (-t27 * ko * t11 * t12 * t14 - t27 * ko * t17 * t12 * t19);
        
    }// if (flag == 12)
    
    if (flag == 13)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t10 = t6 / t7 / GAMMA;
        complex<double> t11 = sin(theta_q);
        complex<double> t12 = cos(theta_p);
        complex<double> t14 = coef[0];
        complex<double> t17 = sin(theta_p);
        complex<double> t19 = coef[1];
        complex<double> t22 = cos(theta_q);
        complex<double> t24 = coef[2];
        complex<double> t28 = coef[3];
        complex<double> t36 = t6 / t7 * j;
        complex<double> t37 = ko * t11;
        complex<double> t44 = ko * t22;
        //
        X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t11 * t17 * t19 - t10 * t22 * t12 * t24 - t10 * t22 * t17 * t28) + K[2] * (-t36 * t37 * t12 * t14 - t36 * t37 * t17 * t19 - t36 * t44 * t12 * t24 - t36 * t44 * t17 * t28);
        
    }// if (flag == 13)
    
    if (flag == 21)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t10 = t6 / t7 / GAMMA;
        complex<double> t11 = sin(theta_p);
        complex<double> t12 = sin(theta_q);
        complex<double> t14 = coef[0];
        complex<double> t17 = cos(theta_q);
        complex<double> t19 = coef[1];
        complex<double> t27 = t6 / t7 * j;
        //
        X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t17 * t11 * t19) + K[2] * (-t27 * ko * t11 * t12 * t14 - t27 * ko * t17 * t11 * t19);
        
    }// if (flag == 21)
    
    if (flag == 22)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = pow(GAMMA, 0.2e1);
        complex<double> t8 = 0.1e1 / t6 / GAMMA;
        complex<double> t9 = sin(theta_p);
        complex<double> t11 = coef[1];
        complex<double> t14 = t2 * t2;
        complex<double> t15 = t14 * t3;
        complex<double> t16 = sin(theta_q);
        complex<double> t18 = coef[2];
        complex<double> t24 = t14 * t4;
        complex<double> t25 = t24 * t8;
        complex<double> t27 = coef[0];
        complex<double> t30 = cos(theta_p);
        complex<double> t32 = coef[3];
        complex<double> t35 = cos(theta_q);
        complex<double> t37 = coef[4];
        complex<double> t40 = 0.1e1 / t6;
        complex<double> t42 = j * ko;
        complex<double> t54 = t24 * t40 * j;
        //
        X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18) + K[1] * (-t25 * t9 * t16 * t27 - t25 * t30 * t16 * t32 - t25 * t35 * t9 * t37 - t5 * t40 * t42 * t9 * t11 - t15 * t40 * t42 * t16 * t18) + K[2] * (-t54 * ko * t9 * t16 * t27 - t54 * ko * t30 * t16 * t32 - t54 * ko * t35 * t9 * t37);
        
    }// if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t3 * t4;
        complex<double> t6 = pow(GAMMA, 0.2e1);
        complex<double> t8 = 0.1e1 / t6 / GAMMA;
        complex<double> t9 = sin(theta_q);
        complex<double> t11 = coef[1];
        complex<double> t14 = t4 * t4;
        complex<double> t15 = t2 * t14;
        complex<double> t16 = sin(theta_p);
        complex<double> t18 = coef[2];
        complex<double> t21 = cos(theta_q);
        complex<double> t23 = coef[3];
        complex<double> t29 = t3 * t14;
        complex<double> t30 = t29 * t8;
        complex<double> t32 = coef[0];
        complex<double> t35 = cos(theta_p);
        complex<double> t37 = coef[4];
        complex<double> t41 = coef[5];
        complex<double> t45 = coef[6];
        complex<double> t48 = 0.1e1 / t6;
        complex<double> t49 = t5 * t48;
        complex<double> t50 = j * ko;
        complex<double> t65 = t29 * t48 * j;
        complex<double> t70 = ko * t21;
        //
        X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18 - t5 * t8 * t21 * t23) + K[1] * (-t30 * t16 * t9 * t32 - t30 * t21 * t35 * t37 - t30 * t9 * t35 * t41 - t30 * t21 * t16 * t45 - t49 * t50 * t9 * t11 - t15 * t48 * t50 * t16 * t18 - t49 * t50 * t21 * t23) + K[2] * (-t65 * ko * t16 * t9 * t32 - t65 * t70 * t35 * t37 - t65 * ko * t9 * t35 * t41 - t65 * t70 * t16 * t45);
        
    }// if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t10 = t6 / t7 / GAMMA;
        complex<double> t11 = sin(theta_q);
        complex<double> t12 = cos(theta_p);
        complex<double> t14 = coef[0];
        complex<double> t17 = sin(theta_p);
        complex<double> t19 = coef[1];
        complex<double> t22 = cos(theta_q);
        complex<double> t24 = coef[2];
        complex<double> t28 = coef[3];
        complex<double> t36 = t6 / t7 * j;
        complex<double> t37 = ko * t11;
        complex<double> t44 = ko * t22;
        //
        X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t11 * t17 * t19 - t10 * t22 * t12 * t24 - t10 * t22 * t17 * t28) + K[2] * (-t36 * t37 * t12 * t14 - t36 * t37 * t17 * t19 - t36 * t44 * t12 * t24 - t36 * t44 * t17 * t28);
        
    }// if (flag == 31)
    
    if (flag == 32)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = pow(GAMMA, 0.2e1);
        complex<double> t8 = 0.1e1 / t6 / GAMMA;
        complex<double> t9 = sin(theta_p);
        complex<double> t11 = coef[1];
        complex<double> t14 = cos(theta_p);
        complex<double> t16 = coef[2];
        complex<double> t19 = t2 * t2;
        complex<double> t20 = t19 * t3;
        complex<double> t21 = sin(theta_q);
        complex<double> t23 = coef[3];
        complex<double> t29 = t19 * t4;
        complex<double> t30 = t29 * t8;
        complex<double> t32 = coef[0];
        complex<double> t35 = cos(theta_q);
        complex<double> t37 = coef[4];
        complex<double> t41 = coef[5];
        complex<double> t45 = coef[6];
        complex<double> t48 = 0.1e1 / t6;
        complex<double> t49 = t5 * t48;
        complex<double> t50 = j * ko;
        complex<double> t65 = t29 * t48 * j;
        complex<double> t70 = ko * t35;
        //
        X = K[0] * (-t5 * t8 * t9 * t11 - t5 * t8 * t14 * t16 - t20 * t8 * t21 * t23) + K[1] * (-t30 * t9 * t21 * t32 - t30 * t35 * t14 * t37 - t30 * t21 * t14 * t41 - t30 * t35 * t9 * t45 - t49 * t50 * t9 * t11 - t49 * t50 * t14 * t16 - t20 * t48 * t50 * t21 * t23) + K[2] * (-t65 * ko * t9 * t21 * t32 - t65 * t70 * t14 * t37 - t65 * ko * t21 * t14 * t41 - t65 * t70 * t9 * t45);
        
    }// if (flag == 32)
    
    if (flag == 33)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t2 * t4;
        complex<double> t6 = pow(GAMMA, 0.2e1);
        complex<double> t8 = 0.1e1 / t6 / GAMMA;
        complex<double> t9 = sin(theta_p);
        complex<double> t11 = coef[0];
        complex<double> t14 = t2 * t2;
        complex<double> t15 = t14 * t3;
        complex<double> t16 = cos(theta_q);
        complex<double> t18 = coef[1];
        complex<double> t21 = sin(theta_q);
        complex<double> t23 = coef[2];
        complex<double> t26 = cos(theta_p);
        complex<double> t28 = coef[3];
        complex<double> t34 = t14 * t4;
        complex<double> t35 = t34 * t8;
        complex<double> t37 = coef[4];
        complex<double> t41 = coef[5];
        complex<double> t45 = coef[6];
        complex<double> t49 = coef[7];
        complex<double> t52 = 0.1e1 / t6;
        complex<double> t53 = t5 * t52;
        complex<double> t54 = j * ko;
        complex<double> t58 = t15 * t52;
        complex<double> t72 = t34 * t52 * j;
        complex<double> t73 = ko * t16;
        complex<double> t80 = ko * t21;
        //
        X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18 - t15 * t8 * t21 * t23 - t5 * t8 * t26 * t28) + K[1] * (-t35 * t16 * t9 * t37 - t35 * t16 * t26 * t41 - t35 * t21 * t9 * t45 - t35 * t21 * t26 * t49 - t53 * t54 * t9 * t11 - t58 * t54 * t16 * t18 - t58 * t54 * t21 * t23 - t53 * t54 * t26 * t28) + K[2] * (-t72 * t73 * t9 * t37 - t72 * t73 * t26 * t41 - t72 * t80 * t9 * t45 - t72 * t80 * t26 * t49);
        
    }// if (flag == 33)
    
	// Final Output
	return X;
}



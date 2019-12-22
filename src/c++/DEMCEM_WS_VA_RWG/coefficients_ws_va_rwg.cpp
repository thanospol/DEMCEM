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
//			IMPLEMENTATION OF void coefficients_va
// ***********************************************************************

void coefficients_wa_va_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const complex<double> ko, complex<double> coef[], int flag )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
    
    if (flag == 0)
    {
        coef[0] = 1.0;
    } // if (flag == 0)
    
    if (flag == 11)
    {
        complex<double> t1 = bb[0] * dd[0];
        complex<double> t3 = bb[1] * dd[1];
        complex<double> t5 = bb[2] * dd[2];
        complex<double> t7 = bb[0] * ee[0];
        complex<double> t9 = bb[1] * ee[1];
        complex<double> t11 = bb[2] * ee[2];
        complex<double> t13 = cc[0] * dd[0];
        complex<double> t15 = cc[1] * dd[1];
        complex<double> t17 = cc[2] * dd[2];
        complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t29 = t1 * t27 / 0.12e2;
        complex<double> t31 = t3 * t27 / 0.12e2;
        complex<double> t35 = t5 * t27 / 0.12e2;
        //
        coef[0] = t25;
        coef[1] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
        coef[2] = -t29 - t31 + t9 * t27 / 0.6e1 - t35 + t7 * t27 / 0.6e1 + t11 * t27 / 0.6e1;
        coef[3] = -t29 - t35 - t31 + t13 * t27 / 0.6e1 + t15 * t27 / 0.6e1 + t17 * t27 / 0.6e1;
        
    } // if (flag == 11)
    
    if (flag == 12)
    {
        complex<double> t1 = bb[1] * dd[1];
        complex<double> t3 = bb[0] * dd[0];
        complex<double> t5 = bb[2] * dd[2];
        complex<double> t7 = bb[1] * ee[1];
        complex<double> t9 = bb[2] * ee[2];
        complex<double> t11 = bb[0] * ee[0];
        complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - cc[0] * dd[0] / 0.6e1 - cc[1] * dd[1] / 0.6e1 - cc[2] * dd[2] / 0.6e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
        complex<double> t26 = -t3 - t1 - t5;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t29 = t27 * bb[2] * dd[2];
        complex<double> t32 = t27 * bb[1] * dd[1];
        complex<double> t35 = t27 * cc[0] * dd[0];
        complex<double> t38 = t27 * bb[0] * dd[0];
        complex<double> t41 = t27 * cc[1] * dd[1];
        complex<double> t44 = t27 * cc[2] * dd[2];
        complex<double> t47 = t32 / 0.12e2;
        complex<double> t49 = t29 / 0.12e2;
        complex<double> t50 = t38 / 0.12e2;
        //
        coef[0] = t25;
        coef[1] = t26 / 0.2e1;
        coef[2] = t29 / 0.6e1 + t32 / 0.6e1 - t35 / 0.3e1 + t38 / 0.6e1 - t41 / 0.3e1 - t44 / 0.3e1;
        coef[3] = -t26 / 0.4e1;
        coef[4] = -t47 + t41 / 0.6e1 - t49 - t50 + t44 / 0.6e1 + t35 / 0.6e1;
        coef[5] = -t50 - t47 + t9 * t27 / 0.6e1 - t49 + t11 * t27 / 0.6e1 + t7 * t27 / 0.6e1;
        
    }// if (flag == 12)
    
    if (flag == 13)
    {
        complex<double> t1 = bb[1] * dd[1];
        complex<double> t3 = bb[0] * dd[0];
        complex<double> t5 = bb[2] * dd[2];
        complex<double> t7 = bb[1] * ee[1];
        complex<double> t9 = bb[2] * ee[2];
        complex<double> t11 = bb[0] * ee[0];
        complex<double> t13 = cc[0] * dd[0];
        complex<double> t15 = cc[1] * dd[1];
        complex<double> t17 = cc[2] * dd[2];
        complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t30 = t27 * bb[2] * ee[2] / 0.6e1;
        complex<double> t33 = t27 * bb[1] * ee[1] / 0.6e1;
        complex<double> t39 = t27 * bb[0] * ee[0] / 0.6e1;
        complex<double> t49 = t1 * t27 / 0.12e2;
        complex<double> t51 = t5 * t27 / 0.12e2;
        complex<double> t55 = t3 * t27 / 0.12e2;
        //
        coef[0] = t25;
        coef[1] = -t11 / 0.2e1 - t7 / 0.2e1 - t9 / 0.2e1;
        coef[2] = t30 + t33 - t27 * cc[0] * ee[0] / 0.3e1 + t39 - t27 * cc[1] * ee[1] / 0.3e1 - t27 * cc[2] * ee[2] / 0.3e1;
        coef[3] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
        coef[4] = -t49 - t51 + t15 * t27 / 0.6e1 - t55 + t13 * t27 / 0.6e1 + t17 * t27 / 0.6e1;
        coef[5] = -t51 - t49 + t39 - t55 + t33 + t30;
        
    }// if (flag == 13)
    
    if (flag == 21)
    {
        complex<double> t3 = bb[0] * dd[0];
        complex<double> t7 = bb[1] * dd[1];
        complex<double> t11 = bb[2] * dd[2];
        complex<double> t13 = cc[1] * dd[1];
        complex<double> t19 = cc[0] * dd[0];
        complex<double> t23 = cc[2] * dd[2];
        complex<double> t25 = cc[0] * ee[0] / 0.3e1 + t3 / 0.12e2 - bb[0] * ee[0] / 0.6e1 + t7 / 0.12e2 - bb[2] * ee[2] / 0.6e1 + t11 / 0.12e2 - t13 / 0.6e1 - bb[1] * ee[1] / 0.6e1 + cc[2] * ee[2] / 0.3e1 - t19 / 0.6e1 + cc[1] * ee[1] / 0.3e1 - t23 / 0.6e1;
        complex<double> t26 = -t3 - t7 - t11;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t28 = t27 * bb[1];
        complex<double> t29 = t28 * dd[1];
        complex<double> t30 = t29 / 0.12e2;
        complex<double> t31 = t27 * bb[0];
        complex<double> t32 = t31 * ee[0];
        complex<double> t34 = t31 * dd[0];
        complex<double> t35 = t34 / 0.12e2;
        complex<double> t36 = t27 * bb[2];
        complex<double> t37 = t36 * ee[2];
        complex<double> t39 = t28 * ee[1];
        complex<double> t41 = t36 * dd[2];
        complex<double> t42 = t41 / 0.12e2;
        //
        coef[0] = t25;
        coef[1] = t26 / 0.2e1;
        coef[2] = -t26 / 0.4e1;
        coef[3] = -t30 + t32 / 0.6e1 - t35 + t37 / 0.6e1 + t39 / 0.6e1 - t42;
        coef[4] = t23 * t27 / 0.6e1 - t30 + t19 * t27 / 0.6e1 - t42 + t13 * t27 / 0.6e1 - t35;
        coef[5] = t29 / 0.6e1 + t34 / 0.6e1 - t39 / 0.3e1 + t41 / 0.6e1 - t32 / 0.3e1 - t37 / 0.3e1;
        
    }// if (flag == 21)
    
    if (flag == 22)
    {
        complex<double> t1 = bb[2] * dd[2];
        complex<double> t2 = bb[0] * dd[0];
        complex<double> t3 = bb[1] * dd[1];
        complex<double> t4 = t1 + t2 + t3;
        complex<double> t26 = -bb[1] * ee[1] / 0.6e1 + t3 / 0.12e2 + t1 / 0.12e2 - bb[2] * ee[2] / 0.6e1 + t2 / 0.12e2 - bb[0] * ee[0] / 0.6e1 - cc[0] * dd[0] / 0.6e1 - cc[1] * dd[1] / 0.6e1 - cc[2] * dd[2] / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t28 = t27 * bb[0];
        complex<double> t29 = t28 * dd[0];
        complex<double> t30 = t29 / 0.6e1;
        complex<double> t31 = t27 * bb[1];
        complex<double> t32 = t31 * dd[1];
        complex<double> t33 = t32 / 0.6e1;
        complex<double> t35 = t27 * cc[0] * dd[0];
        complex<double> t38 = t27 * cc[1] * dd[1];
        complex<double> t40 = t27 * bb[2];
        complex<double> t41 = t40 * dd[2];
        complex<double> t42 = t41 / 0.6e1;
        complex<double> t44 = t27 * cc[2] * dd[2];
        complex<double> t47 = t40 * ee[2];
        complex<double> t49 = t28 * ee[0];
        complex<double> t51 = t31 * ee[1];
        complex<double> t55 = t41 / 0.12e2;
        complex<double> t57 = t29 / 0.12e2;
        complex<double> t58 = t32 / 0.12e2;
        //
        coef[0] = t4;
        coef[1] = t26;
        coef[2] = -t4 / 0.2e1;
        coef[3] = t30 + t33 - t35 / 0.3e1 - t38 / 0.3e1 + t42 - t44 / 0.3e1;
        coef[4] = -t4 / 0.2e1;
        coef[5] = -t47 / 0.3e1 + t30 + t33 + t42 - t49 / 0.3e1 - t51 / 0.3e1;
        coef[6] = t51 / 0.6e1 - t55 + t47 / 0.6e1 - t57 - t58 + t49 / 0.6e1;
        coef[7] = t35 / 0.6e1 + t44 / 0.6e1 - t55 - t57 - t58 + t38 / 0.6e1;
        coef[8] = t4 / 0.4e1;
        
    }// if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> t1 = bb[2] * ee[2];
        complex<double> t2 = bb[0] * ee[0];
        complex<double> t3 = bb[1] * ee[1];
        complex<double> t4 = t1 + t2 + t3;
        complex<double> t6 = bb[2] * dd[2];
        complex<double> t8 = bb[1] * dd[1];
        complex<double> t11 = bb[0] * dd[0];
        complex<double> t14 = cc[0] * dd[0];
        complex<double> t16 = cc[1] * dd[1];
        complex<double> t18 = cc[2] * dd[2];
        complex<double> t26 = -t1 / 0.6e1 + t6 / 0.12e2 + t8 / 0.12e2 - t2 / 0.6e1 + t11 / 0.12e2 - t3 / 0.6e1 - t14 / 0.6e1 - t16 / 0.6e1 - t18 / 0.6e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1 + cc[2] * ee[2] / 0.3e1;
        complex<double> t27 = t8 + t11 + t6;
        complex<double> t28 = sqrt(0.3e1);
        complex<double> t31 = t28 * bb[2];
        complex<double> t32 = t31 * dd[2];
        complex<double> t33 = t32 / 0.12e2;
        complex<double> t34 = t28 * bb[1];
        complex<double> t35 = t34 * dd[1];
        complex<double> t36 = t35 / 0.12e2;
        complex<double> t41 = t28 * bb[0];
        complex<double> t42 = t41 * dd[0];
        complex<double> t43 = t42 / 0.12e2;
        complex<double> t45 = t31 * ee[2];
        complex<double> t46 = t45 / 0.6e1;
        complex<double> t47 = t34 * ee[1];
        complex<double> t48 = t47 / 0.6e1;
        complex<double> t49 = t41 * ee[0];
        complex<double> t50 = t49 / 0.6e1;
        //
        coef[0] = t4;
        coef[1] = t26;
        coef[2] = t27 / 0.4e1;
        coef[3] = t18 * t28 / 0.6e1 - t33 - t36 + t14 * t28 / 0.6e1 + t16 * t28 / 0.6e1 - t43;
        coef[4] = -t43 + t46 + t48 - t36 - t33 + t50;
        coef[5] = t46 + t50 - t28 * cc[0] * ee[0] / 0.3e1 - t28 * cc[1] * ee[1] / 0.3e1 + t48 - t28 * cc[2] * ee[2] / 0.3e1;
        coef[6] = -t27 / 0.2e1;
        coef[7] = -t47 / 0.3e1 + t32 / 0.6e1 + t35 / 0.6e1 + t42 / 0.6e1 - t45 / 0.3e1 - t49 / 0.3e1;
        coef[8] = -t4 / 0.2e1;
        
    }// if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t1 = bb[0] * dd[0];
        complex<double> t3 = bb[1] * dd[1];
        complex<double> t5 = bb[2] * dd[2];
        complex<double> t7 = bb[0] * ee[0];
        complex<double> t9 = bb[1] * ee[1];
        complex<double> t11 = bb[2] * ee[2];
        complex<double> t13 = cc[0] * dd[0];
        complex<double> t15 = cc[1] * dd[1];
        complex<double> t17 = cc[2] * dd[2];
        complex<double> t25 = t1 / 0.12e2 + t3 / 0.12e2 + t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 - t11 / 0.6e1 - t13 / 0.6e1 - t15 / 0.6e1 - t17 / 0.6e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1 + cc[1] * ee[1] / 0.3e1;
        complex<double> t26 = sqrt(0.3e1);
        complex<double> t27 = t26 * cc[1];
        complex<double> t29 = t27 * dd[1] / 0.6e1;
        complex<double> t30 = t26 * cc[0];
        complex<double> t32 = t30 * dd[0] / 0.6e1;
        complex<double> t35 = t26 * cc[2];
        complex<double> t37 = t35 * dd[2] / 0.6e1;
        complex<double> t46 = t1 * t26 / 0.12e2;
        complex<double> t48 = t3 * t26 / 0.12e2;
        complex<double> t52 = t5 * t26 / 0.12e2;
        //
        coef[0] = t25;
        coef[1] = t29 + t32 - t27 * ee[1] / 0.3e1 + t37 - t30 * ee[0] / 0.3e1 - t35 * ee[2] / 0.3e1;
        coef[2] = -t17 / 0.2e1 - t13 / 0.2e1 - t15 / 0.2e1;
        coef[3] = t1 / 0.4e1 + t3 / 0.4e1 + t5 / 0.4e1;
        coef[4] = -t46 - t48 + t9 * t26 / 0.6e1 - t52 + t7 * t26 / 0.6e1 + t11 * t26 / 0.6e1;
        coef[5] = -t46 - t52 - t48 + t32 + t29 + t37;
        
    }// if (flag ==31)
    
    if (flag == 32)
    {
        complex<double> t1 = cc[2] * dd[2];
        complex<double> t2 = cc[0] * dd[0];
        complex<double> t3 = cc[1] * dd[1];
        complex<double> t4 = t1 + t2 + t3;
        complex<double> t5 = bb[2] * dd[2];
        complex<double> t7 = bb[1] * ee[1];
        complex<double> t9 = bb[2] * ee[2];
        complex<double> t11 = bb[0] * dd[0];
        complex<double> t15 = bb[1] * dd[1];
        complex<double> t19 = bb[0] * ee[0];
        complex<double> t26 = t5 / 0.12e2 - t7 / 0.6e1 - t9 / 0.6e1 + t11 / 0.12e2 + cc[1] * ee[1] / 0.3e1 + t15 / 0.12e2 - t2 / 0.6e1 - t3 / 0.6e1 - t19 / 0.6e1 - t1 / 0.6e1 + cc[2] * ee[2] / 0.3e1 + cc[0] * ee[0] / 0.3e1;
        complex<double> t27 = t15 + t11 + t5;
        complex<double> t28 = sqrt(0.3e1);
        complex<double> t30 = t28 * bb[1] * dd[1];
        complex<double> t31 = t30 / 0.12e2;
        complex<double> t35 = t28 * bb[2] * dd[2];
        complex<double> t36 = t35 / 0.12e2;
        complex<double> t42 = t28 * bb[0] * dd[0];
        complex<double> t43 = t42 / 0.12e2;
        complex<double> t45 = t28 * cc[1];
        complex<double> t46 = t45 * dd[1];
        complex<double> t47 = t46 / 0.6e1;
        complex<double> t48 = t28 * cc[0];
        complex<double> t49 = t48 * dd[0];
        complex<double> t50 = t49 / 0.6e1;
        complex<double> t51 = t28 * cc[2];
        complex<double> t52 = t51 * dd[2];
        complex<double> t53 = t52 / 0.6e1;
        //
        coef[0] = t4;
        coef[1] = t26;
        coef[2] = -t4 / 0.2e1;
        coef[3] = t27 / 0.4e1;
        coef[4] = -t31 + t9 * t28 / 0.6e1 - t36 + t7 * t28 / 0.6e1 + t19 * t28 / 0.6e1 - t43;
        coef[5] = -t43 + t47 + t50 - t36 + t53 - t31;
        coef[6] = -t48 * ee[0] / 0.3e1 + t47 + t50 + t53 - t51 * ee[2] / 0.3e1 - t45 * ee[1] / 0.3e1;
        coef[7] = -t27 / 0.2e1;
        coef[8] = t35 / 0.6e1 - t49 / 0.3e1 + t30 / 0.6e1 + t42 / 0.6e1 - t52 / 0.3e1 - t46 / 0.3e1;
        
    }// if (flag == 32)
    
    if (flag == 33)
    {
        complex<double> t1 = cc[2] * ee[2];
        complex<double> t2 = cc[0] * ee[0];
        complex<double> t3 = cc[1] * ee[1];
        complex<double> t5 = bb[1] * ee[1];
        complex<double> t7 = bb[2] * dd[2];
        complex<double> t9 = bb[1] * dd[1];
        complex<double> t11 = bb[2] * ee[2];
        complex<double> t13 = bb[0] * dd[0];
        complex<double> t15 = bb[0] * ee[0];
        complex<double> t17 = cc[0] * dd[0];
        complex<double> t19 = cc[1] * dd[1];
        complex<double> t21 = cc[2] * dd[2];
        complex<double> t26 = -t5 / 0.6e1 + t7 / 0.12e2 + t9 / 0.12e2 - t11 / 0.6e1 + t13 / 0.12e2 - t15 / 0.6e1 - t17 / 0.6e1 - t19 / 0.6e1 - t21 / 0.6e1 + t1 / 0.3e1 + t2 / 0.3e1 + t3 / 0.3e1;
        complex<double> t27 = sqrt(0.3e1);
        complex<double> t28 = t27 * cc[2];
        complex<double> t30 = t28 * ee[2] / 0.3e1;
        complex<double> t32 = t28 * dd[2] / 0.6e1;
        complex<double> t33 = t27 * cc[0];
        complex<double> t35 = t33 * dd[0] / 0.6e1;
        complex<double> t36 = t27 * cc[1];
        complex<double> t38 = t36 * dd[1] / 0.6e1;
        complex<double> t40 = t33 * ee[0] / 0.3e1;
        complex<double> t42 = t36 * ee[1] / 0.3e1;
        complex<double> t46 = t27 * bb[1] * ee[1] / 0.6e1;
        complex<double> t49 = t27 * bb[2] * ee[2] / 0.6e1;
        complex<double> t52 = t27 * bb[0] * ee[0] / 0.6e1;
        complex<double> t56 = t9 * t27 / 0.12e2;
        complex<double> t58 = t7 * t27 / 0.12e2;
        complex<double> t60 = t13 * t27 / 0.12e2;
        //
        coef[0] = t1 + t2 + t3;
        coef[1] = t26;
        coef[2] = -t30 + t32 + t35 + t38 - t40 - t42;
        coef[3] = t46 + t49 - t42 - t30 + t52 - t40;
        coef[4] = t9 / 0.4e1 + t13 / 0.4e1 + t7 / 0.4e1;
        coef[5] = -t56 + t49 - t58 - t60 + t46 + t52;
        coef[6] = -t56 - t60 + t38 - t58 + t35 + t32;
        coef[7] = -t5 / 0.2e1 - t11 / 0.2e1 - t15 / 0.2e1;
        coef[8] = -t21 / 0.2e1 - t17 / 0.2e1 - t19 / 0.2e1;
        
    }// if (flag == 33)
    
        
	
}
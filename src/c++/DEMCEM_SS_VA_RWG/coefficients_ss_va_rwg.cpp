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

using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_va
// ***********************************************************************

void coefficients_ss_va_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const complex<double> ko, complex<double> coef[], int flag )
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

    // cases
    if (flag == 12)
    {
        complex<double> t1 = bb[0] * ee[1];
        complex<double> t4 = bb[0] * ee[2];
        complex<double> t7 = bb[1] * ee[2];
        complex<double> t10 = bb[1] * ee[0];
        complex<double> t13 = bb[2] * ee[0];
        complex<double> t16 = bb[2] * ee[1];
        complex<double> t37 = -t1 * dd[2] / 0.6e1 + t4 * dd[1] / 0.6e1 - t7 * dd[0] / 0.6e1 + t10 * dd[2] / 0.6e1 - t13 * dd[1] / 0.6e1 + t16 * dd[0] / 0.6e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[1] * ee[2] * dd[0] / 0.3e1 - cc[1] * ee[0] * dd[2] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 - cc[2] * ee[1] * dd[0] / 0.3e1;
        complex<double> t38 = sqrt(0.3e1);
        complex<double> t39 = dd[2] * t38;
        complex<double> t41 = dd[1] * t38;
        complex<double> t43 = dd[0] * t38;
        //
        coef[0] = t37;
        coef[1] = t1 * t39 / 0.6e1 - t4 * t41 / 0.6e1 + t7 * t43 / 0.6e1 - t10 * t39 / 0.6e1 + t13 * t41 / 0.6e1 - t16 * t43 / 0.6e1;
        
    }// if (flag == 12)
    
    if (flag == 13)
    {
        complex<double> t1 = bb[0] * ee[1];
        complex<double> t2 = sqrt(0.3e1);
        complex<double> t3 = dd[2] * t2;
        complex<double> t4 = t1 * t3;
        complex<double> t5 = bb[0] * dd[1];
        complex<double> t6 = ee[2] * t2;
        complex<double> t7 = t5 * t6;
        complex<double> t8 = bb[1] * ee[2];
        complex<double> t9 = dd[0] * t2;
        complex<double> t10 = t8 * t9;
        complex<double> t11 = bb[1] * dd[2];
        complex<double> t12 = ee[0] * t2;
        complex<double> t13 = t11 * t12;
        complex<double> t14 = bb[2] * ee[0];
        complex<double> t15 = dd[1] * t2;
        complex<double> t16 = t14 * t15;
        complex<double> t17 = bb[2] * dd[0];
        complex<double> t18 = ee[1] * t2;
        complex<double> t19 = t17 * t18;
        complex<double> t21 = cc[1] * dd[2];
        complex<double> t24 = cc[0] * ee[1];
        complex<double> t27 = t11 * ee[0];
        complex<double> t29 = t14 * dd[1];
        complex<double> t31 = t17 * ee[1];
        complex<double> t33 = t5 * ee[2];
        complex<double> t35 = t8 * dd[0];
        complex<double> t37 = t1 * dd[2];
        complex<double> t39 = cc[1] * ee[2];
        complex<double> t42 = cc[2] * ee[0];
        complex<double> t45 = cc[2] * dd[0];
        complex<double> t48 = cc[0] * dd[1];
        complex<double> t51 = -t21 * ee[0] / 0.6e1 + t24 * dd[2] / 0.6e1 + t27 / 0.12e2 - t29 / 0.12e2 + t31 / 0.12e2 + t33 / 0.12e2 - t35 / 0.12e2 - t37 / 0.12e2 + t39 * dd[0] / 0.6e1 + t42 * dd[1] / 0.6e1 - t45 * ee[1] / 0.6e1 - t48 * ee[2] / 0.6e1;
        complex<double> t71 = t4 / 0.12e2 - t7 / 0.12e2 + t10 / 0.12e2 - t13 / 0.12e2 - t24 * t3 / 0.6e1 + t16 / 0.12e2 - t39 * t9 / 0.6e1 - t19 / 0.12e2 - t42 * t15 / 0.6e1 + t48 * t6 / 0.6e1 + t45 * t18 / 0.6e1 + t21 * t12 / 0.6e1;
        //
        coef[0] = t4 / 0.12e2 - t7 / 0.12e2 + t10 / 0.12e2 - t13 / 0.12e2 + t16 / 0.12e2 - t19 / 0.12e2;
        coef[1] = t51;
        coef[2] = t33 / 0.4e1 + t31 / 0.4e1 - t37 / 0.4e1 - t35 / 0.4e1 - t29 / 0.4e1 + t27 / 0.4e1;
        coef[3] = t71;
        
    }// if (flag == 13)
    
    if (flag == 21)
    {
        complex<double> t1 = bb[0] * cc[2];
        complex<double> t4 = bb[0] * cc[1];
        complex<double> t7 = bb[1] * cc[0];
        complex<double> t10 = bb[1] * cc[2];
        complex<double> t13 = bb[2] * cc[1];
        complex<double> t16 = bb[2] * cc[0];
        complex<double> t31 = -t1 * dd[1] / 0.6e1 + t4 * dd[2] / 0.6e1 - t7 * dd[2] / 0.6e1 + t10 * dd[0] / 0.6e1 - t13 * dd[0] / 0.6e1 + t16 * dd[1] / 0.6e1 - t4 * ee[2] / 0.3e1 + t1 * ee[1] / 0.3e1 - t10 * ee[0] / 0.3e1 + t7 * ee[2] / 0.3e1 - t16 * ee[1] / 0.3e1 + t13 * ee[0] / 0.3e1;
        complex<double> t32 = sqrt(0.3e1);
        complex<double> t33 = dd[1] * t32;
        complex<double> t35 = dd[2] * t32;
        complex<double> t38 = dd[0] * t32;
        //
        coef[0] = t31;
        coef[1] = t1 * t33 / 0.6e1 - t4 * t35 / 0.6e1 + t7 * t35 / 0.6e1 - t10 * t38 / 0.6e1 + t13 * t38 / 0.6e1 - t16 * t33 / 0.6e1;
        
    }// if (flag == 21)
    
    if (flag == 22)
    {
        complex<double> t7 = bb[0] * cc[2];
        complex<double> t10 = bb[1] * cc[2];
        complex<double> t19 = bb[1] * cc[0];
        complex<double> t25 = bb[0] * cc[1];
        complex<double> t31 = bb[2] * cc[0];
        complex<double> t34 = bb[2] * cc[1];
        complex<double> t37 = bb[1] * ee[0] * dd[2] / 0.6e1 - bb[2] * ee[0] * dd[1] / 0.6e1 - t7 * dd[1] / 0.6e1 + t10 * dd[0] / 0.6e1 - bb[1] * ee[2] * dd[0] / 0.6e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - t19 * dd[2] / 0.6e1 + cc[1] * ee[2] * dd[0] / 0.3e1 + t25 * dd[2] / 0.6e1 - bb[0] * ee[1] * dd[2] / 0.6e1 + t31 * dd[1] / 0.6e1 - t34 * dd[0] / 0.6e1;
        complex<double> t68 = bb[0] * ee[2] * dd[1] / 0.6e1 + bb[2] * ee[1] * dd[0] / 0.6e1 - t25 * ee[2] / 0.3e1 - t31 * ee[1] / 0.3e1 + t34 * ee[0] / 0.3e1 - cc[1] * ee[0] * dd[2] / 0.3e1 - t10 * ee[0] / 0.3e1 + t19 * ee[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 - cc[2] * ee[1] * dd[0] / 0.3e1 + t7 * ee[1] / 0.3e1;
        complex<double> t70 = sqrt(0.3e1);
        complex<double> t71 = t70 * bb[0];
        complex<double> t76 = t70 * bb[1];
        complex<double> t81 = t70 * bb[2];
        complex<double> t86 = t71 * cc[1] * dd[2] - t71 * cc[2] * dd[1] + t76 * cc[2] * dd[0] - t76 * cc[0] * dd[2] + t81 * cc[0] * dd[1] - t81 * cc[1] * dd[0];
        complex<double> t99 = -t71 * ee[1] * dd[2] + t71 * ee[2] * dd[1] - t76 * ee[2] * dd[0] + t76 * ee[0] * dd[2] - t81 * ee[0] * dd[1] + t81 * ee[1] * dd[0];
        //
        coef[0] = t37 + t68;
        coef[1] = t86 / 0.3e1;
        coef[2] = t99 / 0.3e1;
        coef[3] = -t99 / 0.6e1;
        coef[4] = -t86 / 0.6e1;
        
    }// if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> t2 = bb[2] * dd[0] * ee[1];
        complex<double> t5 = bb[2] * ee[0] * dd[1];
        complex<double> t8 = bb[0] * ee[1] * dd[2];
        complex<double> t10 = bb[2] * cc[0];
        complex<double> t13 = bb[0] * cc[1];
        complex<double> t16 = cc[2] * ee[0];
        complex<double> t19 = cc[0] * ee[1];
        complex<double> t22 = cc[1] * ee[2];
        complex<double> t26 = bb[1] * ee[2] * dd[0];
        complex<double> t29 = bb[1] * dd[2] * ee[0];
        complex<double> t32 = bb[0] * dd[1] * ee[2];
        complex<double> t34 = bb[1] * cc[0];
        complex<double> t37 = t2 / 0.12e2 - t5 / 0.12e2 - t8 / 0.12e2 + t10 * dd[1] / 0.6e1 + t13 * dd[2] / 0.6e1 + t16 * dd[1] / 0.6e1 + t19 * dd[2] / 0.6e1 + t22 * dd[0] / 0.6e1 - t26 / 0.12e2 + t29 / 0.12e2 + t32 / 0.12e2 - t34 * dd[2] / 0.6e1;
        complex<double> t38 = bb[1] * cc[2];
        complex<double> t41 = cc[1] * dd[2];
        complex<double> t44 = cc[0] * dd[1];
        complex<double> t47 = bb[2] * cc[1];
        complex<double> t50 = bb[0] * cc[2];
        complex<double> t55 = cc[2] * dd[0];
        complex<double> t68 = t38 * dd[0] / 0.6e1 - t41 * ee[0] / 0.6e1 - t44 * ee[2] / 0.6e1 - t47 * dd[0] / 0.6e1 - t50 * dd[1] / 0.6e1 + t47 * ee[0] / 0.3e1 - t55 * ee[1] / 0.6e1 - t13 * ee[2] / 0.3e1 + t50 * ee[1] / 0.3e1 - t10 * ee[1] / 0.3e1 + t34 * ee[2] / 0.3e1 - t38 * ee[0] / 0.3e1;
        complex<double> t70 = sqrt(0.3e1);
        complex<double> t71 = t70 * bb[0];
        complex<double> t73 = t71 * ee[1] * dd[2];
        complex<double> t75 = t71 * dd[1] * ee[2];
        complex<double> t76 = t70 * bb[1];
        complex<double> t78 = t76 * ee[2] * dd[0];
        complex<double> t80 = t76 * dd[2] * ee[0];
        complex<double> t81 = t70 * bb[2];
        complex<double> t83 = t81 * ee[0] * dd[1];
        complex<double> t85 = t81 * dd[0] * ee[1];
        complex<double> t86 = -t73 + t75 - t78 + t80 - t83 + t85;
        complex<double> t97 = t8 - t32 + t26 - t29 + t5 - t2;
        complex<double> t98 = dd[2] * t70;
        complex<double> t108 = dd[1] * t70;
        complex<double> t116 = dd[0] * t70;
        complex<double> t134 = -t13 * t98 / 0.6e1 - t19 * t98 / 0.6e1 + t73 / 0.12e2 + t34 * t98 / 0.6e1 - t75 / 0.12e2 + t78 / 0.12e2 - t10 * t108 / 0.6e1 + t41 * ee[0] * t70 / 0.6e1 - t16 * t108 / 0.6e1 - t38 * t116 / 0.6e1 - t22 * t116 / 0.6e1 + t44 * ee[2] * t70 / 0.6e1 + t83 / 0.12e2 + t50 * t108 / 0.6e1 - t85 / 0.12e2 + t55 * ee[1] * t70 / 0.6e1 - t80 / 0.12e2 + t47 * t116 / 0.6e1;
        //
        coef[0] = t37 + t68;
        coef[1] = t86 / 0.6e1;
        coef[2] = t71 * t22 / 0.3e1 - t71 * cc[2] * ee[1] / 0.3e1 + t76 * t16 / 0.3e1 - t76 * cc[0] * ee[2] / 0.3e1 + t81 * t19 / 0.3e1 - t81 * cc[1] * ee[0] / 0.3e1;
        coef[3] = t97 / 0.2e1;
        coef[4] = -t97 / 0.4e1;
        coef[5] = -t86 / 0.12e2;
        coef[6] = t134;
        
    }// if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t1 = cc[0] * bb[1];
        complex<double> t2 = sqrt(0.3e1);
        complex<double> t3 = dd[2] * t2;
        complex<double> t4 = t1 * t3;
        complex<double> t6 = cc[0] * bb[2];
        complex<double> t7 = ee[1] * t2;
        complex<double> t10 = cc[2] * bb[1];
        complex<double> t11 = dd[0] * t2;
        complex<double> t12 = t10 * t11;
        complex<double> t14 = dd[1] * t2;
        complex<double> t15 = t6 * t14;
        complex<double> t17 = ee[2] * t2;
        complex<double> t20 = cc[1] * bb[2];
        complex<double> t21 = t20 * t11;
        complex<double> t23 = ee[0] * t2;
        complex<double> t26 = cc[1] * bb[0];
        complex<double> t27 = t26 * t3;
        complex<double> t29 = cc[2] * bb[0];
        complex<double> t32 = t29 * t14;
        complex<double> t38 = t4 / 0.12e2 + t6 * t7 / 0.6e1 - t12 / 0.12e2 - t15 / 0.12e2 - t1 * t17 / 0.6e1 + t21 / 0.12e2 - t20 * t23 / 0.6e1 - t27 / 0.12e2 - t29 * t7 / 0.6e1 + t32 / 0.12e2 + t10 * t23 / 0.6e1 + t26 * t17 / 0.6e1;
        complex<double> t39 = t29 * dd[1];
        complex<double> t41 = t1 * dd[2];
        complex<double> t45 = t6 * dd[1];
        complex<double> t49 = t20 * dd[0];
        complex<double> t57 = t10 * dd[0];
        complex<double> t59 = t26 * dd[2];
        complex<double> t63 = -t39 / 0.12e2 - t41 / 0.12e2 + t1 * ee[2] / 0.6e1 + t45 / 0.12e2 - t6 * ee[1] / 0.6e1 - t49 / 0.12e2 - t10 * ee[0] / 0.6e1 + t29 * ee[1] / 0.6e1 + t20 * ee[0] / 0.6e1 + t57 / 0.12e2 + t59 / 0.12e2 - t26 * ee[2] / 0.6e1;
        //
        coef[0] = t38;
        coef[1] = t63;
        coef[2] = -t39 / 0.4e1 + t57 / 0.4e1 + t45 / 0.4e1 + t59 / 0.4e1 - t49 / 0.4e1 - t41 / 0.4e1;
        coef[3] = -t27 / 0.12e2 + t21 / 0.12e2 + t4 / 0.12e2 - t12 / 0.12e2 - t15 / 0.12e2 + t32 / 0.12e2;
        
    }// if (flag ==31)
    
    if (flag == 32)
    {
        complex<double> t4 = cc[0] * bb[1];
        complex<double> t5 = t4 * dd[2];
        complex<double> t7 = cc[0] * bb[2];
        complex<double> t8 = t7 * dd[1];
        complex<double> t10 = cc[2] * bb[1];
        complex<double> t11 = t10 * dd[0];
        complex<double> t13 = bb[0] * ee[1];
        complex<double> t16 = cc[1] * bb[0];
        complex<double> t17 = t16 * dd[2];
        complex<double> t21 = bb[1] * ee[0];
        complex<double> t26 = bb[2] * ee[1];
        complex<double> t29 = bb[0] * ee[2];
        complex<double> t32 = bb[2] * ee[0];
        complex<double> t35 = -cc[2] * ee[1] * dd[0] / 0.3e1 - t5 / 0.12e2 + t8 / 0.12e2 + t11 / 0.12e2 - t13 * dd[2] / 0.6e1 + t17 / 0.12e2 + t4 * ee[2] / 0.6e1 + t21 * dd[2] / 0.6e1 - t16 * ee[2] / 0.6e1 + t26 * dd[0] / 0.6e1 + t29 * dd[1] / 0.6e1 - t32 * dd[1] / 0.6e1;
        complex<double> t36 = cc[1] * bb[2];
        complex<double> t37 = t36 * dd[0];
        complex<double> t43 = cc[2] * bb[0];
        complex<double> t44 = t43 * dd[1];
        complex<double> t65 = bb[1] * ee[2];
        complex<double> t68 = -t37 / 0.12e2 + t36 * ee[0] / 0.6e1 - t7 * ee[1] / 0.6e1 - t44 / 0.12e2 - t10 * ee[0] / 0.6e1 - cc[1] * ee[0] * dd[2] / 0.3e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 + cc[1] * ee[2] * dd[0] / 0.3e1 + t43 * ee[1] / 0.6e1 - t65 * dd[0] / 0.6e1;
        complex<double> t70 = sqrt(0.3e1);
        complex<double> t71 = t70 * cc[0];
        complex<double> t73 = t71 * bb[1] * dd[2];
        complex<double> t75 = t71 * bb[2] * dd[1];
        complex<double> t76 = t70 * cc[1];
        complex<double> t78 = t76 * bb[2] * dd[0];
        complex<double> t80 = t76 * bb[0] * dd[2];
        complex<double> t81 = t70 * cc[2];
        complex<double> t83 = t81 * bb[0] * dd[1];
        complex<double> t85 = t81 * bb[1] * dd[0];
        complex<double> t86 = -t73 + t75 - t78 + t80 - t83 + t85;
        complex<double> t87 = t5 - t8 + t37 - t17 + t44 - t11;
        complex<double> t101 = dd[2] * t70;
        complex<double> t104 = ee[1] * t70;
        complex<double> t110 = dd[0] * t70;
        complex<double> t115 = ee[2] * t70;
        complex<double> t122 = dd[1] * t70;
        complex<double> t125 = ee[0] * t70;
        complex<double> t137 = t13 * t101 / 0.6e1 + t7 * t104 / 0.6e1 - t43 * t104 / 0.6e1 + t73 / 0.12e2 + t65 * t110 / 0.6e1 - t85 / 0.12e2 - t80 / 0.12e2 - t4 * t115 / 0.6e1 - t75 / 0.12e2 - t21 * t101 / 0.6e1 + t83 / 0.12e2 - t29 * t122 / 0.6e1 - t36 * t125 / 0.6e1 + t32 * t122 / 0.6e1 - t26 * t110 / 0.6e1 + t10 * t125 / 0.6e1 + t78 / 0.12e2 + t16 * t115 / 0.6e1;
        //
        coef[0] = t35 + t68;
        coef[1] = t86 / 0.6e1;
        coef[2] = t87 / 0.2e1;
        coef[3] = -t71 * ee[1] * dd[2] / 0.3e1 + t71 * ee[2] * dd[1] / 0.3e1 - t76 * ee[2] * dd[0] / 0.3e1 + t76 * ee[0] * dd[2] / 0.3e1 - t81 * ee[0] * dd[1] / 0.3e1 + t81 * ee[1] * dd[0] / 0.3e1;
        coef[4] = -t87 / 0.4e1;
        coef[5] = t137;
        coef[6] = -t86 / 0.12e2;
        
    }// if (flag == 32)
    
    if (flag == 33)
    {
        complex<double> t1 = sqrt(0.3e1);
        complex<double> t2 = t1 * cc[0];
        complex<double> t3 = bb[1] * ee[2];
        complex<double> t4 = t2 * t3;
        complex<double> t6 = t2 * bb[2] * ee[1];
        complex<double> t7 = t1 * cc[1];
        complex<double> t8 = bb[2] * ee[0];
        complex<double> t9 = t7 * t8;
        complex<double> t11 = t7 * bb[0] * ee[2];
        complex<double> t12 = t1 * cc[2];
        complex<double> t13 = bb[0] * ee[1];
        complex<double> t14 = t12 * t13;
        complex<double> t16 = t12 * bb[1] * ee[0];
        complex<double> t19 = cc[0] * dd[1] * ee[2];
        complex<double> t21 = cc[1] * ee[2] * dd[0];
        complex<double> t23 = cc[0] * ee[1] * dd[2];
        complex<double> t25 = cc[2] * dd[0] * ee[1];
        complex<double> t27 = cc[1] * dd[2] * ee[0];
        complex<double> t29 = cc[2] * ee[0] * dd[1];
        complex<double> t32 = t2 * ee[1] * dd[2];
        complex<double> t34 = t2 * dd[1] * ee[2];
        complex<double> t36 = t7 * ee[2] * dd[0];
        complex<double> t38 = t7 * dd[2] * ee[0];
        complex<double> t40 = t12 * ee[0] * dd[1];
        complex<double> t42 = t12 * dd[0] * ee[1];
        complex<double> t44 = cc[1] * bb[2];
        complex<double> t45 = t44 * ee[0];
        complex<double> t46 = cc[1] * bb[0];
        complex<double> t47 = t46 * ee[2];
        complex<double> t48 = cc[2] * bb[0];
        complex<double> t49 = t48 * ee[1];
        complex<double> t50 = cc[2] * bb[1];
        complex<double> t51 = t50 * ee[0];
        complex<double> t52 = cc[0] * bb[2];
        complex<double> t53 = t52 * ee[1];
        complex<double> t54 = cc[0] * bb[1];
        complex<double> t55 = t54 * ee[2];
        complex<double> t57 = dd[2] * t1;
        complex<double> t59 = t54 * t57 / 0.12e2;
        complex<double> t60 = bb[2] * dd[0];
        complex<double> t63 = t60 * ee[1] * t1 / 0.12e2;
        complex<double> t64 = dd[1] * t1;
        complex<double> t66 = t8 * t64 / 0.12e2;
        complex<double> t68 = t48 * t64 / 0.12e2;
        complex<double> t69 = dd[0] * t1;
        complex<double> t71 = t44 * t69 / 0.12e2;
        complex<double> t73 = t3 * t69 / 0.12e2;
        complex<double> t75 = t52 * t64 / 0.12e2;
        complex<double> t76 = bb[1] * dd[2];
        complex<double> t79 = t76 * ee[0] * t1 / 0.12e2;
        complex<double> t81 = t46 * t57 / 0.12e2;
        complex<double> t83 = t50 * t69 / 0.12e2;
        complex<double> t85 = t13 * t57 / 0.12e2;
        complex<double> t88 = bb[0] * dd[1];
        complex<double> t91 = t88 * ee[2] * t1 / 0.12e2;
        complex<double> t96 = t59 - t63 + t66 + t68 + t71 + t73 - t75 - t79 - t81 - t83 + t85 - t40 / 0.6e1 - t36 / 0.6e1 - t91 + t38 / 0.6e1 - t32 / 0.6e1 + t42 / 0.6e1 + t34 / 0.6e1;
        complex<double> t97 = t54 * dd[2];
        complex<double> t98 = t3 * dd[0];
        complex<double> t99 = t88 * ee[2];
        complex<double> t100 = t13 * dd[2];
        complex<double> t101 = t50 * dd[0];
        complex<double> t102 = t44 * dd[0];
        complex<double> t103 = t8 * dd[1];
        complex<double> t104 = t76 * ee[0];
        complex<double> t105 = t46 * dd[2];
        complex<double> t106 = t52 * dd[1];
        complex<double> t107 = t60 * ee[1];
        complex<double> t108 = t48 * dd[1];
        complex<double> t109 = -t97 - t98 + t99 - t100 + t101 - t102 - t103 + t104 + t105 + t106 + t107 - t108;
        complex<double> t122 = t104 / 0.12e2 - t98 / 0.12e2 - t102 / 0.12e2 + t99 / 0.12e2 + t105 / 0.12e2 - t100 / 0.12e2 - t103 / 0.12e2 - t108 / 0.12e2 + t107 / 0.12e2 - t97 / 0.12e2 + t106 / 0.12e2 + t49 / 0.6e1;
        complex<double> t135 = t101 / 0.12e2 - t27 / 0.6e1 - t47 / 0.6e1 + t23 / 0.6e1 + t55 / 0.6e1 - t25 / 0.6e1 - t51 / 0.6e1 - t19 / 0.6e1 + t45 / 0.6e1 + t29 / 0.6e1 - t53 / 0.6e1 + t21 / 0.6e1;
        complex<double> t143 = t85 - t63 + t66 + t6 / 0.6e1 - t91 - t83 - t4 / 0.6e1 + t71 + t11 / 0.6e1 + t59 - t9 / 0.6e1 - t81 + t16 / 0.6e1 + t68 - t14 / 0.6e1 - t75 - t79 + t73;
        //
        coef[0] = -t4 / 0.6e1 + t6 / 0.6e1 - t9 / 0.6e1 + t11 / 0.6e1 - t14 / 0.6e1 + t16 / 0.6e1;
        coef[1] = -t19 / 0.2e1 + t21 / 0.2e1 + t23 / 0.2e1 - t25 / 0.2e1 - t27 / 0.2e1 + t29 / 0.2e1;
        coef[2] = -t32 / 0.6e1 + t34 / 0.6e1 - t36 / 0.6e1 + t38 / 0.6e1 - t40 / 0.6e1 + t42 / 0.6e1;
        coef[3] = t45 / 0.2e1 - t47 / 0.2e1 + t49 / 0.2e1 - t51 / 0.2e1 - t53 / 0.2e1 + t55 / 0.2e1;
        coef[4] = t96;
        coef[5] = t109 / 0.4e1;
        coef[6] = t122 + t135;
        coef[7] = t143;
        
    }// if (flag == 33)
    
        
	
}
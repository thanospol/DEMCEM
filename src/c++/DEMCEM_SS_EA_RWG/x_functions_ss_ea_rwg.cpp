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
#include "demcem_ss_ea_rwg.h"


// ***********************************************************************
//			IMPLEMENTATION OF complex<double> x_functions
// ***********************************************************************

complex<double> x_functions_ss_ea_rwg (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[], int flag)
{
    complex<double> X;
    //
    
    if (flag == 11)
    {
        double t2  = pow(B,2);
        double t3  = pow(t2,2);
        double t5  = 1.0 / t3 /B;
        double t9  = cos(Psi);
        double t10 = pow(t9,2);
        double t11 = sin(Psi);
        double t13 = sin(theta);
        double t14 = t10 * t11 * t13;
        double t20 = t11 * t13;
        double t39 = pow(Bm,2);
        double t40 = pow(t39,2);
        double t42 = 1.0 / t40 / Bm;
        //
        X = N[0] * t5 * coef[2] * t14 + N[2] / t3 * t10 * t20 * coef[0] + N[3] / t2 / B * t10 * t20 * coef[1]
        + N[1] * t5 * t10 * t20 * coef[3] + Nm[0] * t42 * coefm[0] * t14 + Nm[2] / t40 * t10 * t20 * coefm[2]
        + Nm[3] / t39 / Bm * t10 * t20 * coefm[3] + Nm[1] * t42 * t10 * t20 * coefm[1];
    } // if (flag == 11)
    
    
    if (flag == 13)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = sin(Psi);
        double t9 = t7 * t8;
        double t12 = 0.1e1 / t3 / B;
        double t15 = t7 * t7;
        double t17 = sin(theta);
        double t18 = t15 * t8 * t17;
        double t22 = cos(theta);
        double t24 = t15 * t22 * t8;
        double t40 = 0.1e1 / t2 / B;
        double t79 = Bm * Bm;
        double t80 = t79 * t79;
        double t81 = 0.1e1 / t80;
        double t86 = 0.1e1 / t80 / Bm;
        double t107 = 0.1e1 / t79 / Bm;
        //
        X = N[0] * (t4 * coef[0] * t9 + t12 * coef[9] * t18 + t12 * coef[8] * t24) + N[7] * t4 * t9 * coef[11] + N[5] * t4 * t9 * coef[6]
        + N[3] * (t40 * coef[1] * t18 + t40 * coef[3] * t24) + N[9] * t40 * t9 * coef[2] + N[2] * (t4 * coef[4] * t24 + t4 * coef[5] * t18 + t40 * coef[10] * t9)
        + N[1] * (t12 * coef[12] * t24 + t12 * coef[13] * t18 + t4 * coef[7] * t9) + Nm[0] * (t81 * coefm[0] * t9 + t86 * coefm[2] * t24 + t86 * coefm[3] * t18)
        + Nm[5] * t81 * t9 * coefm[1] + Nm[7] * t81 * t9 * coefm[10] + Nm[3] * (t107 * coefm[8] * t24 + t107 * coefm[9] * t18)
        + Nm[2] * (t81 * coefm[6] * t24 + t81 * coefm[7] * t18 + t107 * coefm[13] * t9) + Nm[1] * (t86 * coefm[11] * t24 + t86 * coefm[12] * t18 + t81 * coefm[4] * t9) + Nm[9] * t107 * t9 * coefm[5];
    } // if (flag == 13)
    
    
    if (flag == 22)
    {
        double t2   = pow(B,2);
        double t3   = pow(t2,2);
        double t5   = 0.1e1 / t3 / B;
        double t9   = cos(Psi);
        double t10  = pow(t9,2);
        double t11  = sin(Psi);
        double t13  = sin(theta);
        double t14  = t10 * t11 * t13;
        double t20  = t10 * t13;
        double t39  = pow(Bm,2);
        double t40  = pow(t39,2);
        double t42  = 0.1e1 / t40 / Bm;
        //
        X = N[0] * t5 * coef[1] * t14 + N[2] / t3 * t11 * t20 * coef[2] + N[3] / t2 / B * t11 * t20 * coef[0] + N[1] * t5 * t11 * t20 * coef[3] + Nm[0] * t42 * coefm[2] * t14 + Nm[2] / t40 * t11 * t20 * coefm[1] + Nm[3] / t39 / Bm * t11 * t20 * coefm[0] + Nm[1] * t42 * t11 * t20 * coefm[3];
        
    } // if (flag == 22)
    
    if (flag == 23)
    {
        double t2   = pow(B,2);
        double t3   = pow(t2,2);
        double t4   = 0.1e1 / t3;
        double t7   = cos(Psi);
        double t8   = sin(Psi);
        double t9   = t7 * t8;
        double t12  = 0.1e1 / t3 / B;
        double t15  = t7 * t7;
        double t17  = sin(theta);
        double t18  = t15 * t8 * t17;
        double t22  = cos(theta);
        double t24  = t15 * t22 * t8;
        double t52  = 0.1e1 / t2 / B;
        double t79  = Bm * Bm;
        double t80  = t79 * t79;
        double t81  = 0.1e1 / t80;
        double t86  = 0.1e1 / t80 / Bm;
        double t107 = 0.1e1 / t79 / Bm;
        //
        X = N[0] * (t4 * coef[0] * t9 + t12 * coef[2] * t18 + t12 * coef[3] * t24) + N[5] * t4 * t9 * coef[1]
        + N[7] * t4 * t9 * coef[11] + N[1] * (t12 * coef[12] * t18 + t12 * coef[13] * t24 + t4 * coef[4] * t9)
        + N[9] * t52 * t9 * coef[5] + N[3] * (t52 * coef[6] * t18 + t52 * coef[7] * t24)
        + N[2] * (t4 * coef[8] * t18 + t4 * coef[9] * t24 + t52 * coef[10] * t9)
        + Nm[0] * (t81 * coefm[8] * t9 + t86 * coefm[1] * t24 + t86 * coefm[2] * t18)
        + Nm[5] * t81 * t9 * coefm[3] + Nm[7] * t81 * t9 * coefm[4] + Nm[3] * (t107 * coefm[11] * t24 + t107 * coefm[12] * t18)
        + Nm[9] * t107 * t9 * coefm[13] + Nm[1] * (t86 * coefm[5] * t18 + t86 * coefm[6] * t24 + t81 * coefm[0] * t9)
        + Nm[2] * (t81 * coefm[9] * t18 + t81 * coefm[10] * t24 + t107 * coefm[7] * t9);
        
    } // if (flag == 23)
    
    if (flag == 31)
    {
        double t2  = pow(B,2);
        double t3  = pow(t2,2);
        double t4  = 1.0 / t3;
        double t7  = cos(Psi);
        double t8  = pow(t7,2);
        double t9  = sin(theta);
        double t10 = t8 * t9;
        double t13 = 1.0 / t3 / B;
        double t16 = sin(Psi);
        double t18 = t8 * t16 * t9;
        double t34 = 1.0 / t2 / B;
        double t64 = pow(Bm,2);
        double t65 = pow(t64,2);
        double t67 = 1.0 / t65 / Bm;
        double t71 = 1.0 / t65;
        double t89 = 1.0 / t64 / Bm;
        //
        X = N[0] * ( t4 * coef[0] * t10 + t13 * coef[4] * t18 )   +   N[7] * t4 * t10 * coef[8]   +   N[5] * t4 * t10 * coef[5]   +   N[9] * t34 * t10 * coef[1]
        +   N[2] * ( t34 * coef[7] * t10 + t4 * coef[3] * t18 )   +   N[3] * t16 * t34 * t10 * coef[2]
        +   N[1] * ( t4 * coef[6] * t10 + t13 * coef[9] * t18 )   +   Nm[0] * ( t67 * coefm[1] * t18 + t71 * coefm[3] * t10 )
        +	 Nm[5] * t71 * t10 * coefm[0]   +   Nm[7] * t71 * t10 * coefm[7]   +    Nm[2] * ( t89 * coefm[9] * t10 + t71 * coefm[6] * t18 )
        +   Nm[3] * t16 * t89 * t10 * coefm[4]   +   Nm[1] * ( t71 * coefm[2] * t10 + t67 * coefm[8] * t18 )   +   Nm[9] * t89 * t10 * coefm[5] ;
        
    } // if (flag == 31)
    
    if (flag == 32)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / B;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = sin(Psi);
        double t12 = sin(theta);
        double t13 = t9 * t10 * t12;
        double t15 = 0.1e1 / t3;
        double t18 = t9 * t12;
        double t29 = 0.1e1 / t2 / B;
        double t64 = Bm * Bm;
        double t65 = t64 * t64;
        double t67 = 0.1e1 / t65 / Bm;
        double t71 = 0.1e1 / t65;
        double t79 = 0.1e1 / t64 / Bm;
        //
        X = N[0] * (t5 * coef[2] * t13 + t15 * coef[3] * t18) + N[7] * t15 * t18 * coef[5] + N[2] * (t29 * coef[4] * t18 + t15 * coef[9] * t13)
        + N[1] * (t15 * coef[1] * t18 + t5 * coef[6] * t13) + N[5] * t15 * t18 * coef[0] + N[9] * t29 * t18 * coef[8]
        + N[3] * t10 * t29 * t18 * coef[7] + Nm[0] * (t67 * coefm[0] * t13 + t71 * coefm[2] * t18)
        + Nm[2] * (t79 * coefm[8] * t18 + t71 * coefm[5] * t13) + Nm[5] * t71 * t18 * coefm[9]
        + Nm[1] * (t71 * coefm[1] * t18 + t67 * coefm[7] * t13) + Nm[9] * t79 * t18 * coefm[3]
        + Nm[7] * t71 * t18 * coefm[6] + Nm[3] * t10 * t79 * t18 * coefm[4];
        
    } // if (flag == 32)
    
    if (flag == 33)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = cos(theta);
        double t10 = t8 * t9;
        double t13 = 0.1e1 / t3 / B;
        double t16 = sin(Psi);
        double t17 = t10 * t16;
        double t24 = sin(theta);
        double t25 = t8 * t24;
        double t29 = t7 * t16;
        double t35 = 0.1e1 / t2 / B;
        double t79 = Bm * Bm;
        double t80 = t79 * t79;
        double t81 = 0.1e1 / t80;
        double t86 = 0.1e1 / t80 / Bm;
        double t112 = 0.1e1 / t79 / Bm;
        //
        X = N[0] * (t4 * coef[0] * t10 + t13 * coef[3] * t17) + N[5] * (t4 * coef[1] * t25 + t4 * coef[2] * t29)
        + N[9] * (t35 * coef[5] * t25 + t35 * coef[6] * t29) + N[3] * t9 * t16 * t35 * t8 * coef[7] + N[2] * (t35 * coef[9] * t10 + t4 * coef[8] * t17)
        + N[1] * (t4 * coef[4] * t10 + t13 * coef[12] * t17) + N[7] * (t4 * coef[10] * t25 + t4 * coef[11] * t29)
        + Nm[0] * (t81 * coefm[4] * t10 + t86 * coefm[7] * t17) + Nm[5] * (t81 * coefm[5] * t25 + t81 * coefm[6] * t29)
        + Nm[7] * (t81 * coefm[1] * t25 + t81 * coefm[2] * t29) + Nm[9] * (t112 * coefm[9] * t25 + t112 * coefm[10] * t29)
        + Nm[3] * t9 * t16 * t112 * t8 * coefm[11] + Nm[2] * (t112 * coefm[0] * t10 + t81 * coefm[12] * t17)
        + Nm[1] * (t81 * coefm[8] * t10 + t86 * coefm[3] * t17);
        
    } // if (flag == 33)
    
    // Final Output
    return X;
	
}



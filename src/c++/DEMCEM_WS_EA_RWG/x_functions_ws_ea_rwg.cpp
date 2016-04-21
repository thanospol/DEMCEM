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
#include "demcem_ws_ea_rwg.h"


// ***********************************************************************
//			IMPLEMENTATION OF complex<double> x_functions
// ***********************************************************************

complex<double> x_functions_ws_ea_rwg (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[], int flag)
{
    complex<double> X;
    //
    if (flag == 0)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t22 = Bm * Bm;
        double t24 = 0.1e1 / t22 / Bm;
        //
        X = N[0] * t4 * coef[2] * t7 + N[1] * t4 * t7 * coef[0] + N[2] / t2 * t7 * coef[1] + Nm[0] * t24 * coefm[1] * t7 + Nm[1] * t24 * t7 * coefm[0] + Nm[2] / t22 * t7 * coefm[2];
    }
    
    if (flag == 11)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = t7 * t7;
        double t14 = cos(theta);
        double t15 = t13 * t14;
        double t18 = 0.1e1 / t9 / B;
        double t21 = sin(Psi);
        double t23 = sin(theta);
        double t24 = t13 * t21 * t23;
        double t28 = t7 * t21;
        double t32 = t13 * t23;
        double t36 = t15 * t21;
        double t58 = 0.1e1 / t2;
        complex<double> t173 = N[0] * (t4 * coef[27] * t7 + t10 * coef[4] * t15 + t18 * coef[19] * t24 + t10 * coef[2] * t28 + t10 * coef[3] * t32 + t18 * coef[18] * t36) + N[5] * (t10 * coef[20] * t32 + t10 * coef[21] * t15 + t10 * coef[17] * t28) + N[6] * t4 * t7 * coef[1] + N[10] * t58 * t7 * coef[29] + N[7] * (t10 * coef[32] * t15 + t10 * coef[38] * t32 + t10 * coef[34] * t28) + N[11] * (t58 * coef[5] * t32 + t58 * coef[6] * t15 + t58 * coef[7] * t28) + N[9] * (t4 * coef[11] * t32 + t4 * coef[12] * t15 + t4 * coef[8] * t28) + N[4] * (t58 * coef[9] * t36 + t58 * coef[10] * t24) + N[2] * (t4 * coef[28] * t32 + t4 * coef[36] * t15 + t10 * coef[13] * t36 + t10 * coef[14] * t24 + t58 * coef[23] * t7 + t4 * coef[35] * t28) + N[3] * (t58 * coef[30] * t15 + t58 * coef[33] * t32 + t4 * coef[15] * t36 + t4 * coef[16] * t24 + t58 * coef[37] * t28) + N[1] * (t10 * coef[24] * t32 + t10 * coef[25] * t15 + t18 * coef[31] * t36 + t18 * coef[39] * t24 + t4 * coef[0] * t7 + t10 * coef[26] * t28) + N[8] * t4 * t7 * coef[22];
        double t175 = Bm * Bm;
        double t176 = t175 * t175;
        double t177 = 0.1e1 / t176;
        double t185 = 0.1e1 / t175 / Bm;
        double t190 = 0.1e1 / t176 / Bm;
        double t232 = 0.1e1 / t175;
        complex<double> t335 = Nm[0] * (t177 * coefm[12] * t28 + t177 * coefm[13] * t15 + t185 * coefm[27] * t7 + t190 * coefm[3] * t24 + t190 * coefm[4] * t36 + t177 * coefm[11] * t32) + Nm[5] * (t177 * coefm[1] * t32 + t177 * coefm[2] * t15 + t177 * coefm[0] * t28) + Nm[8] * t185 * t7 * coefm[8] + Nm[7] * (t177 * coefm[29] * t15 + t177 * coefm[36] * t32 + t177 * coefm[31] * t28) + Nm[10] * t232 * t7 * coefm[32] + Nm[6] * t185 * t7 * coefm[14] + Nm[11] * (t232 * coefm[15] * t32 + t232 * coefm[16] * t15 + t232 * coefm[18] * t28) + Nm[9] * (t185 * coefm[21] * t32 + t185 * coefm[22] * t15 + t185 * coefm[17] * t28) + Nm[4] * (t232 * coefm[19] * t24 + t232 * coefm[20] * t36) + Nm[2] * (t185 * coefm[28] * t15 + t185 * coefm[30] * t32 + t177 * coefm[23] * t24 + t177 * coefm[24] * t36 + t232 * coefm[9] * t7 + t185 * coefm[34] * t28) + Nm[3] * (t232 * coefm[37] * t15 + t232 * coefm[39] * t32 + t185 * coefm[25] * t24 + t185 * coefm[26] * t36 + t232 * coefm[38] * t28) + Nm[1] * (t177 * coefm[5] * t32 + t177 * coefm[6] * t15 + t190 * coefm[33] * t36 + t190 * coefm[35] * t24 + t185 * coefm[10] * t7 + t177 * coefm[7] * t28);
        //
        X = t335 + t173;
    } // if (flag == 11)
    
    if (flag == 12)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = sin(theta);
        double t10 = t8 * t9;
        double t14 = cos(theta);
        double t15 = t8 * t14;
        double t18 = 0.1e1 / t3 / B;
        double t21 = sin(Psi);
        double t22 = t15 * t21;
        double t26 = t7 * t21;
        double t29 = 0.1e1 / t2 / B;
        double t36 = t8 * t21 * t9;
        double t53 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[4] * t10 + t4 * coef[1] * t15 + t18 * coef[5] * t22 + t4 * coef[3] * t26 + t29 * coef[29] * t7 + t18 * coef[6] * t36) + N[2] * (t29 * coef[37] * t10 + t29 * coef[39] * t15 + t4 * coef[24] * t22 + t4 * coef[25] * t36 + t53 * coef[15] * t7 + t29 * coef[30] * t26) + N[5] * (t4 * coef[7] * t15 + t4 * coef[8] * t10 + t29 * coef[2] * t7 + t4 * coef[9] * t26) + N[6] * t29 * t7 * coef[28] + N[11] * (t53 * coef[16] * t15 + t53 * coef[17] * t10 + t53 * coef[18] * t26) + N[4] * (t53 * coef[20] * t22 + t53 * coef[21] * t36) + N[1] * (t4 * coef[10] * t10 + t4 * coef[11] * t15 + t18 * coef[34] * t36 + t18 * coef[35] * t22 + t29 * coef[0] * t7 + t4 * coef[12] * t26) + N[8] * t29 * t7 * coef[13] + N[10] * t53 * t7 * coef[31] + N[7] * (t4 * coef[40] * t10 + t4 * coef[41] * t15 + t29 * coef[14] * t7 + t4 * coef[32] * t26) + N[9] * (t29 * coef[22] * t15 + t29 * coef[23] * t10 + t53 * coef[33] * t7 + t29 * coef[19] * t26) + N[3] * (t53 * coef[38] * t10 + t53 * coef[42] * t15 + t29 * coef[26] * t22 + t29 * coef[27] * t36 + t53 * coef[36] * t26);
        double t184 = Bm * Bm;
        double t186 = 0.1e1 / t184 / Bm;
        double t190 = t184 * t184;
        double t191 = 0.1e1 / t190;
        double t202 = 0.1e1 / t190 / Bm;
        double t227 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t186 * coefm[0] * t7 + t191 * coefm[27] * t15 + t191 * coefm[29] * t10 + t191 * coefm[10] * t26 + t202 * coefm[9] * t36 + t202 * coefm[11] * t22) + Nm[5] * (t191 * coefm[3] * t15 + t191 * coefm[4] * t10 + t186 * coefm[14] * t7 + t191 * coefm[12] * t26) + Nm[11] * (t227 * coefm[15] * t15 + t227 * coefm[16] * t10 + t227 * coefm[18] * t26) + Nm[6] * t186 * t7 * coefm[13] + Nm[2] * (t186 * coefm[32] * t15 + t186 * coefm[39] * t10 + t191 * coefm[22] * t36 + t191 * coefm[23] * t22 + t227 * coefm[8] * t7 + t186 * coefm[34] * t26) + Nm[9] * (t186 * coefm[21] * t10 + t186 * coefm[26] * t15 + t227 * coefm[31] * t7 + t186 * coefm[17] * t26) + Nm[7] * (t191 * coefm[33] * t10 + t191 * coefm[36] * t15 + t186 * coefm[6] * t7 + t191 * coefm[42] * t26) + Nm[10] * t227 * t7 * coefm[38] + Nm[8] * t186 * t7 * coefm[5] + Nm[4] * (t227 * coefm[19] * t36 + t227 * coefm[20] * t22) + Nm[1] * (t191 * coefm[1] * t10 + t191 * coefm[2] * t15 + t202 * coefm[35] * t36 + t202 * coefm[41] * t22 + t186 * coefm[28] * t7 + t191 * coefm[7] * t26) + Nm[3] * (t227 * coefm[30] * t15 + t227 * coefm[37] * t10 + t186 * coefm[24] * t36 + t186 * coefm[25] * t22 + t227 * coefm[40] * t26);
        //
        X = t353 + t182;
        
    } // if (flag == 12)
    
    if (flag == 13)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = cos(theta);
        double t10 = t8 * t9;
        double t14 = sin(Psi);
        double t15 = t7 * t14;
        double t18 = 0.1e1 / t3 / B;
        double t21 = t10 * t14;
        double t26 = sin(theta);
        double t27 = t8 * t14 * t26;
        double t31 = t8 * t26;
        double t34 = 0.1e1 / t2 / B;
        double t73 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[1] * t10 + t4 * coef[2] * t15 + t18 * coef[6] * t21 + t18 * coef[7] * t27 + t4 * coef[5] * t31 + t34 * coef[29] * t7) + N[5] * (t4 * coef[8] * t10 + t4 * coef[9] * t31 + t34 * coef[3] * t7 + t4 * coef[10] * t15) + N[6] * t34 * t7 * coef[0] + N[2] * (t34 * coef[37] * t31 + t34 * coef[42] * t10 + t4 * coef[25] * t21 + t4 * coef[26] * t27 + t73 * coef[16] * t7 + t34 * coef[32] * t15) + N[3] * (t73 * coef[38] * t10 + t73 * coef[40] * t31 + t34 * coef[27] * t21 + t34 * coef[28] * t27 + t73 * coef[31] * t15) + N[11] * (t73 * coef[17] * t10 + t73 * coef[18] * t31 + t73 * coef[19] * t15) + N[4] * (t73 * coef[21] * t21 + t73 * coef[22] * t27) + N[1] * (t4 * coef[11] * t31 + t4 * coef[12] * t10 + t18 * coef[35] * t21 + t18 * coef[36] * t27 + t34 * coef[4] * t7 + t4 * coef[13] * t15) + N[8] * t34 * t7 * coef[14] + N[10] * t73 * t7 * coef[33] + N[9] * (t34 * coef[23] * t10 + t34 * coef[24] * t31 + t73 * coef[34] * t7 + t34 * coef[20] * t15) + N[7] * (t4 * coef[39] * t31 + t4 * coef[41] * t10 + t34 * coef[15] * t7 + t4 * coef[30] * t15);
        double t184 = Bm * Bm;
        double t185 = t184 * t184;
        double t187 = 0.1e1 / t185 / Bm;
        double t191 = 0.1e1 / t185;
        double t202 = 0.1e1 / t184 / Bm;
        double t247 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t187 * coefm[31] * t21 + t191 * coefm[28] * t31 + t191 * coefm[30] * t15 + t191 * coefm[42] * t10 + t202 * coefm[29] * t7 + t187 * coefm[41] * t27) + Nm[5] * (t191 * coefm[38] * t31 + t191 * coefm[39] * t10 + t202 * coefm[26] * t7 + t191 * coefm[40] * t15) + Nm[7] * (t191 * coefm[13] * t31 + t191 * coefm[24] * t10 + t202 * coefm[36] * t7 + t191 * coefm[17] * t15) + Nm[8] * t202 * t7 * coefm[37] + Nm[11] * (t247 * coefm[1] * t10 + t247 * coefm[2] * t31 + t247 * coefm[4] * t15) + Nm[6] * t202 * t7 * coefm[27] + Nm[2] * (t202 * coefm[20] * t10 + t202 * coefm[22] * t31 + t191 * coefm[9] * t27 + t191 * coefm[10] * t21 + t247 * coefm[35] * t7 + t202 * coefm[16] * t15) + Nm[3] * (t247 * coefm[14] * t10 + t247 * coefm[25] * t31 + t202 * coefm[11] * t27 + t202 * coefm[12] * t21 + t247 * coefm[23] * t15) + Nm[9] * (t202 * coefm[7] * t10 + t202 * coefm[8] * t31 + t247 * coefm[15] * t7 + t202 * coefm[3] * t15) + Nm[10] * t247 * t7 * coefm[18] + Nm[4] * (t247 * coefm[5] * t27 + t247 * coefm[6] * t21) + Nm[1] * (t191 * coefm[32] * t31 + t191 * coefm[33] * t10 + t187 * coefm[19] * t27 + t187 * coefm[21] * t21 + t202 * coefm[0] * t7 + t191 * coefm[34] * t15);
        //
        X = t353 + t182;
    } // if (flag == 13)
    
    if (flag == 21)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = t7 * t7;
        double t14 = cos(theta);
        double t15 = t13 * t14;
        double t19 = sin(theta);
        double t20 = t13 * t19;
        double t23 = 0.1e1 / t9 / B;
        double t26 = sin(Psi);
        double t28 = t13 * t26 * t19;
        double t32 = t7 * t26;
        double t36 = t15 * t26;
        double t61 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[1] * t15 + t10 * coef[2] * t20 + t23 * coef[7] * t28 + t10 * coef[5] * t32 + t23 * coef[6] * t36) + N[5] * (t10 * coef[8] * t15 + t10 * coef[9] * t20 + t4 * coef[3] * t7 + t10 * coef[10] * t32) + N[6] * t4 * t7 * coef[4] + N[10] * t61 * t7 * coef[32] + N[9] * (t4 * coef[23] * t15 + t4 * coef[24] * t20 + t61 * coef[33] * t7 + t4 * coef[19] * t32) + N[11] * (t61 * coef[17] * t15 + t61 * coef[18] * t20 + t61 * coef[20] * t32) + N[4] * (t61 * coef[21] * t36 + t61 * coef[22] * t28) + N[1] * (t10 * coef[11] * t20 + t10 * coef[12] * t15 + t23 * coef[37] * t36 + t23 * coef[38] * t28 + t4 * coef[0] * t7 + t10 * coef[13] * t32) + N[8] * t4 * t7 * coef[14] + N[7] * (t10 * coef[30] * t20 + t10 * coef[39] * t15 + t4 * coef[15] * t7 + t10 * coef[34] * t32) + N[3] * (t61 * coef[40] * t15 + t61 * coef[41] * t20 + t4 * coef[27] * t36 + t4 * coef[28] * t28 + t61 * coef[35] * t32) + N[2] * (t4 * coef[36] * t20 + t4 * coef[42] * t15 + t10 * coef[25] * t36 + t10 * coef[26] * t28 + t61 * coef[16] * t7 + t4 * coef[31] * t32);
        double t184 = Bm * Bm;
        double t185 = t184 * t184;
        double t186 = 0.1e1 / t185;
        double t194 = 0.1e1 / t185 / Bm;
        double t202 = 0.1e1 / t184 / Bm;
        double t252 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t186 * coefm[6] * t32 + t186 * coefm[9] * t20 + t194 * coefm[42] * t36 + t186 * coefm[4] * t15 + t202 * coefm[7] * t7 + t194 * coefm[41] * t28) + Nm[5] * (t186 * coefm[37] * t15 + t186 * coefm[38] * t20 + t202 * coefm[8] * t7 + t186 * coefm[0] * t32) + Nm[8] * t202 * t7 * coefm[1] + Nm[7] * (t186 * coefm[26] * t15 + t186 * coefm[35] * t20 + t202 * coefm[2] * t7 + t186 * coefm[32] * t32) + Nm[6] * t202 * t7 * coefm[5] + Nm[11] * (t252 * coefm[10] * t15 + t252 * coefm[11] * t20 + t252 * coefm[13] * t32) + Nm[9] * (t202 * coefm[16] * t15 + t202 * coefm[17] * t20 + t252 * coefm[24] * t7 + t202 * coefm[12] * t32) + Nm[10] * t252 * t7 * coefm[33] + Nm[4] * (t252 * coefm[14] * t28 + t252 * coefm[15] * t36) + Nm[1] * (t186 * coefm[39] * t20 + t186 * coefm[40] * t15 + t194 * coefm[28] * t28 + t194 * coefm[29] * t36 + t202 * coefm[22] * t7 + t186 * coefm[36] * t32) + Nm[3] * (t252 * coefm[27] * t15 + t252 * coefm[34] * t20 + t202 * coefm[20] * t28 + t202 * coefm[21] * t36 + t252 * coefm[31] * t32) + Nm[2] * (t202 * coefm[23] * t15 + t202 * coefm[30] * t20 + t186 * coefm[18] * t28 + t186 * coefm[19] * t36 + t252 * coefm[3] * t7 + t202 * coefm[25] * t32);
        //
        X = t353 + t182;
    } // if (flag == 21)
    
    if (flag == 22)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = t7 * t7;
        double t14 = sin(theta);
        double t15 = t13 * t14;
        double t19 = sin(Psi);
        double t20 = t7 * t19;
        double t23 = 0.1e1 / t9 / B;
        double t26 = cos(theta);
        double t27 = t13 * t26;
        double t28 = t27 * t19;
        double t36 = t13 * t19 * t14;
        double t58 = 0.1e1 / t2;
        complex<double> t173 = N[0] * (t4 * coef[27] * t7 + t10 * coef[0] * t15 + t10 * coef[1] * t20 + t23 * coef[18] * t28 + t10 * coef[4] * t27 + t23 * coef[19] * t36) + N[5] * (t10 * coef[20] * t15 + t10 * coef[21] * t27 + t10 * coef[17] * t20) + N[6] * t4 * t7 * coef[3] + N[10] * t58 * t7 * coef[35] + N[7] * (t10 * coef[28] * t27 + t10 * coef[34] * t15 + t10 * coef[37] * t20) + N[11] * (t58 * coef[5] * t15 + t58 * coef[6] * t27 + t58 * coef[8] * t20) + N[9] * (t4 * coef[11] * t15 + t4 * coef[12] * t27 + t4 * coef[7] * t20) + N[4] * (t58 * coef[9] * t28 + t58 * coef[10] * t36) + N[2] * (t4 * coef[30] * t27 + t4 * coef[31] * t15 + t10 * coef[13] * t28 + t10 * coef[14] * t36 + t58 * coef[24] * t7 + t4 * coef[32] * t20) + N[3] * (t58 * coef[36] * t27 + t58 * coef[38] * t15 + t4 * coef[15] * t28 + t4 * coef[16] * t36 + t58 * coef[33] * t20) + N[1] * (t10 * coef[25] * t15 + t10 * coef[26] * t27 + t23 * coef[29] * t28 + t23 * coef[39] * t36 + t4 * coef[2] * t7 + t10 * coef[22] * t20) + N[8] * t4 * t7 * coef[23];
        double t175 = Bm * Bm;
        double t176 = t175 * t175;
        double t177 = 0.1e1 / t176;
        double t182 = 0.1e1 / t176 / Bm;
        double t193 = 0.1e1 / t175 / Bm;
        double t220 = 0.1e1 / t175;
        complex<double> t335 = Nm[0] * (t177 * coefm[23] * t27 + t182 * coefm[12] * t36 + t182 * coefm[13] * t28 + t177 * coefm[37] * t20 + t193 * coefm[39] * t7 + t177 * coefm[22] * t15) + Nm[5] * (t177 * coefm[18] * t15 + t177 * coefm[19] * t27 + t177 * coefm[17] * t20) + Nm[8] * t193 * t7 * coefm[15] + Nm[11] * (t220 * coefm[0] * t15 + t220 * coefm[1] * t27 + t220 * coefm[3] * t20) + Nm[9] * (t193 * coefm[6] * t15 + t193 * coefm[7] * t27 + t193 * coefm[2] * t20) + Nm[4] * (t220 * coefm[4] * t36 + t220 * coefm[5] * t28) + Nm[2] * (t193 * coefm[28] * t15 + t193 * coefm[34] * t27 + t177 * coefm[8] * t36 + t177 * coefm[9] * t28 + t220 * coefm[16] * t7 + t193 * coefm[26] * t20) + Nm[3] * (t220 * coefm[25] * t27 + t220 * coefm[32] * t15 + t193 * coefm[10] * t36 + t193 * coefm[11] * t28 + t220 * coefm[30] * t20) + Nm[1] * (t177 * coefm[20] * t15 + t177 * coefm[21] * t27 + t182 * coefm[29] * t28 + t182 * coefm[31] * t36 + t193 * coefm[38] * t7 + t177 * coefm[14] * t20) + Nm[7] * (t177 * coefm[27] * t27 + t177 * coefm[35] * t15 + t177 * coefm[33] * t20) + Nm[10] * t220 * t7 * coefm[24] + Nm[6] * t193 * t7 * coefm[36];
        //
        X = t335 + t173;
        
    } // if (flag == 22)
    
    if (flag == 23)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = t7 * t7;
        double t14 = sin(theta);
        double t15 = t13 * t14;
        double t19 = cos(theta);
        double t20 = t13 * t19;
        double t23 = 0.1e1 / t9 / B;
        double t26 = sin(Psi);
        double t27 = t20 * t26;
        double t31 = t7 * t26;
        double t36 = t13 * t26 * t14;
        double t61 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[4] * t15 + t10 * coef[5] * t20 + t23 * coef[7] * t27 + t10 * coef[3] * t31 + t23 * coef[8] * t36) + N[5] * (t10 * coef[9] * t20 + t10 * coef[10] * t15 + t4 * coef[0] * t7 + t10 * coef[11] * t31) + N[6] * t4 * t7 * coef[1] + N[11] * (t61 * coef[17] * t15 + t61 * coef[20] * t20 + t61 * coef[21] * t31) + N[4] * (t61 * coef[18] * t27 + t61 * coef[22] * t36) + N[1] * (t10 * coef[12] * t15 + t10 * coef[13] * t20 + t23 * coef[32] * t36 + t23 * coef[34] * t27 + t4 * coef[2] * t7 + t10 * coef[14] * t31) + N[8] * t4 * t7 * coef[15] + N[10] * t61 * t7 * coef[31] + N[2] * (t4 * coef[37] * t15 + t4 * coef[42] * t20 + t10 * coef[25] * t27 + t10 * coef[26] * t36 + t61 * coef[6] * t7 + t4 * coef[33] * t31) + N[9] * (t4 * coef[23] * t20 + t4 * coef[24] * t15 + t61 * coef[35] * t7 + t4 * coef[19] * t31) + N[3] * (t61 * coef[40] * t15 + t61 * coef[41] * t20 + t4 * coef[27] * t27 + t4 * coef[28] * t36 + t61 * coef[36] * t31) + N[7] * (t10 * coef[38] * t15 + t10 * coef[39] * t20 + t4 * coef[16] * t7 + t10 * coef[30] * t31);
        double t184 = Bm * Bm;
        double t185 = t184 * t184;
        double t186 = 0.1e1 / t185;
        double t191 = 0.1e1 / t185 / Bm;
        double t205 = 0.1e1 / t184 / Bm;
        double t232 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t186 * coefm[13] * t20 + t191 * coefm[26] * t27 + t191 * coefm[27] * t36 + t186 * coefm[16] * t31 + t186 * coefm[12] * t15 + t205 * coefm[29] * t7) + Nm[5] * (t186 * coefm[17] * t15 + t186 * coefm[18] * t20 + t205 * coefm[14] * t7 + t186 * coefm[23] * t31) + Nm[8] * t205 * t7 * coefm[22] + Nm[11] * (t232 * coefm[0] * t20 + t232 * coefm[1] * t15 + t232 * coefm[2] * t31) + Nm[4] * (t232 * coefm[4] * t36 + t232 * coefm[5] * t27) + Nm[1] * (t186 * coefm[24] * t20 + t186 * coefm[25] * t15 + t191 * coefm[31] * t27 + t191 * coefm[33] * t36 + t205 * coefm[28] * t7 + t186 * coefm[19] * t31) + Nm[6] * t205 * t7 * coefm[15] + Nm[7] * (t186 * coefm[30] * t20 + t186 * coefm[34] * t15 + t205 * coefm[21] * t7 + t186 * coefm[42] * t31) + Nm[9] * (t205 * coefm[6] * t20 + t205 * coefm[7] * t15 + t232 * coefm[36] * t7 + t205 * coefm[3] * t31) + Nm[10] * t232 * t7 * coefm[41] + Nm[2] * (t205 * coefm[39] * t20 + t205 * coefm[40] * t15 + t186 * coefm[8] * t36 + t186 * coefm[9] * t27 + t232 * coefm[20] * t7 + t205 * coefm[38] * t31) + Nm[3] * (t232 * coefm[35] * t15 + t232 * coefm[37] * t20 + t205 * coefm[10] * t36 + t205 * coefm[11] * t27 + t232 * coefm[32] * t31);
        //
        X = t353 + t182;
        
    } // if (flag == 23)
    
    if (flag == 31)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = sin(Psi);
        double t14 = t7 * t13;
        double t18 = t7 * t7;
        double t19 = cos(theta);
        double t20 = t18 * t19;
        double t23 = 0.1e1 / t9 / B;
        double t27 = sin(theta);
        double t28 = t18 * t13 * t27;
        double t32 = t18 * t27;
        double t36 = t20 * t13;
        double t61 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[0] * t14 + t10 * coef[1] * t20 + t23 * coef[7] * t28 + t10 * coef[5] * t32 + t23 * coef[6] * t36) + N[5] * (t10 * coef[8] * t20 + t10 * coef[9] * t32 + t4 * coef[2] * t7 + t10 * coef[10] * t14) + N[6] * t4 * t7 * coef[3] + N[11] * (t61 * coef[17] * t20 + t61 * coef[18] * t32 + t61 * coef[19] * t14) + N[4] * (t61 * coef[21] * t36 + t61 * coef[22] * t28) + N[1] * (t10 * coef[11] * t32 + t10 * coef[12] * t20 + t23 * coef[30] * t28 + t23 * coef[34] * t36 + t4 * coef[4] * t7 + t10 * coef[13] * t14) + N[8] * t4 * t7 * coef[14] + N[2] * (t4 * coef[31] * t32 + t4 * coef[35] * t20 + t10 * coef[25] * t36 + t10 * coef[26] * t28 + t61 * coef[16] * t7 + t4 * coef[32] * t14) + N[10] * t61 * t7 * coef[33] + N[3] * (t61 * coef[37] * t20 + t61 * coef[42] * t32 + t4 * coef[27] * t36 + t4 * coef[28] * t28 + t61 * coef[36] * t14) + N[7] * (t10 * coef[38] * t32 + t10 * coef[40] * t20 + t4 * coef[15] * t7 + t10 * coef[41] * t14) + N[9] * (t4 * coef[23] * t20 + t4 * coef[24] * t32 + t61 * coef[39] * t7 + t4 * coef[20] * t14);
        double t184 = Bm * Bm;
        double t186 = 0.1e1 / t184 / Bm;
        double t190 = t184 * t184;
        double t191 = 0.1e1 / t190;
        double t199 = 0.1e1 / t190 / Bm;
        double t233 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t186 * coefm[28] * t7 + t191 * coefm[15] * t14 + t191 * coefm[29] * t20 + t199 * coefm[38] * t28 + t199 * coefm[39] * t36 + t191 * coefm[13] * t32) + Nm[5] * (t191 * coefm[36] * t20 + t191 * coefm[37] * t32 + t186 * coefm[14] * t7 + t191 * coefm[31] * t14) + Nm[9] * (t186 * coefm[22] * t20 + t186 * coefm[23] * t32 + t233 * coefm[4] * t7 + t186 * coefm[19] * t14) + Nm[10] * t233 * t7 * coefm[11] + Nm[11] * (t233 * coefm[16] * t20 + t233 * coefm[17] * t32 + t233 * coefm[18] * t14) + Nm[8] * t186 * t7 * coefm[33] + Nm[7] * (t191 * coefm[6] * t20 + t191 * coefm[12] * t32 + t186 * coefm[34] * t7 + t191 * coefm[3] * t14) + Nm[4] * (t233 * coefm[20] * t28 + t233 * coefm[21] * t36) + Nm[1] * (t191 * coefm[40] * t32 + t191 * coefm[41] * t20 + t199 * coefm[1] * t36 + t199 * coefm[10] * t28 + t186 * coefm[42] * t7 + t191 * coefm[35] * t14) + Nm[6] * t186 * t7 * coefm[30] + Nm[2] * (t186 * coefm[7] * t20 + t186 * coefm[8] * t32 + t191 * coefm[24] * t28 + t191 * coefm[25] * t36 + t233 * coefm[32] * t7 + t186 * coefm[0] * t14) + Nm[3] * (t233 * coefm[2] * t20 + t233 * coefm[5] * t32 + t186 * coefm[26] * t28 + t186 * coefm[27] * t36 + t233 * coefm[9] * t14);
        //
        X = t353 + t182;
        
    } // if (flag == 31)
    
    if (flag == 32)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = t7 * t7;
        double t14 = sin(theta);
        double t15 = t13 * t14;
        double t18 = 0.1e1 / t9 / B;
        double t21 = cos(theta);
        double t22 = t13 * t21;
        double t23 = sin(Psi);
        double t24 = t22 * t23;
        double t28 = t7 * t23;
        double t36 = t13 * t23 * t14;
        double t82 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[29] * t7 + t10 * coef[1] * t15 + t18 * coef[7] * t24 + t10 * coef[3] * t28 + t10 * coef[5] * t22 + t18 * coef[6] * t36) + N[5] * (t10 * coef[8] * t22 + t10 * coef[9] * t15 + t4 * coef[0] * t7 + t10 * coef[10] * t28) + N[6] * t4 * t7 * coef[2] + N[1] * (t10 * coef[11] * t15 + t10 * coef[12] * t22 + t18 * coef[30] * t24 + t18 * coef[33] * t36 + t4 * coef[4] * t7 + t10 * coef[13] * t28) + N[11] * (t82 * coef[17] * t22 + t82 * coef[18] * t15 + t82 * coef[20] * t28) + N[4] * (t82 * coef[21] * t36 + t82 * coef[22] * t24) + N[8] * t4 * t7 * coef[14] + N[10] * t82 * t7 * coef[31] + N[7] * (t10 * coef[41] * t15 + t10 * coef[42] * t22 + t4 * coef[15] * t7 + t10 * coef[32] * t28) + N[3] * (t82 * coef[39] * t22 + t82 * coef[40] * t15 + t4 * coef[27] * t36 + t4 * coef[28] * t24 + t82 * coef[34] * t28) + N[2] * (t4 * coef[37] * t15 + t4 * coef[38] * t22 + t10 * coef[25] * t36 + t10 * coef[26] * t24 + t82 * coef[16] * t7 + t4 * coef[35] * t28) + N[9] * (t4 * coef[23] * t22 + t4 * coef[24] * t15 + t82 * coef[36] * t7 + t4 * coef[19] * t28);
        double t184 = Bm * Bm;
        double t185 = t184 * t184;
        double t186 = 0.1e1 / t185;
        double t194 = 0.1e1 / t184 / Bm;
        double t202 = 0.1e1 / t185 / Bm;
        double t258 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t186 * coefm[16] * t28 + t186 * coefm[17] * t15 + t194 * coefm[1] * t7 + t186 * coefm[0] * t22 + t202 * coefm[37] * t36 + t202 * coefm[38] * t24) + Nm[5] * (t186 * coefm[34] * t22 + t186 * coefm[35] * t15 + t194 * coefm[15] * t7 + t186 * coefm[36] * t28) + Nm[8] * t194 * t7 * coefm[39] + Nm[6] * t194 * t7 * coefm[18] + Nm[1] * (t186 * coefm[32] * t15 + t186 * coefm[33] * t22 + t202 * coefm[4] * t24 + t202 * coefm[6] * t36 + t194 * coefm[19] * t7 + t186 * coefm[42] * t28) + Nm[4] * (t258 * coefm[24] * t36 + t258 * coefm[25] * t24) + Nm[11] * (t258 * coefm[20] * t22 + t258 * coefm[21] * t15 + t258 * coefm[23] * t28) + Nm[3] * (t258 * coefm[7] * t22 + t258 * coefm[13] * t15 + t194 * coefm[30] * t36 + t194 * coefm[31] * t24 + t258 * coefm[3] * t28) + Nm[2] * (t194 * coefm[8] * t15 + t194 * coefm[9] * t22 + t186 * coefm[28] * t36 + t186 * coefm[29] * t24 + t258 * coefm[41] * t7 + t194 * coefm[10] * t28) + Nm[9] * (t194 * coefm[26] * t22 + t194 * coefm[27] * t15 + t258 * coefm[5] * t7 + t194 * coefm[22] * t28) + Nm[7] * (t186 * coefm[11] * t15 + t186 * coefm[14] * t22 + t194 * coefm[40] * t7 + t186 * coefm[2] * t28) + Nm[10] * t258 * t7 * coefm[12];
        //
        X = t353 + t182;
        
    } // if (flag == 32)
    
    if (flag == 33)
    {
        double t2 = B * B;
        double t4 = 0.1e1 / t2 / B;
        double t7 = cos(Psi);
        double t9 = t2 * t2;
        double t10 = 0.1e1 / t9;
        double t13 = sin(Psi);
        double t14 = t7 * t13;
        double t17 = 0.1e1 / t9 / B;
        double t20 = t7 * t7;
        double t21 = cos(theta);
        double t22 = t20 * t21;
        double t23 = t22 * t13;
        double t30 = sin(theta);
        double t31 = t20 * t30;
        double t36 = t20 * t13 * t30;
        double t61 = 0.1e1 / t2;
        complex<double> t182 = N[0] * (t4 * coef[21] * t7 + t10 * coef[1] * t14 + t17 * coef[6] * t23 + t10 * coef[3] * t22 + t10 * coef[4] * t31 + t17 * coef[7] * t36) + N[5] * (t10 * coef[8] * t22 + t10 * coef[35] * t31 + t4 * coef[2] * t7 + t10 * coef[36] * t14) + N[6] * t4 * t7 * coef[5] + N[11] * (t61 * coef[9] * t22 + t61 * coef[10] * t31 + t61 * coef[12] * t14) + N[4] * (t61 * coef[13] * t23 + t61 * coef[14] * t36) + N[1] * (t10 * coef[37] * t31 + t10 * coef[38] * t22 + t17 * coef[27] * t23 + t17 * coef[28] * t36 + t4 * coef[0] * t7 + t10 * coef[39] * t14) + N[8] * t4 * t7 * coef[40] + N[9] * (t4 * coef[15] * t22 + t4 * coef[16] * t31 + t61 * coef[24] * t7 + t4 * coef[11] * t14) + N[3] * (t61 * coef[29] * t22 + t61 * coef[31] * t31 + t4 * coef[19] * t23 + t4 * coef[20] * t36 + t61 * coef[25] * t14) + N[7] * (t10 * coef[32] * t22 + t10 * coef[33] * t31 + t4 * coef[41] * t7 + t10 * coef[26] * t14) + N[2] * (t4 * coef[30] * t31 + t4 * coef[34] * t22 + t10 * coef[17] * t23 + t10 * coef[18] * t36 + t61 * coef[42] * t7 + t4 * coef[22] * t14) + N[10] * t61 * t7 * coef[23];
        double t184 = Bm * Bm;
        double t185 = t184 * t184;
        double t187 = 0.1e1 / t185 / Bm;
        double t194 = 0.1e1 / t185;
        double t199 = 0.1e1 / t184 / Bm;
        double t227 = 0.1e1 / t184;
        complex<double> t353 = Nm[0] * (t187 * coefm[32] * t36 + t187 * coefm[33] * t23 + t194 * coefm[0] * t22 + t199 * coefm[31] * t7 + t194 * coefm[15] * t14 + t194 * coefm[14] * t31) + Nm[5] * (t194 * coefm[36] * t31 + t194 * coefm[37] * t22 + t199 * coefm[28] * t7 + t194 * coefm[34] * t14) + Nm[10] * t227 * t7 * coefm[5] + Nm[9] * (t199 * coefm[22] * t22 + t199 * coefm[23] * t31 + t227 * coefm[9] * t7 + t199 * coefm[18] * t14) + Nm[7] * (t194 * coefm[6] * t22 + t194 * coefm[8] * t31 + t199 * coefm[40] * t7 + t194 * coefm[3] * t14) + Nm[8] * t199 * t7 * coefm[41] + Nm[6] * t199 * t7 * coefm[30] + Nm[4] * (t227 * coefm[20] * t36 + t227 * coefm[21] * t23) + Nm[1] * (t194 * coefm[35] * t31 + t194 * coefm[42] * t22 + t187 * coefm[2] * t23 + t187 * coefm[4] * t36 + t199 * coefm[29] * t7 + t194 * coefm[38] * t14) + Nm[11] * (t227 * coefm[16] * t22 + t227 * coefm[17] * t31 + t227 * coefm[19] * t14) + Nm[3] * (t227 * coefm[7] * t22 + t227 * coefm[12] * t31 + t199 * coefm[26] * t36 + t199 * coefm[27] * t23 + t227 * coefm[1] * t14) + Nm[2] * (t199 * coefm[11] * t22 + t199 * coefm[13] * t31 + t194 * coefm[24] * t36 + t194 * coefm[25] * t23 + t227 * coefm[39] * t7 + t199 * coefm[10] * t14);
        //
        X = t353 + t182;
        
    } // if (flag == 33)
    
    // Final Output
    return X;
	
}



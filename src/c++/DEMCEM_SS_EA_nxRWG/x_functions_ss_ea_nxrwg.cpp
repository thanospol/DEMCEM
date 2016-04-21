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
#include "demcem_ss_ea_nxrwg.h"


// ***********************************************************************
//			IMPLEMENTATION OF complex<double> x_functions
// ***********************************************************************

complex<double> x_functions_ss_ea_nxrwg (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[], int flag)

{
    complex<double> X;
    //
    
    if (flag == 11)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = sin(theta);
        double t10 = t8 * t9;
        double t13 = 0.1e1 / t3 / t2;
        double t19 = t8 * t8;
        double t20 = t19 * t9;
        double t25 = sin(Psi);
        double t28 = 0.1e1 / t3 / B;
        double t40 = 0.1e1 / t2 / B;
        complex<double> t110 = N[0] * (t4 * coef[0] * t10 + t13 * coef[1] * t10 + t13 * coef[2] * t20) + N[5] * t25 * t28 * t10 * coef[3] + N[8] * t4 * t10 * coef[4] + N[2] * (t40 * coef[9] * t10 + t28 * coef[6] * t20 + t28 * coef[5] * t10) + N[3] * (t4 * coef[8] * t20 + t4 * coef[7] * t10) + N[4] * (t40 * coef[10] * t20 + t40 * coef[11] * t10) + N[7] * t25 * t28 * t10 * coef[14] + N[10] * t40 * t10 * coef[15] + N[9] * t25 * t4 * t10 * coef[12] + N[11] * t25 * t40 * t10 * coef[13] + N[1] * (t4 * coef[17] * t10 + t13 * coef[19] * t20 + t13 * coef[18] * t10) + N[6] * t4 * t10 * coef[16];
        double t112 = Bm * Bm;
        double t113 = t112 * t112;
        double t114 = 0.1e1 / t113;
        double t119 = 0.1e1 / t113 / t112;
        double t135 = 0.1e1 / t112 / Bm;
        double t143 = 0.1e1 / t113 / Bm;
        complex<double> t213 = Nm[0] * (t114 * coefm[0] * t10 + t119 * coefm[1] * t10 + t119 * coefm[2] * t20) + Nm[8] * t114 * t10 * coefm[9] + Nm[10] * t135 * t10 * coefm[14] + Nm[7] * t25 * t143 * t10 * coefm[15] + Nm[2] * (t135 * coefm[7] * t10 + t143 * coefm[11] * t20 + t143 * coefm[10] * t10) + Nm[3] * (t114 * coefm[13] * t20 + t114 * coefm[12] * t10) + Nm[4] * (t135 * coefm[5] * t20 + t135 * coefm[6] * t10) + Nm[1] * (t114 * coefm[17] * t10 + t119 * coefm[19] * t20 + t119 * coefm[18] * t10) + Nm[5] * t25 * t143 * t10 * coefm[8] + Nm[6] * t114 * t10 * coefm[16] + Nm[9] * t25 * t114 * t10 * coefm[3] + Nm[11] * t25 * t135 * t10 * coefm[4];
        //
        X = t213 + t110;
    } // if (flag == 11)
    
    if (flag == 12)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / t2;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = sin(theta);
        double t11 = t9 * t10;
        double t15 = t9 * t9;
        double t16 = t15 * t10;
        double t18 = 0.1e1 / t3;
        double t23 = 0.1e1 / t3 / B;
        double t26 = sin(Psi);
        double t28 = t9 * t26 * t10;
        double t63 = 0.1e1 / t2 / B;
        complex<double> t133 = N[0] * (t5 * coef[0] * t11 + t5 * coef[7] * t16 + t18 * coef[8] * t11 + t23 * coef[3] * t28) + N[5] * (t18 * coef[2] * t11 + t23 * coef[17] * t28) + N[6] * t18 * t11 * coef[1] + N[1] * (t18 * coef[6] * t11 + t23 * coef[23] * t28 + t5 * coef[5] * t16 + t5 * coef[4] * t11) + N[4] * (t63 * coef[14] * t16 + t63 * coef[16] * t11) + N[3] * (t63 * coef[11] * t28 + t18 * coef[21] * t16 + t18 * coef[18] * t11) + N[8] * t18 * t11 * coef[19] + N[2] * (t63 * coef[15] * t11 + t18 * coef[13] * t28 + t23 * coef[20] * t16 + t23 * coef[24] * t11) + N[7] * (t18 * coef[22] * t11 + t23 * coef[12] * t28) + N[9] * (t63 * coef[10] * t11 + t18 * coef[26] * t28) + N[11] * t26 * t63 * t11 * coef[25] + N[10] * t63 * t11 * coef[9];
        double t135 = Bm * Bm;
        double t136 = t135 * t135;
        double t137 = 0.1e1 / t136;
        double t142 = 0.1e1 / t136 / t135;
        double t147 = 0.1e1 / t136 / Bm;
        double t167 = 0.1e1 / t135 / Bm;
        complex<double> t257 = Nm[0] * (t137 * coefm[22] * t11 + t142 * coefm[23] * t11 + t147 * coefm[18] * t28 + t142 * coefm[24] * t16) + Nm[5] * (t137 * coefm[19] * t11 + t147 * coefm[11] * t28) + Nm[9] * (t167 * coefm[0] * t11 + t137 * coefm[26] * t28) + Nm[10] * t167 * t11 * coefm[3] + Nm[6] * t137 * t11 * coefm[20] + Nm[7] * (t137 * coefm[9] * t11 + t147 * coefm[1] * t28) + Nm[8] * t137 * t11 * coefm[6] + Nm[1] * (t137 * coefm[16] * t11 + t147 * coefm[8] * t28 + t142 * coefm[21] * t16 + t142 * coefm[17] * t11) + Nm[11] * t26 * t167 * t11 * coefm[25] + Nm[4] * (t167 * coefm[15] * t16 + t167 * coefm[14] * t11) + Nm[3] * (t167 * coefm[4] * t28 + t137 * coefm[12] * t16 + t137 * coefm[7] * t11) + Nm[2] * (t167 * coefm[13] * t11 + t137 * coefm[2] * t28 + t147 * coefm[10] * t16 + t147 * coefm[5] * t11);
        //
        X = t257 + t133;
        
    } // if (flag == 12)
    
    
    if (flag == 13)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / t2;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = t9 * t9;
        double t11 = sin(theta);
        double t12 = t10 * t11;
        double t14 = 0.1e1 / t3;
        double t17 = cos(theta);
        double t18 = t9 * t17;
        double t21 = 0.1e1 / t3 / B;
        double t24 = sin(Psi);
        double t25 = t18 * t24;
        double t30 = t9 * t24 * t11;
        double t34 = t9 * t8;
        double t41 = t8 * t24;
        double t45 = t9 * t11;
        double t78 = t14 * t9;
        double t84 = 0.1e1 / t2 / B;
        double t136 = t84 * t9;
        complex<double> t203 = N[0] * (t5 * coef[47] * t12 + t14 * coef[17] * t18 + t21 * coef[32] * t25 + t21 * coef[19] * t30 + t21 * coef[15] * t34 + t21 * coef[16] * t8 + t14 * coef[45] * t41 + t5 * coef[33] * t45) + N[1] * (t14 * coef[27] * t18 + t21 * coef[2] * t25 + t21 * coef[8] * t30 + t5 * coef[28] * t12 + t14 * coef[29] * t41 + t21 * coef[46] * t34 + t5 * coef[25] * t45 + t21 * coef[22] * t8) + N[8] * t11 * t78 * coef[12] + N[4] * (t84 * coef[43] * t12 + t84 * coef[44] * t45) + N[5] * (t14 * coef[18] * t18 + t14 * coef[24] * t45 + t21 * coef[9] * t30 + t14 * coef[31] * t41) + N[6] * t11 * t78 * coef[23] + N[3] * (t84 * coef[39] * t25 + t84 * coef[40] * t30 + t14 * coef[6] * t12 + t84 * coef[20] * t34 + t14 * coef[3] * t45 + t84 * coef[21] * t8) + N[10] * t11 * t136 * coef[42] + N[9] * (t84 * coef[35] * t45 + t84 * coef[41] * t18 + t14 * coef[0] * t30 + t84 * coef[34] * t41) + N[7] * (t14 * coef[5] * t18 + t14 * coef[13] * t45 + t21 * coef[38] * t30 + t14 * coef[4] * t41) + N[2] * (t84 * coef[11] * t18 + t14 * coef[36] * t25 + t14 * coef[37] * t30 + t21 * coef[14] * t12 + t84 * coef[10] * t41 + t14 * coef[26] * t34 + t21 * coef[7] * t45 + t14 * coef[30] * t8) + N[11] * t11 * t24 * t136 * coef[1];
        double t205 = Bm * Bm;
        double t206 = t205 * t205;
        double t207 = 0.1e1 / t206;
        double t212 = 0.1e1 / t206 / Bm;
        double t223 = 0.1e1 / t206 / t205;
        double t255 = 0.1e1 / t205 / Bm;
        double t293 = t207 * t9;
        double t382 = t255 * t9;
        complex<double> t392 = Nm[0] * (t207 * coefm[38] * t41 + t212 * coefm[23] * t30 + t212 * coefm[35] * t34 + t212 * coefm[34] * t8 + t223 * coefm[36] * t12 + t223 * coefm[37] * t45 + t207 * coefm[33] * t18 + t212 * coefm[19] * t25) + Nm[5] * (t207 * coefm[31] * t18 + t207 * coefm[32] * t45 + t212 * coefm[13] * t30 + t207 * coefm[26] * t41) + Nm[4] * (t255 * coefm[8] * t12 + t255 * coefm[9] * t45) + Nm[1] * (t207 * coefm[22] * t18 + t212 * coefm[11] * t25 + t212 * coefm[15] * t30 + t223 * coefm[25] * t12 + t207 * coefm[30] * t41 + t212 * coefm[18] * t34 + t223 * coefm[21] * t45 + t212 * coefm[2] * t8) + Nm[8] * t11 * t293 * coefm[10] + Nm[6] * t11 * t293 * coefm[20] + Nm[2] * (t255 * coefm[7] * t18 + t207 * coefm[39] * t25 + t207 * coefm[40] * t30 + t212 * coefm[14] * t12 + t255 * coefm[6] * t41 + t207 * coefm[24] * t34 + t212 * coefm[17] * t45 + t207 * coefm[29] * t8) + Nm[3] * (t255 * coefm[44] * t30 + t255 * coefm[47] * t25 + t207 * coefm[5] * t12 + t255 * coefm[27] * t34 + t207 * coefm[12] * t45 + t255 * coefm[28] * t8) + Nm[9] * (t255 * coefm[42] * t18 + t255 * coefm[45] * t45 + t207 * coefm[1] * t30 + t255 * coefm[43] * t41) + Nm[7] * (t207 * coefm[4] * t45 + t207 * coefm[16] * t18 + t212 * coefm[41] * t30 + t207 * coefm[3] * t41) + Nm[10] * t11 * t382 * coefm[46] + Nm[11] * t11 * t24 * t382 * coefm[0];
        //
        X = t392 + t203;
    } // if (flag == 13)
    
    if (flag == 21)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / t2;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = sin(theta);
        double t11 = t9 * t10;
        double t14 = 0.1e1 / t3 / B;
        double t17 = sin(Psi);
        double t19 = t9 * t17 * t10;
        double t21 = 0.1e1 / t3;
        double t27 = t9 * t9;
        double t28 = t27 * t10;
        double t48 = 0.1e1 / t2 / B;
        complex<double> t133 = N[0] * (t5 * coef[0] * t11 + t14 * coef[10] * t19 + t21 * coef[14] * t11 + t5 * coef[15] * t28) + N[5] * (t21 * coef[9] * t11 + t14 * coef[26] * t19) + N[6] * t21 * t11 * coef[8] + N[10] * t48 * t11 * coef[5] + N[7] * (t21 * coef[25] * t11 + t14 * coef[6] * t19) + N[1] * (t21 * coef[13] * t11 + t14 * coef[20] * t19 + t5 * coef[12] * t28 + t5 * coef[11] * t11) + N[3] * (t48 * coef[4] * t19 + t21 * coef[21] * t28 + t21 * coef[19] * t11) + N[2] * (t48 * coef[17] * t11 + t21 * coef[7] * t19 + t14 * coef[24] * t28 + t14 * coef[22] * t11) + N[8] * t21 * t11 * coef[23] + N[4] * (t48 * coef[18] * t28 + t48 * coef[16] * t11) + N[9] * (t48 * coef[3] * t11 + t21 * coef[2] * t19) + N[11] * t17 * t48 * t11 * coef[1];
        double t135 = Bm * Bm;
        double t136 = t135 * t135;
        double t138 = 0.1e1 / t136 / t135;
        double t146 = 0.1e1 / t136 / Bm;
        double t150 = 0.1e1 / t136;
        double t167 = 0.1e1 / t135 / Bm;
        complex<double> t257 = Nm[0] * (t138 * coefm[20] * t11 + t138 * coefm[17] * t28 + t146 * coefm[24] * t19 + t150 * coefm[16] * t11) + Nm[5] * (t150 * coefm[25] * t11 + t146 * coefm[3] * t19) + Nm[10] * t167 * t11 * coefm[14] + Nm[11] * t17 * t167 * t11 * coefm[19] + Nm[9] * (t167 * coefm[13] * t11 + t150 * coefm[18] * t19) + Nm[7] * (t150 * coefm[0] * t11 + t146 * coefm[11] * t19) + Nm[8] * t150 * t11 * coefm[1] + Nm[1] * (t150 * coefm[21] * t11 + t146 * coefm[4] * t19 + t138 * coefm[23] * t28 + t138 * coefm[22] * t11) + Nm[3] * (t167 * coefm[15] * t19 + t150 * coefm[5] * t28 + t150 * coefm[6] * t11) + Nm[2] * (t167 * coefm[9] * t11 + t150 * coefm[12] * t19 + t146 * coefm[7] * t28 + t146 * coefm[2] * t11) + Nm[4] * (t167 * coefm[10] * t28 + t167 * coefm[8] * t11) + Nm[6] * t150 * t11 * coefm[26];
        //
        X = t257 + t133;
        
    } // if (flag == 21)
    
    if (flag == 22)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / t2;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = t9 * t9;
        double t11 = sin(theta);
        double t12 = t10 * t11;
        double t16 = t9 * t11;
        double t18 = 0.1e1 / t3;
        double t25 = sin(Psi);
        double t28 = 0.1e1 / t3 / B;
        double t40 = 0.1e1 / t2 / B;
        complex<double> t110 = N[0] * (t5 * coef[17] * t12 + t5 * coef[2] * t16 + t18 * coef[16] * t16) + N[5] * t25 * t28 * t16 * coef[10] + N[6] * t18 * t16 * coef[3] + N[2] * (t40 * coef[7] * t16 + t28 * coef[13] * t12 + t28 * coef[12] * t16) + N[4] * (t40 * coef[9] * t12 + t40 * coef[8] * t16) + N[8] * t18 * t16 * coef[11] + N[3] * (t18 * coef[15] * t12 + t18 * coef[14] * t16) + N[10] * t40 * t16 * coef[0] + N[7] * t25 * t28 * t16 * coef[1] + N[1] * (t18 * coef[4] * t16 + t5 * coef[6] * t12 + t5 * coef[5] * t16) + N[9] * t25 * t18 * t16 * coef[18] + N[11] * t25 * t40 * t16 * coef[19];
        double t112 = Bm * Bm;
        double t113 = t112 * t112;
        double t114 = 0.1e1 / t113;
        double t119 = 0.1e1 / t113 / t112;
        double t131 = 0.1e1 / t113 / Bm;
        double t138 = 0.1e1 / t112 / Bm;
        complex<double> t213 = Nm[0] * (t114 * coefm[0] * t16 + t119 * coefm[1] * t12 + t119 * coefm[11] * t16) + Nm[5] * t25 * t131 * t16 * coefm[2] + Nm[2] * (t138 * coefm[8] * t16 + t131 * coefm[5] * t12 + t131 * coefm[4] * t16) + Nm[4] * (t138 * coefm[10] * t12 + t138 * coefm[9] * t16) + Nm[3] * (t114 * coefm[7] * t12 + t114 * coefm[6] * t16) + Nm[10] * t138 * t16 * coefm[18] + Nm[7] * t25 * t131 * t16 * coefm[19] + Nm[6] * t114 * t16 * coefm[14] + Nm[9] * t25 * t114 * t16 * coefm[12] + Nm[11] * t25 * t138 * t16 * coefm[13] + Nm[1] * (t114 * coefm[15] * t16 + t119 * coefm[17] * t12 + t119 * coefm[16] * t16) + Nm[8] * t114 * t16 * coefm[3];
        //
        X = t213 + t110;
        
    } // if (flag == 22)
    
    if (flag == 23)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = cos(theta);
        double t10 = t8 * t9;
        double t13 = 0.1e1 / t3 / B;
        double t16 = t8 * t7;
        double t22 = 0.1e1 / t3 / t2;
        double t25 = sin(theta);
        double t26 = t8 * t25;
        double t30 = sin(Psi);
        double t31 = t7 * t30;
        double t35 = t8 * t8;
        double t36 = t35 * t25;
        double t41 = t8 * t30 * t25;
        double t45 = t10 * t30;
        double t51 = 0.1e1 / t2 / B;
        double t62 = t4 * t8;
        double t69 = t51 * t8;
        complex<double> t203 = N[0] * (t4 * coef[14] * t10 + t13 * coef[46] * t16 + t13 * coef[47] * t7 + t22 * coef[15] * t26 + t4 * coef[11] * t31 + t22 * coef[10] * t36 + t13 * coef[7] * t41 + t13 * coef[3] * t45) + N[4] * (t51 * coef[31] * t36 + t51 * coef[32] * t26) + N[8] * t25 * t62 * coef[18] + N[11] * t25 * t30 * t69 * coef[12] + N[5] * (t4 * coef[5] * t10 + t4 * coef[6] * t26 + t13 * coef[24] * t41 + t4 * coef[22] * t31) + N[6] * t25 * t62 * coef[4] + N[7] * (t4 * coef[25] * t26 + t4 * coef[26] * t10 + t13 * coef[37] * t41 + t4 * coef[23] * t31) + N[9] * (t51 * coef[41] * t26 + t51 * coef[44] * t10 + t4 * coef[13] * t41 + t51 * coef[38] * t31) + N[3] * (t51 * coef[40] * t45 + t51 * coef[42] * t41 + t4 * coef[29] * t36 + t51 * coef[1] * t16 + t4 * coef[34] * t26 + t51 * coef[0] * t7) + N[2] * (t51 * coef[19] * t10 + t4 * coef[39] * t45 + t4 * coef[45] * t41 + t13 * coef[33] * t36 + t51 * coef[20] * t31 + t4 * coef[21] * t16 + t13 * coef[30] * t26 + t4 * coef[36] * t7) + N[10] * t25 * t69 * coef[43] + N[1] * (t4 * coef[2] * t10 + t13 * coef[27] * t41 + t13 * coef[28] * t45 + t22 * coef[8] * t36 + t4 * coef[35] * t31 + t13 * coef[16] * t16 + t22 * coef[9] * t26 + t13 * coef[17] * t7);
        double t205 = Bm * Bm;
        double t206 = t205 * t205;
        double t207 = 0.1e1 / t206;
        double t212 = 0.1e1 / t206 / Bm;
        double t217 = 0.1e1 / t206 / t205;
        double t255 = 0.1e1 / t205 / Bm;
        double t266 = t207 * t8;
        double t273 = t255 * t8;
        complex<double> t392 = Nm[0] * (t207 * coefm[2] * t10 + t212 * coefm[28] * t45 + t217 * coefm[34] * t36 + t212 * coefm[32] * t41 + t212 * coefm[46] * t7 + t207 * coefm[45] * t31 + t217 * coefm[1] * t26 + t212 * coefm[47] * t16) + Nm[5] * (t207 * coefm[20] * t26 + t207 * coefm[33] * t10 + t212 * coefm[13] * t41 + t207 * coefm[24] * t31) + Nm[4] * (t255 * coefm[9] * t36 + t255 * coefm[12] * t26) + Nm[8] * t25 * t266 * coefm[14] + Nm[11] * t25 * t30 * t273 * coefm[6] + Nm[6] * t25 * t266 * coefm[19] + Nm[2] * (t255 * coefm[10] * t10 + t207 * coefm[35] * t41 + t207 * coefm[36] * t45 + t212 * coefm[18] * t36 + t255 * coefm[11] * t31 + t207 * coefm[25] * t16 + t212 * coefm[7] * t26 + t207 * coefm[22] * t7) + Nm[3] * (t255 * coefm[40] * t45 + t255 * coefm[43] * t41 + t207 * coefm[5] * t36 + t255 * coefm[31] * t16 + t207 * coefm[8] * t26 + t255 * coefm[30] * t7) + Nm[10] * t25 * t273 * coefm[42] + Nm[9] * (t255 * coefm[39] * t26 + t255 * coefm[41] * t10 + t207 * coefm[0] * t41 + t255 * coefm[38] * t31) + Nm[7] * (t207 * coefm[4] * t10 + t207 * coefm[17] * t26 + t212 * coefm[37] * t41 + t207 * coefm[15] * t31) + Nm[1] * (t207 * coefm[29] * t10 + t212 * coefm[3] * t45 + t212 * coefm[16] * t41 + t217 * coefm[21] * t36 + t207 * coefm[26] * t31 + t212 * coefm[27] * t16 + t217 * coefm[23] * t26 + t212 * coefm[44] * t7);
        //
        X = t392 + t203;
        
    } // if (flag == 23)
    
    if (flag == 31)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t4 = 0.1e1 / t3;
        double t7 = cos(Psi);
        double t8 = t7 * t7;
        double t9 = sin(theta);
        double t10 = t8 * t9;
        double t13 = 0.1e1 / t3 / B;
        double t16 = sin(Psi);
        double t18 = t8 * t16 * t9;
        double t21 = 0.1e1 / t3 / t2;
        double t24 = t8 * t8;
        double t25 = t24 * t9;
        double t39 = 0.1e1 / t2 / B;
        complex<double> t133 = N[0] * (t4 * coef[25] * t10 + t13 * coef[7] * t18 + t21 * coef[26] * t25 + t21 * coef[11] * t10) + N[8] * t4 * t10 * coef[19] + N[4] * (t39 * coef[14] * t25 + t39 * coef[12] * t10) + N[5] * (t4 * coef[6] * t10 + t13 * coef[17] * t18) + N[6] * t4 * t10 * coef[5] + N[11] * t16 * t39 * t10 * coef[23] + N[10] * t39 * t10 * coef[2] + N[3] * (t39 * coef[1] * t18 + t4 * coef[18] * t25 + t4 * coef[15] * t10) + N[2] * (t39 * coef[13] * t10 + t4 * coef[4] * t18 + t13 * coef[16] * t25 + t13 * coef[21] * t10) + N[7] * (t4 * coef[22] * t10 + t13 * coef[3] * t18) + N[9] * (t39 * coef[0] * t10 + t4 * coef[24] * t18) + N[1] * (t4 * coef[9] * t10 + t13 * coef[20] * t18 + t21 * coef[10] * t25 + t21 * coef[8] * t10);
        double t135 = Bm * Bm;
        double t136 = t135 * t135;
        double t138 = 0.1e1 / t136 / t135;
        double t143 = 0.1e1 / t136 / Bm;
        double t147 = 0.1e1 / t136;
        double t158 = 0.1e1 / t135 / Bm;
        complex<double> t257 = Nm[0] * (t138 * coefm[13] * t10 + t143 * coefm[7] * t18 + t147 * coefm[14] * t10 + t138 * coefm[15] * t25) + Nm[4] * (t158 * coefm[26] * t25 + t158 * coefm[24] * t10) + Nm[5] * (t147 * coefm[9] * t10 + t143 * coefm[23] * t18) + Nm[6] * t147 * t10 * coefm[8] + Nm[8] * t147 * t10 * coefm[22] + Nm[11] * t16 * t158 * t10 * coefm[0] + Nm[3] * (t158 * coefm[6] * t18 + t147 * coefm[18] * t25 + t147 * coefm[20] * t10) + Nm[2] * (t158 * coefm[25] * t10 + t147 * coefm[3] * t18 + t143 * coefm[16] * t25 + t143 * coefm[21] * t10) + Nm[1] * (t147 * coefm[12] * t10 + t143 * coefm[17] * t18 + t138 * coefm[11] * t25 + t138 * coefm[10] * t10) + Nm[10] * t158 * t10 * coefm[4] + Nm[9] * (t158 * coefm[5] * t10 + t147 * coefm[1] * t18) + Nm[7] * (t147 * coefm[19] * t10 + t143 * coefm[2] * t18);
        //
        X = t257 + t133;
        
    } // if (flag == 31)
    
    if (flag == 32)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / t2;
        double t8 = cos(Psi);
        double t9 = t8 * t8;
        double t10 = t9 * t9;
        double t11 = sin(theta);
        double t12 = t10 * t11;
        double t16 = t9 * t11;
        double t19 = 0.1e1 / t3 / B;
        double t22 = sin(Psi);
        double t24 = t9 * t22 * t11;
        double t26 = 0.1e1 / t3;
        double t34 = 0.1e1 / t2 / B;
        complex<double> t133 = N[0] * (t5 * coef[15] * t12 + t5 * coef[0] * t16 + t19 * coef[10] * t24 + t26 * coef[14] * t16) + N[2] * (t34 * coef[17] * t16 + t26 * coef[7] * t24 + t19 * coef[26] * t12 + t19 * coef[24] * t16) + N[4] * (t34 * coef[18] * t12 + t34 * coef[16] * t16) + N[8] * t26 * t16 * coef[22] + N[3] * (t34 * coef[4] * t24 + t26 * coef[23] * t12 + t26 * coef[20] * t16) + N[5] * (t26 * coef[9] * t16 + t19 * coef[25] * t24) + N[6] * t26 * t16 * coef[8] + N[11] * t22 * t34 * t16 * coef[1] + N[9] * (t34 * coef[5] * t16 + t26 * coef[2] * t24) + N[7] * (t26 * coef[21] * t16 + t19 * coef[6] * t24) + N[10] * t34 * t16 * coef[3] + N[1] * (t26 * coef[13] * t16 + t19 * coef[19] * t24 + t5 * coef[12] * t12 + t5 * coef[11] * t16);
        double t135 = Bm * Bm;
        double t136 = t135 * t135;
        double t138 = 0.1e1 / t136 / t135;
        double t143 = 0.1e1 / t136 / Bm;
        double t147 = 0.1e1 / t136;
        double t158 = 0.1e1 / t135 / Bm;
        complex<double> t257 = Nm[0] * (t138 * coefm[15] * t16 + t143 * coefm[8] * t24 + t147 * coefm[14] * t16 + t138 * coefm[0] * t12) + Nm[2] * (t158 * coefm[17] * t16 + t147 * coefm[7] * t24 + t143 * coefm[21] * t12 + t143 * coefm[26] * t16) + Nm[4] * (t158 * coefm[18] * t12 + t158 * coefm[16] * t16) + Nm[3] * (t158 * coefm[5] * t24 + t147 * coefm[19] * t12 + t147 * coefm[23] * t16) + Nm[5] * (t147 * coefm[9] * t16 + t143 * coefm[25] * t24) + Nm[10] * t158 * t16 * coefm[4] + Nm[7] * (t147 * coefm[24] * t16 + t143 * coefm[6] * t24) + Nm[9] * (t158 * coefm[3] * t16 + t147 * coefm[2] * t24) + Nm[11] * t22 * t158 * t16 * coefm[1] + Nm[6] * t147 * t16 * coefm[10] + Nm[1] * (t147 * coefm[13] * t16 + t143 * coefm[22] * t24 + t138 * coefm[12] * t12 + t138 * coefm[11] * t16) + Nm[8] * t147 * t16 * coefm[20];
        //
        X = t257 + t133;
        
    } // if (flag == 32)
    
    if (flag == 33)
    {
        double t2 = B * B;
        double t3 = t2 * t2;
        double t5 = 0.1e1 / t3 / B;
        double t8 = cos(Psi);
        double t12 = t8 * t8;
        double t13 = sin(Psi);
        double t15 = sin(theta);
        double t16 = t12 * t13 * t15;
        double t19 = 0.1e1 / t3 / t2;
        double t22 = t12 * t15;
        double t26 = t12 * t8;
        double t30 = t12 * t12;
        double t31 = t30 * t15;
        double t35 = cos(theta);
        double t36 = t12 * t35;
        double t37 = t36 * t13;
        double t39 = 0.1e1 / t3;
        double t45 = t8 * t13;
        double t52 = 0.1e1 / t2 / B;
        double t53 = t52 * t12;
        double t179 = t39 * t12;
        complex<double> t203 = N[0] * (t5 * coef[19] * t8 + t5 * coef[15] * t16 + t19 * coef[0] * t22 + t5 * coef[20] * t26 + t19 * coef[3] * t31 + t5 * coef[9] * t37 + t39 * coef[21] * t36 + t39 * coef[32] * t45) + N[10] * t15 * t53 * coef[30] + N[11] * t15 * t13 * t53 * coef[2] + N[9] * (t52 * coef[25] * t22 + t52 * coef[26] * t36 + t39 * coef[1] * t16 + t52 * coef[31] * t45) + N[1] * (t39 * coef[16] * t36 + t5 * coef[36] * t37 + t5 * coef[37] * t16 + t19 * coef[7] * t31 + t39 * coef[5] * t45 + t5 * coef[18] * t26 + t19 * coef[6] * t22 + t5 * coef[22] * t8) + N[7] * (t39 * coef[38] * t22 + t39 * coef[40] * t36 + t5 * coef[27] * t16 + t39 * coef[47] * t45) + N[2] * (t52 * coef[33] * t36 + t39 * coef[28] * t37 + t39 * coef[29] * t16 + t5 * coef[41] * t31 + t52 * coef[34] * t45 + t39 * coef[14] * t26 + t5 * coef[42] * t22 + t39 * coef[4] * t8) + N[4] * (t52 * coef[35] * t31 + t52 * coef[46] * t22) + N[3] * (t52 * coef[23] * t37 + t52 * coef[24] * t16 + t39 * coef[44] * t31 + t52 * coef[8] * t26 + t39 * coef[43] * t22 + t52 * coef[17] * t8) + N[8] * t15 * t179 * coef[45] + N[5] * (t39 * coef[10] * t22 + t39 * coef[12] * t36 + t5 * coef[39] * t16 + t39 * coef[11] * t45) + N[6] * t15 * t179 * coef[13];
        double t205 = Bm * Bm;
        double t206 = t205 * t205;
        double t208 = 0.1e1 / t206 / Bm;
        double t219 = 0.1e1 / t206 / t205;
        double t226 = 0.1e1 / t206;
        double t242 = 0.1e1 / t205 / Bm;
        double t243 = t242 * t12;
        double t332 = t226 * t12;
        complex<double> t392 = Nm[0] * (t208 * coefm[33] * t26 + t208 * coefm[34] * t8 + t208 * coefm[15] * t37 + t219 * coefm[35] * t31 + t219 * coefm[16] * t22 + t226 * coefm[28] * t45 + t208 * coefm[3] * t16 + t226 * coefm[17] * t36) + Nm[11] * t15 * t13 * t243 * coefm[0] + Nm[9] * (t242 * coefm[20] * t36 + t242 * coefm[26] * t22 + t226 * coefm[1] * t16 + t242 * coefm[19] * t45) + Nm[3] * (t242 * coefm[18] * t16 + t242 * coefm[21] * t37 + t226 * coefm[40] * t31 + t242 * coefm[7] * t26 + t226 * coefm[31] * t22 + t242 * coefm[8] * t8) + Nm[2] * (t242 * coefm[36] * t36 + t226 * coefm[24] * t37 + t226 * coefm[25] * t16 + t208 * coefm[47] * t31 + t242 * coefm[38] * t45 + t226 * coefm[11] * t26 + t208 * coefm[43] * t22 + t226 * coefm[5] * t8) + Nm[7] * (t226 * coefm[30] * t36 + t226 * coefm[45] * t22 + t208 * coefm[23] * t16 + t226 * coefm[46] * t45) + Nm[10] * t15 * t243 * coefm[22] + Nm[6] * t15 * t332 * coefm[6] + Nm[4] * (t242 * coefm[39] * t31 + t242 * coefm[37] * t22) + Nm[8] * t15 * t332 * coefm[42] + Nm[5] * (t226 * coefm[13] * t22 + t226 * coefm[14] * t36 + t208 * coefm[41] * t16 + t226 * coefm[10] * t45) + Nm[1] * (t226 * coefm[4] * t36 + t208 * coefm[32] * t37 + t208 * coefm[44] * t16 + t219 * coefm[9] * t31 + t226 * coefm[2] * t45 + t208 * coefm[29] * t26 + t219 * coefm[12] * t22 + t208 * coefm[27] * t8);
        //
        X = t392 + t203;
        
    } // if (flag == 33)
    
    // Final Output
    return X;
	
}



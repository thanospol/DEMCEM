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
#include "demcem_ss_va_nxrwg.h"
#include "demcem_constants.h"


using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> omega_functions
// ***********************************************************************


complex<double> omega_functions_ss_va_nxrwg (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[], int flag)


{                                 
	complex<double> X;
	//
    complex<double> j = Iunit;
    //
    
    if (flag == 11)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t7 = t3 * t5 * t4;
        complex<double> t8 = pow(GAMMA, 0.2e1);
        complex<double> t10 = 0.1e1 / t8 / GAMMA;
        complex<double> t11 = t7 * t10;
        complex<double> t12 = cos(theta_q);
        complex<double> t13 = sin(theta_p);
        complex<double> t15 = cos(theta_p);
        complex<double> t16 = coef[0];
        complex<double> t20 = sin(theta_q);
        complex<double> t22 = coef[1];
        complex<double> t26 = t15 * t15;
        complex<double> t28 = coef[2];
        complex<double> t32 = coef[3];
        complex<double> t36 = coef[4];
        complex<double> t39 = t12 * t36;
        complex<double> t43 = coef[5];
        complex<double> t46 = t20 * t43;
        complex<double> t52 = 0.1e1 / t8;
        complex<double> t54 = t7 * t52 * j;
        complex<double> t55 = ko * t12;
        complex<double> t56 = t13 * t15;
        complex<double> t60 = ko * t20;
        complex<double> t70 = t7 * t52;
        complex<double> t71 = j * ko;
        //
        X = K[2] * (-t11 * t12 * t13 * t15 * t16 - t11 * t20 * t13 * t15 * t22 - t11 * t12 * t26 * t28 - t11 * t20 * t26 * t32 - t7 * t10 * t12 * t36 + t11 * t39 * t26 - t7 * t10 * t20 * t43 + t11 * t46 * t26) + K[3] * (-t54 * t55 * t56 * t16 - t54 * t60 * t56 * t22 - t54 * t55 * t26 * t28 - t54 * t60 * t26 * t32 - t70 * t71 * t39 + t54 * t55 * t36 * t26 - t70 * t71 * t46 + t54 * t60 * t43 * t26);
    } // if (flag == 11)
    
    if (flag == 12)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t4 * t3;
        complex<double> t6 = t2 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = cos(theta_p);
        complex<double> t11 = t10 * t10;
        complex<double> t13 = coef[4];
        complex<double> t16 = t2 * t2;
        complex<double> t17 = t16 * t4;
        complex<double> t18 = t17 * t9;
        complex<double> t19 = sin(theta_p);
        complex<double> t20 = sin(theta_q);
        complex<double> t21 = t19 * t20;
        complex<double> t22 = coef[1];
        complex<double> t25 = coef[0];
        complex<double> t26 = t9 * t25;
        complex<double> t31 = t19 * t10;
        complex<double> t32 = coef[5];
        complex<double> t36 = coef[10];
        complex<double> t42 = t16 * t5;
        complex<double> t43 = t42 * t9;
        complex<double> t44 = cos(theta_q);
        complex<double> t46 = coef[6];
        complex<double> t50 = coef[2];
        complex<double> t54 = 0.1e1 / t7;
        complex<double> t55 = t54 * j;
        complex<double> t56 = t17 * t55;
        complex<double> t57 = ko * t19;
        complex<double> t61 = coef[9];
        complex<double> t62 = t20 * t61;
        complex<double> t65 = t6 * t54;
        complex<double> t66 = j * ko;
        complex<double> t77 = coef[8];
        complex<double> t87 = coef[3];
        complex<double> t91 = t44 * t77;
        complex<double> t95 = coef[7];
        complex<double> t101 = -t43 * t44 * t11 * t46 - t43 * t44 * t19 * t10 * t50 - t56 * t57 * t20 * t22 + t43 * t62 * t11 - t65 * t66 * t25 - t65 * t66 * t11 * t13 - t6 * t55 * t57 * t10 * t32 - t42 * t9 * t44 * t77 - t42 * t9 * t20 * t61 - t56 * ko * t10 * t20 * t36 - t43 * t21 * t10 * t87 + t43 * t91 * t11 - t43 * t20 * t11 * t95 + t65 * t66 * t25 * t11;
        complex<double> t104 = t42 * t55;
        complex<double> t105 = ko * t20;
        complex<double> t109 = ko * t44;
        complex<double> t122 = t42 * t54;
        //
        X = K[1] * (-t6 * t9 * t11 * t13 - t18 * t21 * t22 - t6 * t26 + t6 * t26 * t11 - t6 * t9 * t31 * t32 - t18 * t10 * t20 * t36) + K[2] * t101 + K[3] * (-t104 * t105 * t31 * t87 - t104 * t109 * t11 * t46 - t104 * t109 * t31 * t50 + t104 * t109 * t77 * t11 + t104 * t105 * t61 * t11 - t122 * t66 * t62 - t122 * t66 * t91 - t104 * t105 * t11 * t95);
        
    }// if (flag == 12)
    
    if (flag == 13)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t4 * t3;
        complex<double> t6 = t2 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = cos(theta_p);
        complex<double> t11 = t10 * t10;
        complex<double> t13 = coef[10];
        complex<double> t17 = sin(theta_p);
        complex<double> t18 = t17 * t10;
        complex<double> t19 = coef[12];
        complex<double> t22 = t2 * t2;
        complex<double> t23 = t22 * t4;
        complex<double> t24 = t23 * t9;
        complex<double> t25 = cos(theta_q);
        complex<double> t27 = coef[7];
        complex<double> t30 = sin(theta_q);
        complex<double> t32 = coef[4];
        complex<double> t35 = t30 * t10;
        complex<double> t36 = coef[6];
        complex<double> t39 = t25 * t10;
        complex<double> t40 = coef[5];
        complex<double> t43 = coef[11];
        complex<double> t44 = t9 * t43;
        complex<double> t51 = t22 * t5;
        complex<double> t52 = t51 * t9;
        complex<double> t54 = coef[1];
        complex<double> t57 = coef[8];
        complex<double> t61 = coef[9];
        complex<double> t65 = coef[2];
        complex<double> t66 = t30 * t65;
        complex<double> t69 = coef[3];
        complex<double> t70 = t25 * t69;
        complex<double> t73 = 0.1e1 / t7;
        complex<double> t74 = t73 * j;
        complex<double> t75 = t23 * t74;
        complex<double> t76 = ko * t30;
        complex<double> t80 = ko * t25;
        complex<double> t93 = t6 * t73;
        complex<double> t94 = j * ko;
        complex<double> t98 = coef[0];
        complex<double> t102 = ko * t17;
        complex<double> t115 = -t52 * t11 * t30 * t54 - t52 * t18 * t30 * t57 - t52 * t39 * t17 * t61 + t52 * t66 * t11 + t52 * t70 * t11 - t75 * t76 * t17 * t32 - t75 * t80 * t10 * t40 - t75 * t76 * t10 * t36 - t75 * t80 * t17 * t27 - t51 * t9 * t30 * t65 - t93 * t94 * t43 - t52 * t25 * t11 * t98 - t6 * t74 * t102 * t10 * t19 - t51 * t9 * t25 * t69 - t93 * t94 * t11 * t13 + t93 * t94 * t43 * t11;
        complex<double> t118 = t51 * t74;
        complex<double> t135 = t51 * t73;
        //
        X = K[1] * (-t6 * t9 * t11 * t13 - t6 * t9 * t18 * t19 - t24 * t25 * t17 * t27 - t24 * t30 * t17 * t32 - t24 * t35 * t36 - t24 * t39 * t40 + t6 * t44 * t11 - t6 * t44) + K[2] * t115 + K[3] * (-t118 * t102 * t35 * t57 - t118 * t80 * t11 * t98 - t118 * ko * t11 * t30 * t54 + t118 * t76 * t65 * t11 + t118 * t80 * t69 * t11 - t135 * t94 * t70 - t135 * t94 * t66 - t118 * t80 * t18 * t61);
        
    }// if (flag == 13)
    
    if (flag == 21)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = t6 * t9;
        complex<double> t11 = sin(theta_p);
        complex<double> t12 = sin(theta_q);
        complex<double> t13 = t11 * t12;
        complex<double> t14 = coef[0];
        complex<double> t17 = cos(theta_q);
        complex<double> t18 = t17 * t11;
        complex<double> t19 = coef[7];
        complex<double> t22 = cos(theta_p);
        complex<double> t24 = coef[8];
        complex<double> t28 = coef[9];
        complex<double> t35 = t3 * t5 * t4;
        complex<double> t36 = t35 * t9;
        complex<double> t37 = t22 * t22;
        complex<double> t39 = coef[4];
        complex<double> t42 = coef[1];
        complex<double> t46 = coef[2];
        complex<double> t51 = coef[3];
        complex<double> t54 = coef[6];
        complex<double> t55 = t12 * t54;
        complex<double> t59 = coef[5];
        complex<double> t62 = t17 * t59;
        complex<double> t68 = 0.1e1 / t7;
        complex<double> t69 = t68 * j;
        complex<double> t70 = t6 * t69;
        complex<double> t75 = ko * t17;
        complex<double> t82 = ko * t12;
        complex<double> t86 = -t36 * t12 * t37 * t39 - t36 * t18 * t22 * t42 - t36 * t13 * t22 * t46 - t36 * t17 * t37 * t51 + t36 * t55 * t37 - t35 * t9 * t17 * t59 + t36 * t62 * t37 - t35 * t9 * t12 * t54 - t70 * ko * t11 * t12 * t14 - t70 * t75 * t11 * t19 - t70 * t75 * t22 * t24 - t70 * t82 * t22 * t28;
        complex<double> t89 = t35 * t69;
        complex<double> t93 = t11 * t22;
        complex<double> t106 = t35 * t68;
        complex<double> t107 = j * ko;
        //
        X = K[1] * (-t10 * t13 * t14 - t10 * t18 * t19 - t10 * t17 * t22 * t24 - t10 * t12 * t22 * t28) + K[2] * t86 + K[3] * (-t89 * t82 * t37 * t39 - t89 * t75 * t93 * t42 - t89 * t82 * t93 * t46 - t89 * t75 * t37 * t51 + t89 * t82 * t54 * t37 - t106 * t107 * t62 + t89 * t75 * t59 * t37 - t106 * t107 * t55);
        
    }// if (flag == 21)
    
    if (flag == 22)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t4 * t3;
        complex<double> t6 = t2 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = coef[7];
        complex<double> t11 = t9 * t10;
        complex<double> t13 = t2 * t2;
        complex<double> t14 = t13 * t4;
        complex<double> t15 = t14 * t9;
        complex<double> t16 = sin(theta_p);
        complex<double> t17 = sin(theta_q);
        complex<double> t18 = t16 * t17;
        complex<double> t19 = coef[0];
        complex<double> t22 = cos(theta_p);
        complex<double> t23 = t22 * t22;
        complex<double> t25 = coef[6];
        complex<double> t28 = t13 * t3;
        complex<double> t29 = 0.1e1 / t7;
        complex<double> t31 = j * ko;
        complex<double> t32 = coef[4];
        complex<double> t37 = coef[14];
        complex<double> t42 = cos(theta_q);
        complex<double> t43 = t42 * t16;
        complex<double> t44 = coef[13];
        complex<double> t47 = t2 * t4;
        complex<double> t48 = t47 * t29;
        complex<double> t49 = coef[5];
        complex<double> t54 = coef[15];
        complex<double> t57 = coef[3];
        complex<double> t62 = t16 * t22;
        complex<double> t63 = coef[8];
        complex<double> t66 = -t6 * t11 - t15 * t18 * t19 - t6 * t9 * t23 * t25 - t28 * t29 * t31 * t17 * t32 - t15 * t22 * t17 * t37 + t6 * t11 * t23 - t15 * t43 * t44 - t48 * t31 * t22 * t49 - t15 * t42 * t22 * t54 - t48 * t31 * t16 * t57 - t6 * t9 * t62 * t63;
        complex<double> t69 = t13 * t5;
        complex<double> t70 = t69 * t9;
        complex<double> t71 = coef[12];
        complex<double> t72 = t17 * t71;
        complex<double> t75 = t6 * t29;
        complex<double> t78 = t29 * j;
        complex<double> t80 = ko * t16;
        complex<double> t84 = coef[1];
        complex<double> t88 = coef[2];
        complex<double> t98 = t14 * t78;
        complex<double> t102 = t9 * t17;
        complex<double> t106 = coef[9];
        complex<double> t110 = coef[11];
        complex<double> t113 = ko * t42;
        complex<double> t125 = coef[10];
        complex<double> t128 = t42 * t110;
        complex<double> t131 = t70 * t72 * t23 - t75 * t31 * t10 - t6 * t78 * t80 * t22 * t63 - t70 * t43 * t22 * t84 - t70 * t18 * t22 * t88 + t75 * t31 * t10 * t23 - t75 * t31 * t23 * t25 - t98 * t80 * t17 * t19 - t69 * t102 * t71 - t70 * t42 * t23 * t106 - t69 * t9 * t42 * t110 - t98 * t113 * t16 * t44 - t98 * ko * t22 * t17 * t37 - t98 * t113 * t22 * t54 - t70 * t17 * t23 * t125 + t70 * t128 * t23;
        complex<double> t134 = t69 * t78;
        complex<double> t135 = ko * t17;
        complex<double> t154 = t69 * t29;
        //
        X = K[1] * t66 + K[2] * t131 + K[3] * (t134 * t135 * t71 * t23 + t134 * t113 * t110 * t23 - t134 * t135 * t62 * t88 - t134 * t113 * t23 * t106 - t134 * t135 * t23 * t125 - t134 * t113 * t62 * t84 - t154 * t31 * t72 - t154 * t31 * t128) + K[0] * (-t28 * t102 * t32 - t47 * t9 * t22 * t49 - t47 * t9 * t16 * t57);
        
    }// if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t4 * t3;
        complex<double> t6 = t2 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = coef[16];
        complex<double> t11 = t9 * t10;
        complex<double> t13 = cos(theta_p);
        complex<double> t14 = t13 * t13;
        complex<double> t17 = t2 * t2;
        complex<double> t18 = t17 * t4;
        complex<double> t19 = t18 * t9;
        complex<double> t20 = cos(theta_q);
        complex<double> t22 = coef[15];
        complex<double> t26 = coef[3];
        complex<double> t29 = t17 * t3;
        complex<double> t30 = 0.1e1 / t7;
        complex<double> t31 = t29 * t30;
        complex<double> t32 = j * ko;
        complex<double> t33 = coef[4];
        complex<double> t37 = sin(theta_p);
        complex<double> t38 = sin(theta_q);
        complex<double> t39 = t37 * t38;
        complex<double> t40 = coef[0];
        complex<double> t43 = t2 * t4;
        complex<double> t44 = t43 * t30;
        complex<double> t45 = coef[7];
        complex<double> t50 = coef[14];
        complex<double> t53 = coef[5];
        complex<double> t57 = coef[6];
        complex<double> t62 = t37 * t13;
        complex<double> t63 = coef[8];
        complex<double> t66 = t20 * t37;
        complex<double> t67 = coef[13];
        complex<double> t70 = -t6 * t11 + t6 * t11 * t14 - t19 * t20 * t13 * t22 - t6 * t9 * t14 * t26 - t31 * t32 * t20 * t33 - t19 * t39 * t40 - t44 * t32 * t13 * t45 - t19 * t38 * t13 * t50 - t31 * t32 * t38 * t53 - t44 * t32 * t37 * t57 - t6 * t9 * t62 * t63 - t19 * t66 * t67;
        complex<double> t79 = t9 * t20;
        complex<double> t82 = t9 * t38;
        complex<double> t88 = t6 * t30;
        complex<double> t92 = t17 * t5;
        complex<double> t93 = t92 * t9;
        complex<double> t94 = coef[2];
        complex<double> t98 = coef[12];
        complex<double> t99 = t38 * t98;
        complex<double> t105 = coef[11];
        complex<double> t108 = t30 * j;
        complex<double> t109 = t18 * t108;
        complex<double> t110 = ko * t20;
        complex<double> t114 = coef[1];
        complex<double> t123 = coef[9];
        complex<double> t127 = ko * t37;
        complex<double> t131 = t20 * t105;
        complex<double> t137 = ko * t38;
        complex<double> t142 = coef[10];
        complex<double> t148 = t88 * t32 * t10 * t14 - t93 * t39 * t13 * t94 + t93 * t99 * t14 - t88 * t32 * t14 * t26 - t92 * t79 * t105 - t109 * t110 * t13 * t22 - t93 * t66 * t13 * t114 - t92 * t82 * t98 - t88 * t32 * t10 - t93 * t20 * t14 * t123 - t6 * t108 * t127 * t13 * t63 + t93 * t131 * t14 - t109 * t110 * t37 * t67 - t109 * t137 * t13 * t50 - t93 * t38 * t14 * t142 - t109 * t127 * t38 * t40;
        complex<double> t151 = t92 * t30;
        complex<double> t154 = t92 * t108;
        //
        X = K[1] * t70 + K[0] * (-t43 * t9 * t37 * t57 - t43 * t9 * t13 * t45 - t29 * t79 * t33 - t29 * t82 * t53) + K[2] * t148 + K[3] * (-t151 * t32 * t99 + t154 * t110 * t105 * t14 + t154 * t137 * t98 * t14 - t154 * t137 * t14 * t142 - t154 * t110 * t14 * t123 - t154 * t110 * t62 * t114 - t154 * t137 * t62 * t94 - t151 * t32 * t131);
        
    }// if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t3 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = t6 * t9;
        complex<double> t11 = sin(theta_p);
        complex<double> t12 = sin(theta_q);
        complex<double> t13 = t11 * t12;
        complex<double> t14 = coef[0];
        complex<double> t17 = cos(theta_q);
        complex<double> t18 = cos(theta_p);
        complex<double> t20 = coef[7];
        complex<double> t24 = coef[8];
        complex<double> t27 = t17 * t11;
        complex<double> t28 = coef[9];
        complex<double> t35 = t3 * t5 * t4;
        complex<double> t36 = t35 * t9;
        complex<double> t37 = t18 * t18;
        complex<double> t39 = coef[4];
        complex<double> t42 = coef[1];
        complex<double> t46 = coef[2];
        complex<double> t51 = coef[3];
        complex<double> t54 = coef[6];
        complex<double> t55 = t12 * t54;
        complex<double> t59 = coef[5];
        complex<double> t62 = t17 * t59;
        complex<double> t68 = 0.1e1 / t7;
        complex<double> t69 = t68 * j;
        complex<double> t70 = t6 * t69;
        complex<double> t75 = ko * t17;
        complex<double> t79 = ko * t12;
        complex<double> t86 = -t36 * t12 * t37 * t39 - t36 * t27 * t18 * t42 - t36 * t13 * t18 * t46 - t36 * t17 * t37 * t51 + t36 * t55 * t37 - t35 * t9 * t17 * t59 + t36 * t62 * t37 - t35 * t9 * t12 * t54 - t70 * ko * t11 * t12 * t14 - t70 * t75 * t18 * t20 - t70 * t79 * t18 * t24 - t70 * t75 * t11 * t28;
        complex<double> t89 = t35 * t69;
        complex<double> t93 = t11 * t18;
        complex<double> t106 = t35 * t68;
        complex<double> t107 = j * ko;
        //
        X = K[1] * (-t10 * t13 * t14 - t10 * t17 * t18 * t20 - t10 * t12 * t18 * t24 - t10 * t27 * t28) + K[2] * t86 + K[3] * (-t89 * t79 * t37 * t39 - t89 * t75 * t93 * t42 - t89 * t79 * t93 * t46 - t89 * t75 * t37 * t51 + t89 * t79 * t54 * t37 - t106 * t107 * t62 + t89 * t75 * t59 * t37 - t106 * t107 * t55);
        
    }// if (flag == 31)
    
    if (flag == 32)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = cos(Psi);
        complex<double> t4 = t3 * t3;
        complex<double> t5 = t4 * t3;
        complex<double> t6 = t2 * t5;
        complex<double> t7 = pow(GAMMA, 0.2e1);
        complex<double> t9 = 0.1e1 / t7 / GAMMA;
        complex<double> t10 = coef[7];
        complex<double> t11 = t9 * t10;
        complex<double> t12 = cos(theta_p);
        complex<double> t13 = t12 * t12;
        complex<double> t16 = t2 * t4;
        complex<double> t17 = 0.1e1 / t7;
        complex<double> t18 = t16 * t17;
        complex<double> t19 = j * ko;
        complex<double> t20 = sin(theta_p);
        complex<double> t21 = coef[3];
        complex<double> t26 = coef[6];
        complex<double> t31 = t20 * t12;
        complex<double> t32 = coef[8];
        complex<double> t35 = t2 * t2;
        complex<double> t36 = t35 * t4;
        complex<double> t37 = t36 * t9;
        complex<double> t38 = sin(theta_q);
        complex<double> t40 = coef[14];
        complex<double> t43 = cos(theta_q);
        complex<double> t45 = coef[13];
        complex<double> t48 = coef[4];
        complex<double> t52 = t20 * t38;
        complex<double> t53 = coef[0];
        complex<double> t56 = t35 * t3;
        complex<double> t58 = coef[5];
        complex<double> t62 = t43 * t20;
        complex<double> t63 = coef[15];
        complex<double> t66 = t6 * t11 * t13 - t18 * t19 * t20 * t21 - t6 * t9 * t13 * t26 - t6 * t11 - t6 * t9 * t31 * t32 - t37 * t38 * t12 * t40 - t37 * t43 * t12 * t45 - t18 * t19 * t12 * t48 - t37 * t52 * t53 - t56 * t17 * t19 * t38 * t58 - t37 * t62 * t63;
        complex<double> t69 = t6 * t17;
        complex<double> t73 = t35 * t5;
        complex<double> t74 = t73 * t9;
        complex<double> t76 = coef[9];
        complex<double> t79 = coef[12];
        complex<double> t80 = t38 * t79;
        complex<double> t84 = coef[11];
        complex<double> t87 = t17 * j;
        complex<double> t89 = ko * t20;
        complex<double> t93 = coef[1];
        complex<double> t97 = coef[2];
        complex<double> t101 = t36 * t87;
        complex<double> t102 = ko * t43;
        complex<double> t110 = coef[10];
        complex<double> t113 = t43 * t84;
        complex<double> t121 = ko * t38;
        complex<double> t128 = t9 * t38;
        complex<double> t131 = t69 * t19 * t10 * t13 - t74 * t43 * t13 * t76 + t74 * t80 * t13 - t73 * t9 * t43 * t84 - t6 * t87 * t89 * t12 * t32 - t74 * t62 * t12 * t93 - t74 * t52 * t12 * t97 - t101 * t102 * t20 * t63 - t101 * t102 * t12 * t45 - t74 * t38 * t13 * t110 + t74 * t113 * t13 - t69 * t19 * t10 - t101 * t89 * t38 * t53 - t101 * t121 * t12 * t40 - t69 * t19 * t13 * t26 - t73 * t128 * t79;
        complex<double> t134 = t73 * t87;
        complex<double> t153 = t73 * t17;
        //
        X = K[1] * t66 + K[2] * t131 + K[3] * (t134 * t121 * t79 * t13 + t134 * t102 * t84 * t13 - t134 * t121 * t31 * t97 - t134 * t102 * t13 * t76 - t134 * t121 * t13 * t110 - t134 * t102 * t31 * t93 - t153 * t19 * t80 - t153 * t19 * t113) + K[0] * (-t16 * t9 * t12 * t48 - t56 * t128 * t58 - t16 * t9 * t20 * t21);
        
    }// if (flag == 32)
    
    if (flag == 33)
    {
        complex<double> t2 = sin(Psi);
        complex<double> t3 = t2 * t2;
        complex<double> t4 = cos(Psi);
        complex<double> t5 = t4 * t4;
        complex<double> t6 = t5 * t4;
        complex<double> t7 = t3 * t6;
        complex<double> t8 = pow(GAMMA, 0.2e1);
        complex<double> t10 = 0.1e1 / t8 / GAMMA;
        complex<double> t11 = t7 * t10;
        complex<double> t12 = sin(theta_q);
        complex<double> t13 = cos(theta_p);
        complex<double> t14 = t13 * t13;
        complex<double> t16 = coef[10];
        complex<double> t19 = t3 * t5;
        complex<double> t20 = 0.1e1 / t8;
        complex<double> t21 = t20 * j;
        complex<double> t22 = t19 * t21;
        complex<double> t23 = cos(theta_q);
        complex<double> t24 = ko * t23;
        complex<double> t25 = coef[13];
        complex<double> t29 = ko * t12;
        complex<double> t30 = coef[14];
        complex<double> t34 = t2 * t6;
        complex<double> t35 = t34 * t20;
        complex<double> t36 = j * ko;
        complex<double> t37 = coef[2];
        complex<double> t41 = sin(theta_p);
        complex<double> t42 = ko * t41;
        complex<double> t43 = coef[3];
        complex<double> t47 = coef[6];
        complex<double> t51 = t10 * t23;
        complex<double> t52 = coef[11];
        complex<double> t55 = t10 * t12;
        complex<double> t56 = coef[12];
        complex<double> t60 = coef[16];
        complex<double> t65 = coef[9];
        complex<double> t68 = t12 * t56;
        complex<double> t73 = t23 * t41;
        complex<double> t74 = coef[4];
        complex<double> t78 = t23 * t52;
        complex<double> t81 = coef[15];
        complex<double> t85 = t12 * t41;
        complex<double> t86 = coef[5];
        complex<double> t90 = -t11 * t12 * t14 * t16 - t22 * t24 * t13 * t25 - t22 * t29 * t13 * t30 - t35 * t36 * t14 * t37 - t22 * t42 * t12 * t43 + t35 * t36 * t47 * t14 - t7 * t51 * t52 - t7 * t55 * t56 - t34 * t21 * t42 * t13 * t60 - t11 * t23 * t14 * t65 + t11 * t68 * t14 - t35 * t36 * t47 - t11 * t73 * t13 * t74 + t11 * t78 * t14 - t22 * t24 * t41 * t81 - t11 * t85 * t13 * t86;
        complex<double> t93 = t7 * t21;
        complex<double> t94 = t41 * t13;
        complex<double> t104 = t7 * t20;
        complex<double> t121 = t2 * t5;
        complex<double> t123 = coef[1];
        complex<double> t127 = coef[0];
        complex<double> t130 = t3 * t4;
        complex<double> t131 = coef[7];
        complex<double> t134 = coef[8];
        complex<double> t140 = t10 * t47;
        complex<double> t143 = t19 * t10;
        complex<double> t147 = t121 * t20;
        complex<double> t154 = t130 * t20;
        complex<double> t175 = t34 * t140 * t14 - t143 * t23 * t13 * t25 - t147 * t36 * t41 * t127 - t147 * t36 * t13 * t123 - t154 * t36 * t12 * t134 - t154 * t36 * t23 * t131 - t34 * t140 - t143 * t85 * t43 - t143 * t12 * t13 * t30 - t143 * t73 * t81 - t34 * t10 * t94 * t60 - t34 * t10 * t14 * t37;
        //
        X = K[2] * t90 + K[3] * (-t93 * t24 * t94 * t74 - t93 * t29 * t94 * t86 + t93 * t24 * t52 * t14 - t104 * t36 * t68 - t93 * t24 * t14 * t65 - t104 * t36 * t78 - t93 * t29 * t14 * t16 + t93 * t29 * t56 * t14) + K[0] * (-t121 * t10 * t13 * t123 - t121 * t10 * t41 * t127 - t130 * t51 * t131 - t130 * t55 * t134) + K[1] * t175;
        
    }// if (flag == 33)
    
	// Final Output
	return X;
}



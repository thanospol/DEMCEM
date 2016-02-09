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
#include "demcem_constants.h"

using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_ss_ea
// ***********************************************************************

void coefficients_ss_ea_rwg (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[], int flag )
{
    double aa[3], bb[3], cc[3], dd[3];
    
    for (int i = 0; i < 3; i++)
    {
        aa[i] = r1[i] - r1[i];
        bb[i] = r2[i] - r1[i];
        cc[i] = r3[i] - r1[i];
        dd[i] = r4[i] - r1[i];
    }
    //
    complex<double> j   = Iunit;
    //
    if (flag == 11)
    {
        complex<double> c1  = cc[0] * dd[1] * bb[2] / 0.3e1 - cc[0] * dd[2] * bb[1] / 0.3e1 + cc[1] * dd[2] * bb[0] / 0.3e1
        - cc[1] * dd[0] * bb[2] / 0.3e1 + cc[2] * dd[0] * bb[1] / 0.3e1 - cc[2] * dd[1] * bb[0] / 0.3e1;
        //
        complex<double> t6  = pow(j,2);
        complex<double> t9  = pow(ko,2);
        complex<double> t12 = 3.0 * c1 / t6 / t9;
        //
        coef[0] = 3.0 * c1 / j / ko;
        coef[1] = c1;
        coef[2] = -t12;
        coef[3] = t12;
        //
        complex<double> cm1 = cc[0] * bb[2] * dd[1] / 0.3e1 - cc[0] * bb[1] * dd[2] / 0.3e1 + cc[1] * bb[0] * dd[2] / 0.3e1
        - cc[1] * bb[2] * dd[0] / 0.3e1 + cc[2] * bb[1] * dd[0] / 0.3e1 - cc[2] * bb[0] * dd[1] / 0.3e1;
        //
        complex<double> t1 = pow(j,2);
        complex<double> t4 = pow(ko,2);
        complex<double> t7 = 3.0 * cm1 / t1 / t4;
        //
        coefm[0] = -t7;
        coefm[1] = t7;
        coefm[2] = 3.0 * cm1 / j / ko;
        coefm[3] = cm1;
	
    } // if (flag == 11)
    
    
    
    if (flag == 13)
    {
        complex<double> c[4], cm[4];
        //
        complex<double> t1  = cc[0] * bb[1];
        complex<double> t2  = sqrt(0.3e1);
        complex<double> t3  = dd[2] * t2;
        complex<double> t5  = cc[0] * bb[2];
        complex<double> t6  = dd[1] * t2;
        complex<double> t8  = cc[1] * bb[2];
        complex<double> t9  = dd[0] * t2;
        complex<double> t11 = cc[1] * bb[0];
        complex<double> t13 = cc[2] * bb[0];
        complex<double> t15 = cc[2] * bb[1];
        complex<double> t17 = -t1 * t3 + t5 * t6 - t8 * t9 + t11 * t3 - t13 * t6 + t15 * t9;
        //
        c[0] = t17 / 0.6e1;
        c[1] = -t8 * dd[0] / 0.6e1 + t11 * dd[2] / 0.6e1 - t1 * dd[2] / 0.6e1 - t13 * dd[1] / 0.6e1 + t5 * dd[1] / 0.6e1 + t15 * dd[0] / 0.6e1;
        c[2] = -t17 / 0.6e1;
        c[3] = -t17 / 0.6e1;
        //
        t3  = 0.1e1 / j / ko;
        t5  = 2.0 * t3 * c[2];
        t11 = 2.0 * t3 * c[3];
        complex<double> t12 = pow(ko,2);
        complex<double> t14 = pow(j,2);
        complex<double> t16 = 0.1e1 / t12 / t14;
        complex<double> t18 = 3.0 * t16 * c[0];
        complex<double> t20 = 3.0 * t16 * c[1];
        //
        coef[0]  = -t5;
        coef[1]  = c[1];
        coef[2]  = c[3];
        coef[3]  = c[0];
        coef[4]  = 3.0 * t3 * c[0];
        coef[5]  = 3.0 * t3 * c[1];
        coef[6]  = -t11;
        coef[7]  = t5;
        coef[8]  = -t18;
        coef[9]  = -t20;
        coef[10] = c[2];
        coef[11] = t11;
        coef[12] = t18;
        coef[13] = t20;
        //
        t1  = cc[0] * bb[1];
        t2  = sqrt(0.3e1);
        t3  = dd[2] * t2;
        t5  = cc[2] * bb[2];
        t6  = dd[1] * t2;
        t8  = cc[1] * bb[2];
        t9  = dd[0] * t2;
        t11 = cc[1] * bb[0];
        t13 = cc[2] * bb[0];
        t15 = cc[2] * bb[1];
        t17 = -t1 * t3 + t5 * t6 - t8 * t9 + t11 * t3 - t13 * t6 + t15 * t9;
        //
        cm[0] = t17 / 0.6e1;
        cm[1] = -t17 / 0.6e1;
        cm[2] = -t17 / 0.6e1;
        cm[3] = -t8 * dd[0] / 0.6e1 + t11 * dd[2] / 0.6e1 - t1 * dd[2] / 0.6e1 - t13 * dd[1] / 0.6e1 + t5 * dd[1] / 0.6e1 + t15 * dd[0] / 0.6e1;;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * cm[1];
        complex<double> t7 = 2.0 * t3 * cm[0];
        t8  = ko * ko;
        complex<double> t10 = j * j;
        t12 = 0.1e1 / t8 / t10;
        t14 = 3.0 * t12 * cm[2];
        t16 = 3.0 * t12 * cm[3];
        //
        coefm[0]  = -t5;
        coefm[1]  = -t7;
        coefm[2]  = -t14;
        coefm[3]  = -t16;
        coefm[4]  = t5;
        coefm[5]  = cm[0];
        coefm[6]  = 3.0 * t3 * cm[2];
        coefm[7]  = 3.0 * t3 * cm[3];
        coefm[8]  = cm[2];
        coefm[9]  = cm[3];
        coefm[10] = t7;
        coefm[11] = t14;
        coefm[12] = t16;
        coefm[13] = cm[1];
        
    } // if (flag == 13)
    
    
    if (flag == 22)
    {
        complex<double> c0  = -bb[0] * cc[1] * dd[2] / 0.3e1 + bb[0] * cc[2] * dd[1] / 0.3e1 - bb[1] * cc[2] * dd[0] / 0.3e1 + bb[1] * cc[0] * dd[2] / 0.3e1 - bb[2] * cc[0] * dd[1] / 0.3e1 + bb[2] * cc[1] * dd[0] / 0.3e1;
        complex<double> t1  = pow(j,2);
        complex<double> t4  = pow(ko,2);
        complex<double> t7  = 3.0 * c0 / t1 / t4;
        //
        coef[0] = c0;
        coef[1] = -t7;
        coef[2] = 3.0 * c0 / j / ko;
        coef[3] = t7;
        //
        complex<double> cm0 = -bb[0] * cc[1] * dd[2] / 0.3e1 + bb[0] * cc[2] * dd[1] / 0.3e1 - bb[1] * cc[2] * dd[0] / 0.3e1 + bb[1] * cc[0] * dd[2] / 0.3e1 - bb[2] * cc[0] * dd[1] / 0.3e1 + bb[2] * cc[1] * dd[0] / 0.3e1;
        complex<double> t6  = pow(j,2);
        complex<double> t9  = pow(ko,2);
        complex<double> t12  = 3.0 * cm0 / t6 / t9;
        //
        coefm[0] = cm0;
        coefm[1] = 3.0 * cm0 / j / ko;
        coefm[2] = -t12;
        coefm[3] = t12;
        
    } // if (flag == 22)
    
    if (flag == 23)
    {
        complex<double> c[4], cm[4];
        //
        complex<double> t1 = cc[0] * bb[2];
        complex<double> t2 = sqrt(0.3e1);
        complex<double> t3 = dd[1] * t2;
        complex<double> t5 = cc[1] * bb[2];
        complex<double> t6 = dd[0] * t2;
        complex<double> t8 = cc[0] * bb[1];
        complex<double> t9 = dd[2] * t2;
        complex<double> t11 = cc[2] * bb[1];
        complex<double> t13 = cc[2] * bb[0];
        complex<double> t15 = cc[1] * bb[0];
        complex<double> t17 = t1 * t3 - t5 * t6 - t8 * t9 + t11 * t6 - t13 * t3 + t15 * t9;
        //
        c[0] = t17 / 0.6e1;
        c[1] = t17 / 0.6e1;
        c[2] = t5 * dd[0] / 0.6e1 - t11 * dd[0] / 0.6e1 + t13 * dd[1] / 0.6e1 - t1 * dd[1] / 0.6e1 + t8 * dd[2] / 0.6e1 - t15 * dd[2] / 0.6e1;
        c[3] = -t17 / 0.6e1;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * c[0];
        complex<double> t7 = 2.0 * t3 * c[3];
        t8 = j * j;
        complex<double> t10 = ko * ko;
        complex<double> t12 = 0.1e1 / t8 / t10;
        complex<double> t14 = 3.0 * t12 * c[2];
        complex<double> t16 = 3.0 * t12 * c[1];
        //
        coef[0] = -t5;
        coef[1] = -t7;
        coef[2] = -t14;
        coef[3] = -t16;
        coef[4] = t5;
        coef[5] = c[3];
        coef[6] = c[2];
        coef[7] = c[1];
        coef[8] = 3.0 * t3 * c[2];
        coef[9] = 3.0 * t3 * c[1];
        coef[10] = c[0];
        coef[11] = t7;
        coef[12] = t14;
        coef[13] = t16;
        //
        t1 = cc[0] * bb[2];
        t2 = sqrt(0.3e1);
        t3 = dd[1] * t2;
        t5 = cc[1] * bb[2];
        t6 = dd[0] * t2;
        t8 = cc[0] * bb[1];
        t9 = dd[2] * t2;
        t11 = cc[2] * bb[1];
        t13 = cc[2] * bb[0];
        t15 = cc[1] * bb[0];
        t17 = t1 * t3 - t5 * t6 - t8 * t9 + t11 * t6 - t13 * t3 + t15 * t9;
        //
        cm[0] = t17 / 0.6e1;
        cm[1] = -t17 / 0.6e1;
        cm[2] = t17 / 0.6e1;
        cm[3] = t5 * dd[0] / 0.6e1 - t11 * dd[0] / 0.6e1 + t13 * dd[1] / 0.6e1 - t1 * dd[1] / 0.6e1 + t8 * dd[2] / 0.6e1 - t15 * dd[2] / 0.6e1;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * cm[0];
        t6 = j * j;
        t8 = ko * ko;
        t10 = 0.1e1 / t6 / t8;
        t12 = 3.0 * t10 * cm[1];
        t14 = 3.0 * t10 * cm[3];
        t16 = 2.0 * t3 * cm[2];
        //
        coefm[0]  = t5;
        coefm[1]  = -t12;
        coefm[2]  = -t14;
        coefm[3]  = -t16;
        coefm[4]  = t16;
        coefm[5]  = t14;
        coefm[6]  = t12;
        coefm[7]  = cm[0];
        coefm[8]  = -t5;
        coefm[9]  = 3.0 * t3 * cm[3];
        coefm[10] = 3.0 * t3 * cm[1];
        coefm[11] = cm[1];
        coefm[12] = cm[3];
        coefm[13] = cm[2];
        
    } // if (flag == 23)
    
    if (flag == 31)
    {
        complex<double> t1  = cc[1] * dd[2];
        complex<double> t2  = sqrt(3.0);
        complex<double> t3  = bb[0] * t2;
        complex<double> t5  = cc[0] * dd[2];
        complex<double> t6  = bb[2] * t2;
        complex<double> t8  = cc[0] * dd[2];
        complex<double> t9  = bb[1] * t2;
        complex<double> t11 = cc[2] * dd[1];
        complex<double> t13 = cc[2] * dd[0];
        complex<double> t15 = cc[1] * dd[0];
        complex<double> t17 = t1 * t3 + t5 * t6 - t8 * t9 - t11 * t3 + t13 * t9 - t15 * t6;
        //
        complex<double> c1 = t17 / 6.0;
        complex<double> c2 = -t8 * bb[1] / 6.0 + t5 * bb[2] / 6.0 - t15 * bb[2] / 6.0 + t1 * bb[0] / 6.0 - t11 * bb[0] / 6.0 + t13 * bb[1] / 6.0;
        complex<double> c3 = -t17 / 6.0;
        //
        t3 = (1.0 / ko / Iunit);
        t5 = 2.0 * t3 * c3;
        t8 = pow(ko,2);
        //
        complex<double> t10 = pow(Iunit,2);
        complex<double> t14 = 3.0 / t8 / t10 * c2;
        complex<double> t16 = 2.0 * t3 * c1;
        //
        coef[0]  = -t5;
        coef[1]  = c1;
        coef[2]  = c2;
        coef[3]  = 3.0 * t3 * c2;
        coef[4]  = -t14;
        coef[5]  = -t16;
        coef[6]  = t5;
        coef[7]  = c3;
        coef[8]  = t16;
        coef[9]  = t14;
        //
        t1   = cc[1] * bb[2];
        t3   = cc[1] * bb[0];
        t5   = cc[2] * bb[1];
        complex<double> t7  = cc[0] * bb[2];
        t9   = cc[2] * bb[0];
        t11  = cc[0] * bb[1];
        t14  = sqrt(3.0);
        t15  = dd[0] * t14;
        t17  = dd[2] * t14;
        complex<double> t20  = dd[1] * t14;
        complex<double> t24  = t1 * t15 - t3 * t17 - t5 * t15 - t7 * t20 + t9 * t20 + t11 * t17;
        //
        complex<double> cm1 = -t1 * dd[0] / 6.0 + t3 * dd[2] / 6.0 + t5 * dd[0] / 6.0 + t7 * dd[1] / 6.0 - t9 * dd[1] / 6.0 - t11 * dd[2] / 6.0;
        complex<double> cm2 = t24 / 6.0;
        complex<double> cm3 = t24 / 6.0;
        //
        t3  = (1.0 / ko / Iunit);
        t5  = 2.0 * t3 * cm2;
        t6  = pow(ko,2);
        t8  = pow(Iunit,2);
        complex<double> t12 = 3.0 / t6 / t8 * cm1;
        t14 = 2.0 * t3 * cm3;
        //
        coefm[0]  = -t5;
        coefm[1]  = -t12;
        coefm[2]  = t14;
        coefm[3]  = -t14;
        coefm[4]  = cm1;
        coefm[5]  = cm2;
        coefm[6]  = 3.0 * t3 * cm1;
        coefm[7]  = t5;
        coefm[8]  = t12;
        coefm[9]  = cm3;
        
    } // if (flag == 31)
    
    if (flag == 32)
    {
        complex<double> c[3], cm[3];
        //
        complex<double> t1 = cc[2] * bb[0];
        complex<double> t3 = cc[1] * bb[0];
        complex<double> t5 = cc[2] * bb[1];
        complex<double> t7 = cc[0] * bb[2];
        complex<double> t9 = cc[0] * bb[1];
        complex<double> t11 = cc[1] * bb[2];
        complex<double> t14 = sqrt(0.3e1);
        complex<double> t15 = dd[2] * t14;
        complex<double> t18 = dd[0] * t14;
        complex<double> t20 = dd[1] * t14;
        complex<double> t24 = -t9 * t15 + t3 * t15 + t5 * t18 + t7 * t20 - t1 * t20 - t11 * t18;
        //
        c[0] = t1 * dd[1] / 0.6e1 - t3 * dd[2] / 0.6e1 - t5 * dd[0] / 0.6e1 - t7 * dd[1] / 0.6e1 + t9 * dd[2] / 0.6e1 + t11 * dd[0] / 0.6e1;
        c[1] = t24 / 0.6e1;
        c[2] = t24 / 0.6e1;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * c[1];
        t7 = 2.0 * t3 * c[2];
        complex<double> t8 = j * j;
        complex<double> t10 = ko * ko;
        t14 = 3.0 / t8 / t10 * c[0];
        //
        coef[0] = -t5;
        coef[1] = t7;
        coef[2] = -t14;
        coef[3] = -t7;
        coef[4] = c[2];
        coef[5] = t5;
        coef[6] = t14;
        coef[7] = c[0];
        coef[8] = c[1];
        coef[9] = 3.0 * t3 * c[0];
        //
        t1 = cc[2] * bb[0];
        t3 = cc[1] * bb[0];
        t5 = cc[2] * bb[1];
        t7 = cc[0] * bb[2];
        t9 = cc[0] * bb[1];
        t11 = cc[1] * bb[2];
        t14 = sqrt(0.3e1);
        t15 = dd[2] * t14;
        complex<double> t17 = dd[1] * t14;
        t20 = dd[0] * t14;
        t24 = -t3 * t15 + t1 * t17 + t9 * t15 + t11 * t20 - t5 * t20 - t7 * t17;
        //
        cm[0] = t1 * dd[1] / 0.6e1 - t3 * dd[2] / 0.6e1 - t5 * dd[0] / 0.6e1 - t7 * dd[1] / 0.6e1 + t9 * dd[2] / 0.6e1 + t11 * dd[0] / 0.6e1;
        cm[1] = t24 / 0.6e1;
        cm[2] = -t24 / 0.6e1;
        //
        t1 = j * j;
        t3 = ko * ko;
        t7 = 3.0 / t1 / t3 * cm[0];
        t10 = 0.1e1 / j / ko;
        complex<double> t12 = 2.0 * t10 * cm[2];
        complex<double> t16 = 2.0 * t10 * cm[1];
        //
        coefm[0] = -t7;
        coefm[1] = t12;
        coefm[2] = -t12;
        coefm[3] = cm[1];
        coefm[4] = cm[0];
        coefm[5] = 3.0 * t10 * cm[0];
        coefm[6] = t16;
        coefm[7] = t7;
        coefm[8] = cm[2];
        coefm[9] = -t16;
        
    } // if (flag == 32)

    if (flag == 33)
    {
        complex<double> c[4], cm[4];
        //
        complex<double> t1 = cc[0] * bb[1];
        complex<double> t2 = sqrt(0.3e1);
        complex<double> t3 = dd[2] * t2;
        complex<double> t5 = cc[0] * bb[2];
        complex<double> t6 = dd[1] * t2;
        complex<double> t8 = cc[1] * bb[2];
        complex<double> t9 = dd[0] * t2;
        complex<double> t11 = cc[1] * bb[0];
        complex<double> t13 = cc[2] * bb[0];
        complex<double> t15 = cc[2] * bb[1];
        complex<double> t17 = t1 * t3 - t5 * t6 + t8 * t9 - t11 * t3 + t13 * t6 - t15 * t9;
        //
        c[0] = t17 / 0.6e1;
        c[1] = -t17 / 0.6e1;
        c[2] = -t17 / 0.6e1;
        c[3] = t1 * dd[2] / 0.2e1 - t5 * dd[1] / 0.2e1 + t8 * dd[0] / 0.2e1 - t11 * dd[2] / 0.2e1 + t13 * dd[1] / 0.2e1 - t15 * dd[0] / 0.2e1;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * c[3];
        complex<double> t7 = 2.0 * t3 * c[1];
        t9 = 2.0 * t3 * c[0];
        complex<double> t10 = j * j;
        complex<double> t12 = ko * ko;
        complex<double> t16 = 3.0 / t10 / t12 * c[2];
        //
        coef[0] = -t5;
        coef[1] = -t7;
        coef[2] = -t9;
        coef[3] = -t16;
        coef[4] = t5;
        coef[5] = c[1];
        coef[6] = c[0];
        coef[7] = c[2];
        coef[8] = 3.0 * t3 * c[2];
        coef[9] = c[3];
        coef[10] = t7;
        coef[11] = t9;
        coef[12] = t16;
        //
        t1 = cc[1] * bb[2];
        t2 = sqrt(0.3e1);
        t3 = dd[0] * t2;
        t5 = cc[1] * bb[0];
        t6 = dd[2] * t2;
        t8 = cc[2] * bb[0];
        t9 = dd[1] * t2;
        t11 = cc[2] * bb[1];
        t13 = cc[0] * bb[2];
        t15 = cc[0] * bb[1];
        t17 = t1 * t3 - t5 * t6 + t8 * t9 - t11 * t3 - t13 * t9 + t15 * t6;
        //
        cm[0] = t17 / 0.6e1;
        cm[1] = t17 / 0.6e1;
        cm[2] = t11 * dd[0] / 0.2e1 - t8 * dd[1] / 0.2e1 - t15 * dd[2] / 0.2e1 + t5 * dd[2] / 0.2e1 - t1 * dd[0] / 0.2e1 + t13 * dd[1] / 0.2e1;
        cm[3] = -t17 / 0.6e1;
        //
        t3 = 0.1e1 / j / ko;
        t5 = 2.0 * t3 * cm[1];
        t7 = 2.0 * t3 * cm[3];
        t8 = j * j;
        t10 = ko * ko;
        complex<double> t14 = 3.0 / t8 / t10 * cm[0];
        t16 = 2.0 * t3 * cm[2];
        //
        coefm[0] = cm[2];
        coefm[1] = t5;
        coefm[2] = t7;
        coefm[3] = t14;
        coefm[4] = -t16;
        coefm[5] = -t5;
        coefm[6] = -t7;
        coefm[7] = -t14;
        coefm[8] = t16;
        coefm[9] = cm[1];
        coefm[10] = cm[3];
        coefm[11] = cm[0];
        coefm[12] = 3.0 * t3 * cm[0];
        
    } // if (flag == 33)
}
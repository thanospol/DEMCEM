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

#include "demcem_ws_st_rwg.h"
#include "demcem_constants.h"

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> PHI_1_1
// ***********************************************************************

complex<double> phi_ws_st_rwg ( complex<double> A, complex<double> B, complex<double> C, int argument, double psi, double Apsi, const double ko , int flag)
{
    complex<double> PHI_out;
    complex<double> j   = Iunit;
    
    if (flag == 11)
    {
        
        complex<double> PSI_1 ;
        complex<double> PSI_2 ;
        complex<double> PSI_3 ;
        complex<double> PSI_ ;
        
        double beta ;
        double gamma ;
        double epsilon ;
        double delta ;
        
        // Int_1_1
        complex<double> D = sin(psi) / ( pow(j*ko,2.0) * pow(Apsi,3.0) );
        //
        switch(argument)
        {
            case 1:
                // PHI_a
                PSI_1 = A / (2.0 * B);
                PSI_2 = -1.0;
                PSI_3 = (B / A) * (1.0 - exp(-A/B));
                break;
            case 2:
                // PHI_b
                PSI_1 =  A / (2.0 * C);
                PSI_2 = -1.0;
                PSI_3 = (C / A) * (1.0 - exp(-A/C));
                break;
            case 3:
                // PHI_c
                beta    = tan(M_PI - psi) / sqrt(3.0);
                gamma   = (1.0-beta) / (1.0+beta);
                
                PSI_1 =  (A/C) * (1.0/2.0 -(gamma-pow(gamma,2.0) / 2.0) );
                PSI_2 = - (1.0 - gamma);
                PSI_3 = (C / A) * (1.0 - exp(-A * (1-gamma) / C) );
                break;
            case 4:
                // PHI_d
                beta    = tan(M_PI - psi) / sqrt(3.0);
                gamma   = (1.0-beta) / (1.0+beta);
                
                PSI_1 = -A * (gamma + pow(gamma,2.0) / 2.0) / B;
                PSI_2 = -gamma;
                PSI_3 = (B / A) * (exp(A * (1+gamma) / B) - exp(A/B));
                break;
            case 5:
                // PHI_e
                epsilon = tan(psi) / sqrt(3.0);
                delta   = -(1-epsilon) / (1+epsilon);
                
                PSI_1 = A * (pow(delta,2.0) / 2.0 - delta) / B;
                PSI_2 = delta;
                PSI_3 = (B / A) * (exp(-A/B) - exp(-A * (1-delta) / B) );
                break;
            case 6:
                // PHI_f
                epsilon = tan(psi) / sqrt(3.0);
                delta   = -(1-epsilon) / (1+epsilon);
                
                PSI_1 = (A / C) * (delta + pow(delta,2.0) / 2.0 + 1.0/2.0);
                PSI_2 = -(1.0 + delta);
                PSI_3 = (C / A) * (1.0 - exp(-A * (1 + delta) / C));
                break;
            case 7:
                // PHI_g
                PSI_1 = A / (2.0 * C);
                PSI_2 = -1.0;
                PSI_3 = (C / A) * (1.0 - exp(-A / C) );
                break;
            case 8:
                // PHI_h
                PSI_1 = -A / (2.0 * B);
                PSI_2 = -1.0;
                PSI_3 = (B / A) * (exp(A/B) - 1.0);
                break;
        }// end switch argument
        
        PSI_  = PSI_1 + PSI_2 + PSI_3;
        //
        PHI_out   = D * PSI_;
    } // if (flag = 11)
    
    if (flag == 12)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2_1 ;
            complex<double> PSI_2_2 ;
            complex<double> PSI_2_3 ;
            complex<double> PSI_2 ;
            complex<double> PSI_3 ;
            complex<double> PSI_4 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_2
            complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_0 = 1.0 / (6.0 * B);
                    PSI_1 = 1.0 / (2.0 * B);
                    
                    PSI_2_1 = 1.0;
                    PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
                    PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
                    
                    PSI_3  = 1.0 - B / A * (1.0 - exp(-A / B) );
                    
                    PSI_4  = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
                    
                    break;
                case 2:
                    // PHI_b
                    PSI_0 = 1.0 / (6.0 * C);
                    PSI_1 = 1.0 / (2.0 * C);
                    
                    PSI_2_1 = 1.0;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
                    
                    PSI_3  = 1.0 - C / A * (1.0 - exp(-A / C) );
                    
                    PSI_4  = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
                    
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
                    PSI_1 = (1.0 / C) * (1.0 / 2.0 -(gamma - pow(gamma,2.0) / 2.0) );
                    
                    PSI_2_1 = 1.0 - gamma;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 -gamma) / C) );
                    PSI_2_3 = -1.0 /A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
                    
                    PSI_3  = (1.0 - gamma) - C / A * (1.0 - exp(-A * (1.0-gamma) / C) );
                    
                    PSI_4  = 1.0 /2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
                    PSI_1 = -1.0 / B * (gamma + pow(gamma,2.0) / 2.0 );
                    
                    PSI_2_1 = gamma;
                    PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B));
                    PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
                    
                    PSI_3  = gamma - B / A * (exp(A * (1.0 + gamma) / B) - exp(A /B) );
                    
                    PSI_4  = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
                    PSI_1 = 1.0 / B * (pow(delta,2.0) / 2.0 - delta );
                    
                    PSI_2_1 = -delta;
                    PSI_2_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
                    
                    PSI_3  = -delta - B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    
                    PSI_4  = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 /3.0 * (pow(delta,3.0) + 1.0) );
                    PSI_1 = 1.0 / C * (delta + pow(delta,2.0) / 2.0 + 1.0 / 2.0);
                    
                    PSI_2_1 = 1.0 + delta;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    PSI_2_3 = 1.0 / A * (( A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C );
                    
                    PSI_3  = (1.0 + delta) - C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    
                    PSI_4  = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    
                    break;
                case 7:
                    // PHI_g
                    PSI_0 = -1.0 / (6.0 * C);
                    PSI_1 = 1.0 / (2.0 * C);
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
                    
                    PSI_3  = 1.0 - C / A * (1.0 - exp(-A / C) );
                    
                    PSI_4  = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    
                    break;
                case 8:
                    // PHI_h
                    PSI_0 = 1.0 / (6.0 * B);
                    PSI_1 = -1.0 / (2.0 * B);
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -B / A * (exp(A / B) - 1.0);
                    PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
                    
                    PSI_3  = 1.0 - B / A * (exp(A / B) - 1.0);
                    
                    PSI_4  = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B) );
                    
                    break;
            }// end switch argument
            
            PSI_2   = PSI_2_1 + PSI_2_2 + PSI_2_3;
            //
            PSI_ = j * ko *PSI_0 + (cos(psi) / Apsi) * PSI_1 - ((cos(psi) / (j * ko * pow(Apsi,2.0) ) ) * (PSI_2 + PSI_3) + (1.0 / Apsi) * PSI_4);
            //
            PHI_out    = D * PSI_;
        
        } // if (flag = 12)
    
    if (flag == 13)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_3
            complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_0 = 1.0 / (3.0 * pow(B,2.0));
                    PSI_1 = 1.0 / (2.0 * B);
                    PSI_2 = 1.0 - B / A * (1.0 - exp(-A / B) );
                    break;
                case 2:
                    // PHI_b
                    PSI_0 =  1.0 / (3.0 * pow(C,2.0));
                    PSI_1 = 1.0 / (2.0 * C);
                    PSI_2 = 1.0 - C / A * (1.0 - exp(-A / C) );
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 =  pow(1.0-gamma,3.0) / (3.0 * pow(C,2.0) );
                    PSI_1 = (1.0 / C) * (1.0 / 2.0 - (gamma - pow(gamma,2.0) / 2.0 ) );
                    PSI_2 = (1.0 - gamma) - C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = (pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0) );
                    PSI_1 = -1.0 / B * (gamma + pow(gamma,2.0) / 2.0);
                    PSI_2 = gamma - B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B));
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 =  (pow(1.0-delta,3.0) - 1.0) / (3.0 * pow(B,2.0) );
                    PSI_1 = (1.0 / B) * ( pow(delta,2.0) / 2.0  - delta );
                    PSI_2 = -delta - B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B));
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));
                    PSI_1 = 1.0 / C * (delta + pow(delta,2.0) / 2.0 + 1.0 / 2.0);
                    PSI_2 = (1.0 + delta) - C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    break;
                case 7:
                    // PHI_g
                    PSI_0 = 1.0 / (3.0 * pow(C,2.0) );
                    PSI_1 = 1.0 / (2.0 * C);
                    PSI_2 = 1.0 - C / A * (1.0 - exp(-A / C) );
                    break;
                case 8:
                    // PHI_h
                    PSI_0 = 1.0 / (3.0 * pow(B,2.0) );
                    PSI_1 = -1.0 / (2.0 * B);
                    PSI_2 = 1.0 - B / A * (exp(A / B) - 1.0);
                    break;
            }// end switch argument
            
            PSI_  = (j * ko * sin(psi) / 2.0) * PSI_0 - (sin(psi) / Apsi) * PSI_1 + (sin(psi) / (j * ko * pow(Apsi,2.0) ) ) * PSI_2;
            //
            PHI_out    = D * PSI_;
            
         } // if (flag = 13)
    
    if (flag == 21)
        {
            complex<double> PSI_1 ;
            complex<double> PSI_2 ;
            complex<double> PSI_3 ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_3
            complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,3.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_1 = A / (6.0 * B);
                    PSI_2 = -1.0 / 2.0;
                    PSI_3 = B / A - (pow(B,2.0) / pow(A,2.0) ) * (1.0 - exp(-A / B) );
                    break;
                case 2:
                    // PHI_b
                    PSI_1 = A / (6.0 * C);
                    PSI_2 = -1.0 / 2.0;
                    PSI_3 = C / A - (pow(C,2.0) / pow(A,2.0) ) * (1.0 - exp(-A / C) );
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_1 = (A / C) * (1.0 / 2.0 * (1.0 - pow(gamma,2.0) ) - 1.0 / 3.0 * (1.0 - pow(gamma,3.0) ) );
                    PSI_2 = -1.0 / 2.0 * (1.0 - pow(gamma,2.0) );
                    PSI_3 = -(C / pow(A,2.0) ) * (C - A + (A * gamma - C) * exp(A * (gamma - 1.0) / C ) );
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_1 = -(A / B) * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0 );
                    PSI_2 = -pow(gamma,2.0) / 2.0;
                    PSI_3 = (B / pow(A,2.0)) * ((A * gamma - B) * exp(A * (1.0 + gamma) / B) + B * exp(A / B) );
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_1 = A / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
                    PSI_2 = pow(delta,2.0) / 2.0;
                    PSI_3 = -(B / pow(A,2.0)) * (B * exp(-A / B) + (A * delta - B) * exp(A * (delta - 1.0) / B ) );
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_1 = (A / C ) * (1.0 / 2.0 *(pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0 ) );
                    PSI_2 = -1.0 / 2.0 * (pow(delta,2.0) - 1.0);
                    PSI_3 = -(C / pow(A,2.0) ) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C ) );
                    break;
                case 7:
                    // PHI_g
                    PSI_1 = -A / (6.0 * C);
                    PSI_2 = 1.0 / 2.0;
                    PSI_3 = -(C / pow(A,2.0)) * (C * exp(-A / C) - C + A);
                    break;
                case 8:
                    // PHI_h
                    PSI_1 = A / (6.0 * B);
                    PSI_2 = 1.0 / 2.0;
                    PSI_3 = (B / pow(A,2.0)) * (A + B - B * exp(A / B) );
                    break;
            }// end switch argument
            
            //
            PHI_out    = D * (PSI_1 + PSI_2 + PSI_3);
        } // if (flag = 21)
    
    if (flag == 22)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2 ;
            complex<double> PSI_3 ;
            complex<double> PSI_3_1 ;
            complex<double> PSI_3_2 ;
            complex<double> PSI_4 ;
            complex<double> PSI_5 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_2_2
            complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_0 = 1.0 / (6.0 * B);
                    PSI_1 = 1.0 / (12.0 * B);
                    PSI_2 = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
                    
                    PSI_5  = B / A - 2.0 * (pow(B,2.0) / pow(A,2.0)) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (1.0 - exp(-A / B) );
                    
                    PSI_3_1 = B / A - (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B));
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / B * (PSI_3_1 - PSI_3_2);
                    
                    PSI_4  = 1.0 / 3.0 ;
                    
                    break;
                case 2:
                    // PHI_b
                    PSI_0 = 1.0 / (6.0 * C);
                    PSI_1 = 1.0 / (12.0 * C);
                    PSI_2 = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
                    
                    PSI_5  = C / A - 2.0 * (pow(C,2.0) / pow(A,2.0)) + 2.0 * (pow(C,3.0) / pow(A,3.0)) * (1.0 - exp(-A / C) );
                    
                    PSI_3_1 = C / A - (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C));
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / C * (PSI_3_1 - PSI_3_2);
                    
                    PSI_4  = 1.0 / 3.0 ;
                    
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
                    PSI_1 = (1.0 / C) * (1.0 / 12.0 - (pow(gamma,3.0) / 3.0 - pow(gamma,4.0) / 4.0) );
                    PSI_2 = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1 - gamma) / C) );
                    
                    PSI_5  = C / A * (1.0 - pow(gamma,2.0) * exp(-A * (1.0 - gamma) / C) ) - 2.0 * (pow(C,2.0) / pow(A,2.0)) * (1.0- gamma * exp(-A * (1.0 - gamma) / C) ) + 2.0 * (pow(C,3.0) / pow(A,3.0)) * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    
                    PSI_3_1 = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / C * (PSI_3_1 - PSI_3_2);
                    
                    PSI_4  = (1.0 - pow(gamma,3.0) ) / 3.0 ;
                    
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
                    PSI_1 = -1.0 / B * (pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0);
                    PSI_2 = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    
                    PSI_5  = B / A * pow(gamma,2.0) * exp(A * (1.0 + gamma) / B) - 2.0 * (pow(B,2.0) / pow(A,2.0)) * gamma * exp(A * (1.0 + gamma) / B) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    
                    PSI_3_1 = B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    PSI_3_2 = PSI_5;
                    PSI_3   = -1.0 / B * (PSI_3_1 + PSI_3_2);
                    
                    PSI_4  = pow(gamma,3.0) / 3.0;
                    
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
                    PSI_1 = 1.0 / B * (pow(delta,4.0) / 4.0 - pow(delta,3.0) / 3.0);
                    PSI_2 = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    
                    PSI_5  = -B / A * pow(delta,2.0) * exp(-A * (1.0 - delta) / B) + 2.0 * (pow(B,2.0) / pow(A,2.0)) * delta * exp(-A * (1.0 - delta) / B) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    
                    PSI_3_1 = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / B * (PSI_3_1 - PSI_3_2);
                    
                    PSI_4  = -pow(delta,3.0) / 3.0;
                    
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0) );
                    PSI_1 = 1.0 / C * (pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0 + 1.0 / 12.0);
                    PSI_2 = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    
                    PSI_5  = C / A * (1.0 - pow(delta,2.0) * exp(-A * (1.0 + delta) / C)) - 2.0 * (pow(C,2.0) / pow(A,2.0)) * (delta * exp(-A * (1.0 + delta) / C) + 1.0) - 2.0 * (pow(C,3.0) / pow(A,3.0)) * (exp(-A * (1.0 + delta) / C) - 1.0);
                    
                    PSI_3_1 = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / C * (PSI_3_1 + PSI_3_2);
                    
                    PSI_4  = (pow(delta,3.0) + 1.0) / 3.0;
                    
                    break;
                case 7:
                    // PHI_g
                    PSI_0 = -1.0 / (6.0 * C);
                    PSI_1 = 1.0 / (12.0 * C);
                    PSI_2 = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    
                    PSI_5  = C / A - 2.0 * (pow(C,2.0) / pow(A,2.0)) - 2.0 * (pow(C,3.0) / pow(A,3.0)) * (exp(-A / C) - 1.0);
                    
                    PSI_3_1 = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    PSI_3_2 = PSI_5;
                    PSI_3   = 1.0 / C * (PSI_3_1 + PSI_3_2);
                    
                    PSI_4  = 1.0 / 3.0 ;
                    
                    break;
                case 8:
                    // PHI_h
                    PSI_0 = 1.0 / (6.0 * B);
                    PSI_1 = -1.0 / (12.0 * B);
                    PSI_2 = -1.0 / 2.0 - B / pow(A,2.0) * ( A + B - B * exp(A / B) );
                    
                    PSI_5  = -B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * (pow(B,3.0) / pow(A,3.0)) * (exp(A / B) - 1.0);
                    
                    PSI_3_1 = B / pow(A,2.0) * (A + B - B * exp(A / B) );
                    PSI_3_2 = PSI_5;
                    PSI_3   = -1.0 / B * (PSI_3_1 + PSI_3_2);
                    
                    PSI_4  = 1.0 / 3.0 ;
                    
                    break;
            }// end switch argument
            
            PSI_ = (cos(psi) / Apsi) * PSI_0 + j * ko * PSI_1 - (2.0 * cos(psi) / (j * ko * pow(Apsi,2.0))) * PSI_2 + (cos(psi) / Apsi) * PSI_3 - (1 / Apsi) * (PSI_4 - PSI_5);
            //
            PHI_out    = D * PSI_;
        } // if (flag = 22)
    
    if (flag == 23)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_3
            complex<double> D = sin(psi) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_0 = 1.0 / (12.0 * pow(B,2.0));
                    PSI_1 = 1.0 / (6.0 * B);
                    PSI_2 = 1.0/ 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
                    break;
                case 2:
                    // PHI_b
                    PSI_0 =  1.0 / (12.0 * pow(C,2.0));
                    PSI_1 = 1.0 / (6.0 * C);
                    PSI_2 = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 =  (1.0 / pow(C,2.0)) * (1.0 / 12.0 - ( pow(gamma,2.0) / 2.0 - 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );
                    PSI_1 = (1.0 / C) * ((1.0 / 2.0) * (1.0 - pow(gamma,2.0)) - (1.0 / 3.0) * (1.0 - pow(gamma,3.0) ) );
                    PSI_2 = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = (1.0 / pow(B,2.0)) * ( pow(gamma,2.0) / 2.0 + 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0);
                    PSI_1 = -1.0 / B * (pow(gamma,2.0) / 2.0 + pow(gamma,3.0) / 3.0);
                    PSI_2 = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 =  -(1.0 / pow(B,2.0)) * (pow(delta,2.0) / 2.0 -2.0* pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0);
                    PSI_1 = 1.0 / B * (pow(delta,3.0) / 3.0 - pow(delta,2.0) / 2.0);
                    PSI_2 = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = (1.0 / pow(C,2.0)) * ((pow(delta,2.0) / 2.0 + 2.0* pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0) - 1.0 / 12.0);
                    PSI_1 = 1.0 / C * (1.0 / 2.0 * (pow(delta,2.0) - 1.0) + 1.0 / 3.0 * (pow(delta,3.0) + 1.0) );
                    PSI_2 = 1.0 / 2.0 *(pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    break;
                case 7:
                    // PHI_g
                    PSI_0 = -1.0 / (12.0 * pow(C,2.0) );
                    PSI_1 = -1.0 / (6.0 * C);
                    PSI_2 = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    break;
                case 8:
                    // PHI_h
                    PSI_0 = -1.0 / (12.0 * pow(B,2.0) );
                    PSI_1 = 1.0 / (6.0 * B);
                    PSI_2 = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B) );
                    break;
            }// end switch argument
            
            PSI_  = (j * ko * sin(psi) / 2.0) * PSI_0 - (sin(psi) / Apsi) * PSI_1 + (sin(psi) / (j * ko * pow(Apsi,2.0)) ) * PSI_2;
            //
            PHI_out    = D * PSI_;
        } // if (flag = 23)
    
    if (flag == 31)
        {
            complex<double> PSI_2 ;
            complex<double> PSI_5_1 ;
            complex<double> PSI_5_2 ;
            complex<double> PSI_5_3 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_3
            complex<double> D = pow(sin(psi),2.0) / ( j * ko * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_2 = 1.0 / (3.0 * pow(B,2.0));
                    
                    PSI_5_1 = 1.0 ;
                    PSI_5_2 = -B / A * (1.0 - exp(-A / B) );
                    PSI_5_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
                    break;
                case 2:
                    // PHI_b
                    PSI_2 = 1.0 / (3.0 * pow(C,2.0));
                    
                    PSI_5_1 = 1.0 ;
                    PSI_5_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_5_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_2 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));
                    
                    PSI_5_1 = 1.0 - gamma ;
                    PSI_5_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    PSI_5_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_2 = (pow(1.0 + gamma,3.0) - 1.0 ) / (3.0 * pow(B,2.0) );
                    
                    PSI_5_1 = gamma ;
                    PSI_5_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    PSI_5_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_2 = (pow(1.0 - delta,3.0) - 1.0 ) / (3.0 * pow(B,2.0) );
                    
                    PSI_5_1 = -delta ;
                    PSI_5_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    PSI_5_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_2 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0) );
                    
                    PSI_5_1 = 1.0 + delta ;
                    PSI_5_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    PSI_5_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
                    break;
                case 7:
                    // PHI_g
                    PSI_2 = 1.0 / (3.0 * pow(C,2.0) );
                    
                    PSI_5_1 = 1.0 ;
                    PSI_5_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_5_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
                    break;
                case 8:
                    // PHI_h
                    PSI_2 = 1.0 / (3.0 * pow(B,2.0) );
                    
                    PSI_5_1 = 1.0 ;
                    PSI_5_2 = -B / A * (exp(A / B) - 1.0);
                    PSI_5_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
                    break;
            }// end switch argument
            
            PSI_  = PSI_2 / 2.0 - (PSI_5_1 + PSI_5_2 + PSI_5_3) / pow(j * ko * Apsi,2.0);
            //
            PHI_out    = D * PSI_;
        } // if (flag = 31)
    
    if (flag == 32)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2_1 ;
            complex<double> PSI_2_2 ;
            complex<double> PSI_2_3 ;
            complex<double> PSI_2 ;
            complex<double> PSI_3 ;
            complex<double> PSI_3_1 ;
            complex<double> PSI_3_2 ;
            complex<double> PSI_3_3 ;
            complex<double> PSI_4 ;
            complex<double> PSI_5_1 ;
            complex<double> PSI_5_2 ;
            complex<double> PSI_5 ;
            complex<double> PSI_5_2_2 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_2
            complex<double> D = pow(sin(psi),2.0) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_5_2_2 = B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * pow(B,3.0) / pow(A,3.0) * (1.0 - exp(-A / B) );
                    
                    PSI_0 = 1.0 / (3.0 * pow(B,2.0));
                    PSI_1 = 1.0 / (12.0 * pow(B,2.0));
                    
                    PSI_2_1 = 1.0;
                    PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
                    PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
                    
                    PSI_3_1  = B / A * (1.0 - exp(-A / B) );
                    PSI_3_2  = B / A - pow(B,2.0) / pow(A,2.0) * (1.0 - exp(-A / B) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = 1.0 / 2.0 - B / A + (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
                    
                    PSI_5_1  = B / A - (pow(B,2.0) / pow(A,2.0)) * (1.0 - exp(-A / B) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / B * (PSI_5_1 - PSI_5_2);
                    
                    break;
                case 2:
                    // PHI_b
                    PSI_5_2_2 = C / A - 2.0 * pow(C,2.0) / pow(A,2.0) + 2.0 * pow(C,3.0) / pow(A,3.0) * (1.0 - exp(-A / C) );
                    
                    PSI_0 = 1.0 / (3.0 * pow(C,2.0));
                    PSI_1 = 1.0 / (12.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
                    
                    PSI_3_1  = C / A * (1.0 - exp(-A / C) );
                    PSI_3_2  = C / A - pow(C,2.0) / pow(A,2.0) * (1.0 - exp(-A / C) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = 1.0 / 2.0 - C / A + (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
                    
                    PSI_5_1  = C / A - (pow(C,2.0) / pow(A,2.0)) * (1.0 - exp(-A / C) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / C * (PSI_5_1 - PSI_5_2);
                    
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_5_2_2 = C / A * (1.0 - pow(gamma,2.0) * exp(-A * (1.0 - gamma) / C) ) - 2.0 * pow(C,2.0) / pow(A,2.0) * (1.0 - gamma * exp(-A * (1.0 - gamma) / C) ) + 2.0 * pow(C,3.0) / pow(A,3.0) * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    
                    PSI_0 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));
                    PSI_1 = (1.0 / pow(C,2.0) ) * (1.0 / 12.0 - (pow(gamma,2.0) / 2.0 - 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );
                    
                    PSI_2_1 = 1.0 - gamma;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    PSI_2_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
                    
                    PSI_3_1  = C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    PSI_3_2  = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = 1.0 / 2.0 * (1.0 - pow(gamma,2.0)) + C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    
                    PSI_5_1  = -C / pow(A,2.0) * (C - A + (A * gamma - C) * exp(-A * (1.0 - gamma) / C) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / C * (PSI_5_1 - PSI_5_2);
                    
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_5_2_2 = B / A * pow(gamma,2.0) * exp(A * (1.0 + gamma) / B) - 2.0 * pow(B,2.0) / pow(A,2.0) * gamma * exp(A * (1.0 + gamma) / B) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    
                    PSI_0 = (pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0));
                    PSI_1 = (1.0 / pow(B,2.0) ) * ((pow(gamma,2.0) / 2.0 + 2.0 * pow(gamma,3.0) / 3.0 + pow(gamma,4.0) / 4.0) );
                    
                    PSI_2_1 = gamma;
                    PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
                    
                    PSI_3_1  = B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    PSI_3_2  = B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = pow(gamma,2.0) / 2.0 - B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    
                    PSI_5_1  =  B / pow(A,2.0) * (B * exp(A / B) + (A * gamma - B) * exp(A * (1.0 + gamma) / B) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = -1.0 / B * (PSI_5_1 + PSI_5_2);
                    
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_5_2_2 = -B / A * pow(delta,2.0) * exp(-A * (1.0 - delta) / B) + 2.0 * pow(B,2.0) / pow(A,2.0) * delta * exp(-A * (1.0 - delta) / B) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    
                    PSI_0 = (pow(1.0 - delta,3.0) - 1.0) / (3.0 * pow(B,2.0));
                    PSI_1 = -(1.0 / pow(B,2.0) ) * ((pow(delta,2.0) / 2.0 - 2.0 * pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0) );
                    
                    PSI_2_1 = -delta;
                    PSI_2_2 = -B / A * (exp(- A / B) - exp(-A * (1.0 - delta) / B) );
                    PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
                    
                    PSI_3_1  = B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    PSI_3_2  = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 - 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = -pow(delta,2.0) / 2.0 + B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    
                    PSI_5_1  = -B / pow(A,2.0) * (B * exp(-A / B) + (A * delta - B) * exp(-A * (1.0 - delta) / B) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / B * (PSI_5_1 - PSI_5_2);
                    
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_5_2_2 = C / A * (1.0 - pow(delta,2.0) * exp(-A * (1.0 + delta) / C)) - 2.0 * pow(C,2.0) / pow(A,2.0) * (delta * exp(-A * (1.0 + delta) / C) + 1.0) - 2.0 * pow(C,3.0) / pow(A,3.0) * (exp(-A * (1.0 + delta) / C) - 1.0);
                    
                    PSI_0 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));
                    PSI_1 = 1.0 / pow(C,2.0) * (pow(delta,2.0) / 2.0 + 2.0 * pow(delta,3.0) / 3.0 + pow(delta,4.0) / 4.0 - 1.0 / 12.0);
                    
                    PSI_2_1 = 1.0 + delta;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    PSI_2_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
                    
                    PSI_3_1  = C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    PSI_3_2  = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = 1.0 / 2.0 * (pow(delta,2.0) - 1.0) + C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    
                    PSI_5_1  = -C / pow(A,2.0) * (A - C + (A * delta + C) * exp(-A * (1.0 + delta) / C) );
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / C * (PSI_5_1 + PSI_5_2);
                    
                    break;
                case 7:
                    // PHI_g
                    PSI_5_2_2 = C / A - 2.0 * pow(C,2.0) / pow(A,2.0) - 2.0 * pow(C,3.0) / pow(A,3.0) * (exp(-A / C) - 1.0);
                    
                    PSI_0 = 1.0 / (3.0 * pow(C,2.0));
                    PSI_1 = -1.0 / (12.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
                    
                    PSI_3_1  = C / A * (1.0 - exp(-A / C));
                    PSI_3_2  = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(C,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = -1.0 / 2.0 + C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    
                    PSI_5_1  = -C / pow(A,2.0) * (C * exp(-A / C) - C + A);
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = 1.0 / C * (PSI_5_1 + PSI_5_2);
                    
                    break;
                case 8:
                    // PHI_h
                    PSI_5_2_2 = -B / A - 2.0 * pow(B,2.0) / pow(A,2.0) + 2.0 * pow(B,3.0) / pow(A,3.0) * (exp(A / B) - 1.0);
                    
                    PSI_0 = 1.0 / (3.0 * pow(B,2.0));
                    PSI_1 = -1.0 / (12.0 * pow(B,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -B / A * ( exp(A / B) - 1.0 );
                    PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
                    
                    PSI_3_1  = B / A * (exp(A / B) -1.0);
                    PSI_3_2  = B / pow(A,2.0) * (A +B -B * exp(A / B));
                    PSI_3_3  = PSI_5_2_2;
                    
                    PSI_3  = 1.0 / pow(B,2.0) * (PSI_3_1 + 2.0 * PSI_3_2 + PSI_3_3);
                    
                    PSI_4  = -1.0 / 2.0 - B / pow(A,2.0) * (A + B - B * exp(A / B));
                    
                    PSI_5_1  = B / pow(A,2.0) * (A + B - B * exp(A / B));
                    PSI_5_2  = PSI_5_2_2;
                    PSI_5  = -1.0 / B * (PSI_5_1 + PSI_5_2);
                    
                    break;
            }// end switch argument
            
            PSI_2 = PSI_2_1 + PSI_2_2 + PSI_2_3;
            //
            PSI_ = (cos(psi) / (2.0 * Apsi)) * PSI_0 + (j * ko / 2.0) * PSI_1 - (3.0 * cos(psi) / (pow(j * ko,2.0) * pow(Apsi,3.0)) ) * PSI_2 + (cos(psi) / Apsi) * PSI_3 - (1.0 / (j * ko * pow(Apsi,2.0)) ) * PSI_4 + (1.0 / Apsi) * PSI_5;
            //
            PHI_out    = D * PSI_;
        } // if (flag = 32)
    
    if (flag == 33)
        {
            complex<double> PSI_0 ;
            complex<double> PSI_1 ;
            complex<double> PSI_2_1 ;
            complex<double> PSI_2_2 ;
            complex<double> PSI_2_3 ;
            complex<double> PSI_2 ;
            complex<double> PSI_ ;
            
            double beta ;
            double gamma ;
            double epsilon ;
            double delta ;
            
            // Int_1_3
            complex<double> D = pow(sin(psi),2.0) / ( pow(j * ko,2.0) * pow(Apsi,2.0) );
            //
            switch(argument)
            {
                case 1:
                    // PHI_a
                    PSI_0 = 1.0 / (4.0 * pow(B,3.0));
                    PSI_1 = 1.0 / (3.0 * pow(B,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -B / A * (1.0 - exp(-A / B) );
                    PSI_2_3 = -1.0 / A * (B - (A + B) * exp(-A / B) );
                    break;
                case 2:
                    // PHI_b
                    PSI_0 = 1.0 / (4.0 * pow(C,3.0));
                    PSI_1 = 1.0 / (3.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = -1.0 / A * (C - (A + C) * exp(-A / C) );
                    break;
                case 3:
                    // PHI_c
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = pow(1.0 - gamma,4.0) / (4.0 * pow(C,3.0));
                    PSI_1 = pow(1.0 - gamma,3.0) / (3.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0 - gamma ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 - gamma) / C) );
                    PSI_2_3 = -1.0 / A * (C + (A * gamma - C - A) * exp(-A * (1.0 - gamma) / C) );
                    break;
                case 4:
                    // PHI_d
                    beta    = tan(M_PI - psi) / sqrt(3.0);
                    gamma   = (1.0-beta) / (1.0+beta);
                    
                    PSI_0 = (1.0 - pow(1.0 + gamma,4.0)) / (4.0 * pow(B,3.0));
                    PSI_1 = ( pow(1.0 + gamma,3.0) - 1.0) / (3.0 * pow(B,2.0));
                    
                    PSI_2_1 = gamma ;
                    PSI_2_2 = -B / A * (exp(A * (1.0 + gamma) / B) - exp(A / B) );
                    PSI_2_3 = -1.0 / A * ((A - B) * exp(A / B) + (B - A - A * gamma) * exp(A * (1.0 + gamma) / B) );
                    break;
                case 5:
                    // PHI_e
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = (pow(1.0 - delta,4.0) - 1.0) / (4.0 * pow(B,3.0));
                    PSI_1 = (pow(1.0 - delta,3.0) - 1.0) / (3.0 * pow(B,2.0));
                    
                    PSI_2_1 = -delta ;
                    PSI_2_2 = -B / A * (exp(-A / B) - exp(-A * (1.0 - delta) / B) );
                    PSI_2_3 = -1.0 / A * ((A + B) * exp(-A / B) + (A * delta - A - B) * exp(-A * (1.0 - delta) / B) );
                    break;
                case 6:
                    // PHI_f
                    epsilon = tan(psi) / sqrt(3.0);
                    delta   = -(1-epsilon) / (1+epsilon);
                    
                    PSI_0 = pow(1.0 + delta,4.0) / (4.0 * pow(C,3.0));
                    PSI_1 = pow(1.0 + delta,3.0) / (3.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0 + delta ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A * (1.0 + delta) / C) );
                    PSI_2_3 = 1.0 / A * ((A + A * delta + C) * exp(-A * (1.0 + delta) / C) - C);
                    break;
                case 7:
                    // PHI_g
                    PSI_0 =  1.0 / (4.0 * pow(C,3.0));
                    PSI_1 =  1.0 / (3.0 * pow(C,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -C / A * (1.0 - exp(-A / C) );
                    PSI_2_3 = 1.0 / A * ((A + C) * exp(-A / C) - C);
                    break;
                case 8:
                    // PHI_h
                    PSI_0 = - 1.0 / (4.0 * pow(B,3.0));
                    PSI_1 =  1.0 / (3.0 * pow(B,2.0));
                    
                    PSI_2_1 = 1.0 ;
                    PSI_2_2 = -B / A * (exp(A / B) - 1.0);
                    PSI_2_3 = 1.0 / A * ((A - B) * exp(A / B) + B);
                    break;
            }// end switch argument
            
            PSI_2  = PSI_2_1+PSI_2_2+PSI_2_3;
            //
            PSI_ = (j * ko * sin(psi) / 3.0) * PSI_0 - (sin(psi) / (2.0 * Apsi)) *PSI_1 + (sin(psi) / (pow(j * ko,2.0) * pow(Apsi,3.0))) *PSI_2;
            //
            PHI_out    = D * PSI_;
        } // if (flag = 33)

	// Final Output
	return PHI_out;
}
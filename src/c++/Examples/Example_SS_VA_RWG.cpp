//
//  Example_WS_VA_RWG.cpp
//
//  Created by Athanasios Polimeridis on 2/17/15.
//  Copyright (c) 2015 Athanasios Polimeridis. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <complex>
#include <ctime>

#include "demcem_ss_va_rwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"

void demcem_ss_va_rwg (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const complex<double> ko, const int N_theta_p, const int N_theta_q, const int N_psi, complex<double> I_DE[] );

void Example_SS_VA_RWG (const int Counter_start,  const int Counter_end )
{
    
    int N_theta_p, N_theta_q, N_psi;
    
    
    //
    double tstart, tstop, ttime;
    double RunTime;
    
    //
    complex<double> I_DE[9];
    
    const double r1[] = {0.1 , 0.1 , 0.1};
    const double r2[] = {0.2 , 0.1 , 0.1};
    const double r3[] = {0.1 , 0.2 , 0.1};
    const double r4[] = {0.0 , 0.1 , 0.2};
    const double r5[] = {0.0 , 0.2 , 0.2};
    
    const double ko = 2.0 * M_PI;
    
    //////////////////////////////////////////////////
    //                      VA
    //////////////////////////////////////////////////
    for (int points = Counter_start; points < Counter_end+1; points++)
    {
        
        N_theta_p = points;
        N_theta_q = points;
        N_psi = points;
        // Time keeper
        tstart = (double)clock()/CLOCKS_PER_SEC;
        
        //
        for (int i=0; i<100; i++)
        {
            demcem_ss_va_rwg ( r1, r2, r3, r4, r5, ko, N_theta_p, N_theta_q, N_psi, I_DE );
            
        } // for
        //
        tstop = (double)clock()/CLOCKS_PER_SEC;
        ttime= tstop-tstart; /*ttime is how long your code run */
        //
        RunTime = ttime * 1000 / 100;
        //
        cout << endl;
        cout << "SS - VA" << endl;
        
        cout << "Points: " << points << endl;
        cout << "Runtime: " << setprecision (6) << RunTime << " [msec]" << endl << endl;
        cout << "I_f1_f1 = " << setprecision (20) << I_DE[0] << endl ;
        cout << "I_f1_f2 = " << setprecision (20) << I_DE[1] << endl ;
        cout << "I_f1_f3 = " << setprecision (20) << I_DE[2] << endl ;
        cout << "I_f2_f1 = " << setprecision (20) << I_DE[3] << endl ;
        cout << "I_f2_f2 = " << setprecision (20) << I_DE[4] << endl ;
        cout << "I_f2_f3 = " << setprecision (20) << I_DE[5] << endl ;
        cout << "I_f3_f1 = " << setprecision (20) << I_DE[6] << endl ;
        cout << "I_f3_f2 = " << setprecision (20) << I_DE[7] << endl ;
        cout << "I_f3_f3 = " << setprecision (20) << I_DE[8] << endl ;
    } // for (int points = 10; points < Counter+1; points++)
    
    
}




//
//  Example_WS_ST_RWG.c
//
//  Created by Athanasios Polimeridis on 2/17/15.
//  Copyright (c) 2015 Athanasios Polimeridis. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <complex>
#include <ctime>

#include "demcem_ws_st_rwg.h"
#include "demcem_inline.h"
#include "demcem_constants.h"

void demcem_ws_st_rwg (const double r1[],const double r2[],const double r3[], const double ko, const int Np_1D, complex<double> I_DE[] );

void Example_WS_ST_RWG (const int Counter_start,  const int Counter_end)
{
    
    int N_p;
    

    //
    double tstart, tstop, ttime;
    double RunTime;

    //
    complex<double> I_DE[9];
    
    const double r1[] = {0.0 , 0.0 , 0.5};
    const double r2[] = {0.0 , 0.0 , 0.0};
    const double r3[] = {0.5 , 0.0 , 0.0};
    
    const double ko = 2.0 * M_PI;
    
    //////////////////////////////////////////////////
    //                      ST
    //////////////////////////////////////////////////
    for (int points = Counter_start; points < Counter_end+1; points++)
    {
        
        N_p = points;
        // Time keeper
        tstart = (double)clock()/CLOCKS_PER_SEC;
        
        //
        for (int i=0; i<100; i++)
        {
            demcem_ws_st_rwg ( r1, r2, r3, ko, N_p, I_DE );
            
        } // for
        //
        tstop = (double)clock()/CLOCKS_PER_SEC;
        ttime= tstop-tstart; /*ttime is how long your code run */
        //
        RunTime = ttime * 1000 / 100;
        //
        cout << endl;
        cout << "WS - ST" << endl;
        
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


//
//  inline_functions.h
//
//  Created by Athanasios Polimeridis on 2/17/15.
//  Copyright (c) 2015 Athanasios Polimeridis. All rights reserved.
//

#ifndef _demcem_inline_h
#define _demcem_inline_h

//#include <cmath>
//#include <complex>

//using namespace std;


// **************************************
//			Inline functions
// **************************************

inline
double vector_dot(double x[], double y[]) {
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}
//
inline
void vector_cross(double x[], double y[], double z[])
{
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}


#endif

//
//  main_test.cpp
//  Examples for DEMCEM
//
//  Created by Athanasios Polimeridis on 2/17/15.
//  Copyright (c) 2015 Athanasios Polimeridis. All rights reserved.
//

void Example_WS_ST_RWG (const int Counter_start,  const int Counter_end);
void Example_WS_EA_RWG (const int Counter_start,  const int Counter_end);
void Example_WS_VA_RWG (const int Counter_start,  const int Counter_end);

void Example_SS_EA_RWG (const int Counter_start,  const int Counter_end);
void Example_SS_VA_RWG (const int Counter_start,  const int Counter_end);

void Example_SS_EA_nxRWG (const int Counter_start,  const int Counter_end);
void Example_SS_VA_nxRWG (const int Counter_start,  const int Counter_end);

//
int main(int argc, char *argv[]){

	//
    
    int Np_min = 5;
    int Np_max = 5;
    
    
	//////////////////////////////////////////////////
    //                      WS-ST
    //////////////////////////////////////////////////
	Example_WS_ST_RWG ( Np_min, Np_max );
    
    //////////////////////////////////////////////////
    //                      WS-EA
    //////////////////////////////////////////////////
    Example_WS_EA_RWG ( Np_min, Np_max );
    
    
    //////////////////////////////////////////////////
    //                      WS-VA
    //////////////////////////////////////////////////
    Example_WS_VA_RWG ( Np_min, Np_max );
    
    //////////////////////////////////////////////////
    //                      SS-EA
    //////////////////////////////////////////////////
    Example_SS_EA_RWG ( Np_min, Np_max );
    
    //////////////////////////////////////////////////
    //                      SS-VA
    //////////////////////////////////////////////////
    Example_SS_VA_RWG ( Np_min, Np_max );
    
    //////////////////////////////////////////////////
    //                      SS-EA-nxRWG
    //////////////////////////////////////////////////
    Example_SS_EA_nxRWG ( Np_min, Np_max );
    
    //////////////////////////////////////////////////
    //                      SS-VA-nxRWG
    //////////////////////////////////////////////////
    Example_SS_VA_nxRWG ( Np_min, Np_max );
    
 
    return 0;
}

/**************************************************************************************************************************
           
		           DIRECT_SS_VA_RWG.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the vertex adjacent 4-D
strongly singular integrals over planar triangular elements.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  21 September 2011

  Author:    Athanasios Polimeridis

  References

  A. G. Polimeridis and T. V. Yioultsis, “On the direct evaluation of weakly singular
  integrals in Galerkin mixed potential integral equation formulations,” IEEE Trans.
  Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

  A. G. Polimeridis and J. R. Mosig, “Complete semi-analytical treatment of weakly
  singular integrals on planar triangles via the direct evaluation method,” Int. J.
  Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

  A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, “Fast and accurate
  computation of hyper-singular integrals in Galerkin surface integral equation
  formulations via the direct evaluation method,” IEEE Trans.
  Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

  A. G. Polimeridis and J. R. Mosig, “On the direct evaluation of surface integral
  equation impedance matrix elements involving point singularities,” IEEE Antennas
  Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

  INPUT DATA
  r1,r2,r3,r4,r5 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r1,r2,r5)
  N_theta_p = order of the GL quadrature for the 1-D smooth integral over theta_p
  N_theta_q = order of the GL quadrature for the 1-D smooth integral over theta_q
  N_psi     = order of the GL quadrature for the 1-D smooth integral over Psi
  ko = wavenumber

  OUTPUT DATA
  I_DE[10]  = 4-D strongly singular integrals
              I_DE[0]  = I_f1_f1
              I_DE[1]  = I_f1_f2
              I_DE[2]  = I_f1_f3
              I_DE[3]  = I_f2_f1
              I_DE[4]  = I_f2_f2
              I_DE[5]  = I_f2_f3
              I_DE[6]  = I_f3_f1
              I_DE[7]  = I_f3_f2
              I_DE[8]  = I_f3_f3

**************************************************************************************************************************/

#include "DIRECT_SS_VA_RWG.h"
#include <math.h>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void DIRECT_SS_VA_RWG (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, const int N_theta_p, const int N_theta_q, const int N_psi, const double z_theta_p[], const double w_theta_p[], const double z_theta_q[], const double w_theta_q[], const double z_psi[], const double w_psi[], complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double rp1[3], rp2[3], rp3[3], rp21[3], rp31[3], Ap_cross[3];
	double Ap, Jp;
	double rq1[3], rq2[3], rq3[3], rq21[3], rq31[3], Aq_cross[3];
	double Aq, Jq;
	double r21[3], r31[3], r32[3], r41[3], r54[3], r51[3];
	//
	complex<double> Ising_f1_f2, Ising_f1_f3, Ising_f2_f1, Ising_f2_f2, Ising_f2_f3, Ising_f3_f1, Ising_f3_f2, Ising_f3_f3;
	//
	double theta_p_A, theta_p_B, THETA_p, J_theta_p, Lp;
	double theta_q_A, theta_q_B, THETA_q, J_theta_q, Lq;

	double b_v0[3], b_v1[3];

	double alpha_vo[3], alpha_v1[3], alpha_v2[3], alpha_v3[3];

	double a_v;
	double b_v;
	double c_v;
		
	double psi_1_A, psi_1_B, PSI_1, J_psi_1, psi_2_A, psi_2_B, PSI_2, J_psi_2;

	double B_1;
	double B_2;

	double L1, L2;

	complex<double> K_1[4], K_2[4];

	complex<double> Omega_f1_f2_1, Omega_f1_f3_1, Omega_f2_f1_1, Omega_f2_f2_1, Omega_f2_f3_1, Omega_f3_f1_1, Omega_f3_f2_1, Omega_f3_f3_1;
	complex<double> Omega_f1_f2_2, Omega_f1_f3_2, Omega_f2_f1_2, Omega_f2_f2_2, Omega_f2_f3_2, Omega_f3_f1_2, Omega_f3_f2_2, Omega_f3_f3_2;
	//
	double WPSI, WTHETA_p, WTHETA_q;
	//
	complex<double> I_theta_p_1_1, I_theta_p_1_2, I_theta_p_1_3, I_theta_p_2_1, I_theta_p_2_2, I_theta_p_2_3, I_theta_p_3_1, I_theta_p_3_2, I_theta_p_3_3; 
	complex<double> I_theta_q_1_1, I_theta_q_1_2, I_theta_q_1_3, I_theta_q_2_1, I_theta_q_2_2, I_theta_q_2_3, I_theta_q_3_1, I_theta_q_3_2, I_theta_q_3_3; 
	complex<double> I_psi_1_2, I_psi_1_3, I_psi_2_1, I_psi_2_2, I_psi_2_3, I_psi_3_1, I_psi_3_2, I_psi_3_3; 

	// 2. Coefficients' parameters

	complex<double> coef_f1_f2[2];
	complex<double> coef_f1_f3[4];
	complex<double> coef_f2_f1[2];
	complex<double> coef_f2_f2[5];
	complex<double> coef_f2_f3[7];
	complex<double> coef_f3_f1[4];
	complex<double> coef_f3_f2[7];
	complex<double> coef_f3_f3[8];

	// ************************************************
	//			MAIN CODE
	// ************************************************
	// Evaluate Jacobians for equilateral parameter space
	for (int i = 0; i < 3; i++)
	{
			rp1[i] =  r1[i];
			rp2[i] =  r2[i];
			rp3[i] =  r3[i];
			rp21[i]   = rp2[i] - rp1[i];
	        rp31[i]   = rp3[i] - rp1[i];
	}
	//
    vector_cross(rp21,rp31,Ap_cross);
	Ap = 0.5 * sqrt(vector_dot(Ap_cross,Ap_cross));
	Jp = Ap / sqrt(3.0);
	//
	for (int i = 0; i < 3; i++)
	{
			rq1[i] =  r1[i];
			rq2[i] =  r4[i];
			rq3[i] =  r5[i];
			rq21[i]   = rq2[i] - rq1[i];
	        rq31[i]   = rq3[i] - rq1[i];
	}
	//
    vector_cross(rq21,rq31,Aq_cross);
	Aq = 0.5 * sqrt(vector_dot(Aq_cross,Aq_cross));
	Jq = Aq / sqrt(3.0);
	// Evaluate alpha, beta, gamma r21, r31, r32, r41, r42 parameters
	for (int i = 0; i < 3; i++)
	{
			alpha_vo[i] = (r2[i] - r1[i]) / 2;
			alpha_v1[i] = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			alpha_v2[i] = (r4[i] - r1[i]) / 2;
			alpha_v3[i] = (2 * r5[i] - r1[i] - r4[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r51[i]   = r5[i] - r1[i];
			r54[i]   = r5[i] - r4[i];
	}
	// Get the coefficients
	coefficients_f1_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f1_f2 );
	coefficients_f1_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f1_f3 );
	coefficients_f2_f1 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f1 );
	coefficients_f2_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f2 );
	coefficients_f2_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f2_f3 );
	coefficients_f3_f1 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f1 );
	coefficients_f3_f2 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f2 );
	coefficients_f3_f3 ( r1,  r2,  r3,  r4, r5,  ko,  coef_f3_f3 );
     // Initialization of I_theta_p
		I_theta_p_1_2 = (0.0,0.0); 
		I_theta_p_1_3 = (0.0,0.0);
		I_theta_p_2_1 = (0.0,0.0); 
		I_theta_p_2_2 = (0.0,0.0);
		I_theta_p_2_3 = (0.0,0.0); 
		I_theta_p_3_1 = (0.0,0.0); 
		I_theta_p_3_2 = (0.0,0.0);
		I_theta_p_3_3 = (0.0,0.0);
	 //
	 for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
	 {
		theta_p_A = 0.0;
        theta_p_B = M_PI / 3.0;
        //
        THETA_p = ((theta_p_B - theta_p_A) / 2.0) * z_theta_p[n_theta_p] + (theta_p_B + theta_p_A) / 2.0;
        //
        J_theta_p = (theta_p_B - theta_p_A) / 2.0;
        //
        Lp = (2.0 * sqrt(3.0) ) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p) );
        //
			for (int i = 0; i < 3; i++)
				{
					b_v0[i]   = alpha_vo[i] * cos(THETA_p) + alpha_v1[i] * sin(THETA_p);
				}
		// Initialization of I_theta_q
		I_theta_q_1_2 = (0.0,0.0); 
		I_theta_q_1_3 = (0.0,0.0);
		I_theta_q_2_1 = (0.0,0.0); 
		I_theta_q_2_2 = (0.0,0.0);
		I_theta_q_2_3 = (0.0,0.0); 
		I_theta_q_3_1 = (0.0,0.0); 
		I_theta_q_3_2 = (0.0,0.0);
		I_theta_q_3_3 = (0.0,0.0);
		//
		 for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 {
			theta_q_A = 0.0;
			theta_q_B = M_PI / 3.0;
			//
			THETA_q = ((theta_q_B - theta_q_A) / 2.0) * z_theta_q[n_theta_q] + (theta_q_B + theta_q_A) / 2.0;
			//
			J_theta_q = (theta_q_B - theta_q_A) / 2.0;
			//
			Lq = (2.0 * sqrt(3.0) ) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q) );
			//
			for (int i = 0; i < 3; i++)
				{
					b_v1[i]   = alpha_v2[i] * cos(THETA_q) + alpha_v3[i] * sin(THETA_q);
				}
			//
			a_v = vector_dot(b_v0,b_v0);
			b_v = -2.0 * vector_dot(b_v0,b_v1);
			c_v = vector_dot(b_v1,b_v1);
			// Initialization of I_psi
			I_psi_1_2 = (0.0,0.0); 
			I_psi_1_3 = (0.0,0.0);
			I_psi_2_1 = (0.0,0.0); 
			I_psi_2_2 = (0.0,0.0);
			I_psi_2_3 = (0.0,0.0); 
			I_psi_3_1 = (0.0,0.0); 
			I_psi_3_2 = (0.0,0.0);
			I_psi_3_3 = (0.0,0.0);
			//
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 // psi_A =< PSI <= psi_B
				 psi_1_A = 0.0;
				 psi_1_B = atan(Lq / Lp);
				 //
				 PSI_1 = ((psi_1_B-psi_1_A) / 2.0) * z_psi[n_psi] + (psi_1_B + psi_1_A) / 2.0;
				 //
				 J_psi_1 = (psi_1_B - psi_1_A) / 2.0;
                 //
				 psi_2_A = atan(Lq/Lp);
				 psi_2_B = M_PI / 2.0;
				 //
				 PSI_2 = ((psi_2_B - psi_2_A) / 2.0) * z_psi[n_psi] + (psi_2_B + psi_2_A) / 2.0;
				 //
				 J_psi_2 = (psi_2_B - psi_2_A) / 2.0;
				 //
				 B_1 = sqrt(a_v * pow(cos(PSI_1),2.0) + b_v * cos(PSI_1) * sin(PSI_1) + c_v * pow(sin(PSI_1),2.0));
				 //
				 B_2 = sqrt(a_v * pow(cos(PSI_2),2.0) + b_v * cos(PSI_2) * sin(PSI_2) + c_v * pow(sin(PSI_2),2.0));
				 //
				 Lp = (2.0 * sqrt(3.0)) / (sin(THETA_p) + sqrt(3.0) * cos(THETA_p));
				 L1 = Lp / cos(PSI_1);

				 Lq = (2.0 * sqrt(3.0)) / (sin(THETA_q) + sqrt(3.0) * cos(THETA_q));
				 L2 = Lq / sin(PSI_2);
				 //

				 K_functions(B_1, L1, K_1, ko);
				 K_functions(B_2, L2, K_2, ko);
				 //
				 Omega_f1_f2_1 = Omega_function_f1_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f1_f2, ko, K_1);
				 Omega_f1_f2_2 = Omega_function_f1_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f1_f2, ko, K_2);
				 //
				 Omega_f1_f3_1 = Omega_function_f1_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f1_f3, ko, K_1);
				 Omega_f1_f3_2 = Omega_function_f1_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f1_f3, ko, K_2);
				 //
				 Omega_f2_f1_1 = Omega_function_f2_f1(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f1, ko, K_1);
				 Omega_f2_f1_2 = Omega_function_f2_f1(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f1, ko, K_2);
				 //
				 Omega_f2_f2_1 = Omega_function_f2_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f2, ko, K_1);
				 Omega_f2_f2_2 = Omega_function_f2_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f2, ko, K_2);
				 //
				 Omega_f2_f3_1 = Omega_function_f2_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f2_f3, ko, K_1);
				 Omega_f2_f3_2 = Omega_function_f2_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f2_f3, ko, K_2);
				 //
				 Omega_f3_f1_1 = Omega_function_f3_f1(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f1, ko, K_1);
				 Omega_f3_f1_2 = Omega_function_f3_f1(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f1, ko, K_2);
				 //
				 Omega_f3_f2_1 = Omega_function_f3_f2(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f2, ko, K_1);
				 Omega_f3_f2_2 = Omega_function_f3_f2(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f2, ko, K_2);
				 //
				 Omega_f3_f3_1 = Omega_function_f3_f3(THETA_p, THETA_q, PSI_1, B_1, coef_f3_f3, ko, K_1);
				 Omega_f3_f3_2 = Omega_function_f3_f3(THETA_p, THETA_q, PSI_2, B_2, coef_f3_f3, ko, K_2);
				 //
         		 WPSI = w_psi[n_psi];
				 //
				 I_psi_1_2   +=  WPSI * (J_psi_1 * Omega_f1_f2_1 + J_psi_2 * Omega_f1_f2_2);
				 I_psi_1_3   +=  WPSI * (J_psi_1 * Omega_f1_f3_1 + J_psi_2 * Omega_f1_f3_2);
				 I_psi_2_1   +=  WPSI * (J_psi_1 * Omega_f2_f1_1 + J_psi_2 * Omega_f2_f1_2);
				 I_psi_2_2   +=  WPSI * (J_psi_1 * Omega_f2_f2_1 + J_psi_2 * Omega_f2_f2_2);
				 I_psi_2_3   +=  WPSI * (J_psi_1 * Omega_f2_f3_1 + J_psi_2 * Omega_f2_f3_2);
				 I_psi_3_1   +=  WPSI * (J_psi_1 * Omega_f3_f1_1 + J_psi_2 * Omega_f3_f1_2);
				 I_psi_3_2   +=  WPSI * (J_psi_1 * Omega_f3_f2_1 + J_psi_2 * Omega_f3_f2_2);
				 I_psi_3_3   +=  WPSI * (J_psi_1 * Omega_f3_f3_1 + J_psi_2 * Omega_f3_f3_2);
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 WTHETA_q = w_theta_q[n_theta_q];
			 //
			 I_theta_q_1_2     +=  WTHETA_q * I_psi_1_2;
			 I_theta_q_1_3     +=  WTHETA_q * I_psi_1_3;
			 I_theta_q_2_1     +=  WTHETA_q * I_psi_2_1;
			 I_theta_q_2_2     +=  WTHETA_q * I_psi_2_2;
			 I_theta_q_2_3     +=  WTHETA_q * I_psi_2_3;
			 I_theta_q_3_1     +=  WTHETA_q * I_psi_3_1;
			 I_theta_q_3_2     +=  WTHETA_q * I_psi_3_2;
			 I_theta_q_3_3     +=  WTHETA_q * I_psi_3_3;
			 
		 } //end for ( int n_theta_q = 0 ; n_theta_q <  N_theta_q ; n_theta_q++ )
		 //
		 I_theta_q_1_2     *=  J_theta_q ;
		 I_theta_q_1_3     *=  J_theta_q ;
		 I_theta_q_2_1     *=  J_theta_q ;
		 I_theta_q_2_2     *=  J_theta_q ;
		 I_theta_q_2_3     *=  J_theta_q ;
		 I_theta_q_3_1     *=  J_theta_q ;
		 I_theta_q_3_2     *=  J_theta_q ;
		 I_theta_q_3_3     *=  J_theta_q ;
		 //
		 WTHETA_p = w_theta_p[n_theta_p];
		 //
		 I_theta_p_1_2    += WTHETA_p * I_theta_q_1_2;
		 I_theta_p_1_3    += WTHETA_p * I_theta_q_1_3;
		 I_theta_p_2_1    += WTHETA_p * I_theta_q_2_1;
		 I_theta_p_2_2    += WTHETA_p * I_theta_q_2_2;
		 I_theta_p_2_3    += WTHETA_p * I_theta_q_2_3;
		 I_theta_p_3_1    += WTHETA_p * I_theta_q_3_1;
		 I_theta_p_3_2    += WTHETA_p * I_theta_q_3_2;
		 I_theta_p_3_3    += WTHETA_p * I_theta_q_3_3;

	 } //end for ( int n_theta_p = 0 ; n_theta_p <  N_theta_p ; n_theta_p++ )
	 //
	 I_theta_p_1_2     *=  J_theta_p ;
	 I_theta_p_1_3     *=  J_theta_p ;
	 I_theta_p_2_1     *=  J_theta_p ;
	 I_theta_p_2_2     *=  J_theta_p ;
	 I_theta_p_2_3     *=  J_theta_p ;
	 I_theta_p_3_1     *=  J_theta_p ;
	 I_theta_p_3_2     *=  J_theta_p ;
	 I_theta_p_3_3     *=  J_theta_p ;
	 //
	 double PRECOMP = Jp * Jq / (4 * Ap * Aq);
	 // FINAL OUTPUT
	 I_DE[0] = (0.0, 0.0);
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r51,r51)) * I_theta_p_1_2;
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * I_theta_p_1_3;
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r54,r54)) * I_theta_p_2_1;
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r51,r51)) * I_theta_p_2_2;
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * I_theta_p_2_3;
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r54,r54)) * I_theta_p_3_1;
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r51,r51)) * I_theta_p_3_2;
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * I_theta_p_3_3;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f2
// ***********************************************************************

void coefficients_f1_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[0] * ee[1];
	complex<double> t4 = bb[0] * ee[2];
	complex<double> t7 = bb[1] * ee[2];
	complex<double> t10 = bb[1] * ee[0];
	complex<double> t13 = bb[2] * ee[0];
	complex<double> t16 = bb[2] * ee[1];
	complex<double> t37 = -t1 * dd[2] / 0.6e1 + t4 * dd[1] / 0.6e1 - t7 * dd[0] / 0.6e1 + t10 * dd[2] / 0.6e1 - t13 * dd[1] / 0.6e1 + t16 * dd[0] / 0.6e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[1] * ee[2] * dd[0] / 0.3e1 - cc[1] * ee[0] * dd[2] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 - cc[2] * ee[1] * dd[0] / 0.3e1;
	complex<double> t38 = sqrt(0.3e1);
	complex<double> t39 = dd[2] * t38;
	complex<double> t41 = dd[1] * t38;
	complex<double> t43 = dd[0] * t38;
	//
	coef[0] = t37;
	coef[1] = t1 * t39 / 0.6e1 - t4 * t41 / 0.6e1 + t7 * t43 / 0.6e1 - t10 * t39 / 0.6e1 + t13 * t41 / 0.6e1 - t16 * t43 / 0.6e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f3
// ***********************************************************************

void coefficients_f1_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = bb[0] * ee[1];
	complex<double> t2 = sqrt(0.3e1);
	complex<double> t3 = dd[2] * t2;
	complex<double> t4 = t1 * t3;
	complex<double> t5 = bb[0] * dd[1];
	complex<double> t6 = ee[2] * t2;
	complex<double> t7 = t5 * t6;
	complex<double> t8 = bb[1] * ee[2];
	complex<double> t9 = dd[0] * t2;
	complex<double> t10 = t8 * t9;
	complex<double> t11 = bb[1] * dd[2];
	complex<double> t12 = ee[0] * t2;
	complex<double> t13 = t11 * t12;
	complex<double> t14 = bb[2] * ee[0];
	complex<double> t15 = dd[1] * t2;
	complex<double> t16 = t14 * t15;
	complex<double> t17 = bb[2] * dd[0];
	complex<double> t18 = ee[1] * t2;
	complex<double> t19 = t17 * t18;
	complex<double> t21 = cc[1] * dd[2];
	complex<double> t24 = cc[0] * ee[1];
	complex<double> t27 = t11 * ee[0];
	complex<double> t29 = t14 * dd[1];
	complex<double> t31 = t17 * ee[1];
	complex<double> t33 = t5 * ee[2];
	complex<double> t35 = t8 * dd[0];
	complex<double> t37 = t1 * dd[2];
	complex<double> t39 = cc[1] * ee[2];
	complex<double> t42 = cc[2] * ee[0];
	complex<double> t45 = cc[2] * dd[0];
	complex<double> t48 = cc[0] * dd[1];
	complex<double> t51 = -t21 * ee[0] / 0.6e1 + t24 * dd[2] / 0.6e1 + t27 / 0.12e2 - t29 / 0.12e2 + t31 / 0.12e2 + t33 / 0.12e2 - t35 / 0.12e2 - t37 / 0.12e2 + t39 * dd[0] / 0.6e1 + t42 * dd[1] / 0.6e1 - t45 * ee[1] / 0.6e1 - t48 * ee[2] / 0.6e1;
	complex<double> t71 = t4 / 0.12e2 - t7 / 0.12e2 + t10 / 0.12e2 - t13 / 0.12e2 - t24 * t3 / 0.6e1 + t16 / 0.12e2 - t39 * t9 / 0.6e1 - t19 / 0.12e2 - t42 * t15 / 0.6e1 + t48 * t6 / 0.6e1 + t45 * t18 / 0.6e1 + t21 * t12 / 0.6e1;
	//
	coef[0] = t4 / 0.12e2 - t7 / 0.12e2 + t10 / 0.12e2 - t13 / 0.12e2 + t16 / 0.12e2 - t19 / 0.12e2;
	coef[1] = t51;
	coef[2] = t33 / 0.4e1 + t31 / 0.4e1 - t37 / 0.4e1 - t35 / 0.4e1 - t29 / 0.4e1 + t27 / 0.4e1;
	coef[3] = t71;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f1
// ***********************************************************************

void coefficients_f2_f1 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double aa[3], bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	complex<double> j   = Iunit;
	//
	complex<double> t1 = bb[0] * cc[2];
	complex<double> t4 = bb[0] * cc[1];
	complex<double> t7 = bb[1] * cc[0];
	complex<double> t10 = bb[1] * cc[2];
	complex<double> t13 = bb[2] * cc[1];
	complex<double> t16 = bb[2] * cc[0];
	complex<double> t31 = -t1 * dd[1] / 0.6e1 + t4 * dd[2] / 0.6e1 - t7 * dd[2] / 0.6e1 + t10 * dd[0] / 0.6e1 - t13 * dd[0] / 0.6e1 + t16 * dd[1] / 0.6e1 - t4 * ee[2] / 0.3e1 + t1 * ee[1] / 0.3e1 - t10 * ee[0] / 0.3e1 + t7 * ee[2] / 0.3e1 - t16 * ee[1] / 0.3e1 + t13 * ee[0] / 0.3e1;
	complex<double> t32 = sqrt(0.3e1);
	complex<double> t33 = dd[1] * t32;
	complex<double> t35 = dd[2] * t32;
	complex<double> t38 = dd[0] * t32;
	//
	coef[0] = t31;
	coef[1] = t1 * t33 / 0.6e1 - t4 * t35 / 0.6e1 + t7 * t35 / 0.6e1 - t10 * t38 / 0.6e1 + t13 * t38 / 0.6e1 - t16 * t33 / 0.6e1;
	}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f2
// ***********************************************************************

void coefficients_f2_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t7 = bb[0] * cc[2];
	complex<double> t10 = bb[1] * cc[2];
	complex<double> t19 = bb[1] * cc[0];
	complex<double> t25 = bb[0] * cc[1];
	complex<double> t31 = bb[2] * cc[0];
	complex<double> t34 = bb[2] * cc[1];
	complex<double> t37 = bb[1] * ee[0] * dd[2] / 0.6e1 - bb[2] * ee[0] * dd[1] / 0.6e1 - t7 * dd[1] / 0.6e1 + t10 * dd[0] / 0.6e1 - bb[1] * ee[2] * dd[0] / 0.6e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - t19 * dd[2] / 0.6e1 + cc[1] * ee[2] * dd[0] / 0.3e1 + t25 * dd[2] / 0.6e1 - bb[0] * ee[1] * dd[2] / 0.6e1 + t31 * dd[1] / 0.6e1 - t34 * dd[0] / 0.6e1;
	complex<double> t68 = bb[0] * ee[2] * dd[1] / 0.6e1 + bb[2] * ee[1] * dd[0] / 0.6e1 - t25 * ee[2] / 0.3e1 - t31 * ee[1] / 0.3e1 + t34 * ee[0] / 0.3e1 - cc[1] * ee[0] * dd[2] / 0.3e1 - t10 * ee[0] / 0.3e1 + t19 * ee[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 - cc[2] * ee[1] * dd[0] / 0.3e1 + t7 * ee[1] / 0.3e1;
	complex<double> t70 = sqrt(0.3e1);
	complex<double> t71 = t70 * bb[0];
	complex<double> t76 = t70 * bb[1];
	complex<double> t81 = t70 * bb[2];
	complex<double> t86 = t71 * cc[1] * dd[2] - t71 * cc[2] * dd[1] + t76 * cc[2] * dd[0] - t76 * cc[0] * dd[2] + t81 * cc[0] * dd[1] - t81 * cc[1] * dd[0];
	complex<double> t99 = -t71 * ee[1] * dd[2] + t71 * ee[2] * dd[1] - t76 * ee[2] * dd[0] + t76 * ee[0] * dd[2] - t81 * ee[0] * dd[1] + t81 * ee[1] * dd[0];
	//
	coef[0] = t37 + t68;
	coef[1] = t86 / 0.3e1;
	coef[2] = t99 / 0.3e1;
	coef[3] = -t99 / 0.6e1;
	coef[4] = -t86 / 0.6e1;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f3
// ***********************************************************************

void coefficients_f2_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t2 = bb[2] * dd[0] * ee[1];
	complex<double> t5 = bb[2] * ee[0] * dd[1];
	complex<double> t8 = bb[0] * ee[1] * dd[2];
	complex<double> t10 = bb[2] * cc[0];
	complex<double> t13 = bb[0] * cc[1];
	complex<double> t16 = cc[2] * ee[0];
	complex<double> t19 = cc[0] * ee[1];
	complex<double> t22 = cc[1] * ee[2];
	complex<double> t26 = bb[1] * ee[2] * dd[0];
	complex<double> t29 = bb[1] * dd[2] * ee[0];
	complex<double> t32 = bb[0] * dd[1] * ee[2];
	complex<double> t34 = bb[1] * cc[0];
	complex<double> t37 = t2 / 0.12e2 - t5 / 0.12e2 - t8 / 0.12e2 + t10 * dd[1] / 0.6e1 + t13 * dd[2] / 0.6e1 + t16 * dd[1] / 0.6e1 + t19 * dd[2] / 0.6e1 + t22 * dd[0] / 0.6e1 - t26 / 0.12e2 + t29 / 0.12e2 + t32 / 0.12e2 - t34 * dd[2] / 0.6e1;
	complex<double> t38 = bb[1] * cc[2];
	complex<double> t41 = cc[1] * dd[2];
	complex<double> t44 = cc[0] * dd[1];
	complex<double> t47 = bb[2] * cc[1];
	complex<double> t50 = bb[0] * cc[2];
	complex<double> t55 = cc[2] * dd[0];
	complex<double> t68 = t38 * dd[0] / 0.6e1 - t41 * ee[0] / 0.6e1 - t44 * ee[2] / 0.6e1 - t47 * dd[0] / 0.6e1 - t50 * dd[1] / 0.6e1 + t47 * ee[0] / 0.3e1 - t55 * ee[1] / 0.6e1 - t13 * ee[2] / 0.3e1 + t50 * ee[1] / 0.3e1 - t10 * ee[1] / 0.3e1 + t34 * ee[2] / 0.3e1 - t38 * ee[0] / 0.3e1;
	complex<double> t70 = sqrt(0.3e1);
	complex<double> t71 = t70 * bb[0];
	complex<double> t73 = t71 * ee[1] * dd[2];
	complex<double> t75 = t71 * dd[1] * ee[2];
	complex<double> t76 = t70 * bb[1];
	complex<double> t78 = t76 * ee[2] * dd[0];
	complex<double> t80 = t76 * dd[2] * ee[0];
	complex<double> t81 = t70 * bb[2];
	complex<double> t83 = t81 * ee[0] * dd[1];
	complex<double> t85 = t81 * dd[0] * ee[1];
	complex<double> t86 = -t73 + t75 - t78 + t80 - t83 + t85;
	complex<double> t97 = t8 - t32 + t26 - t29 + t5 - t2;
	complex<double> t98 = dd[2] * t70;
	complex<double> t108 = dd[1] * t70;
	complex<double> t116 = dd[0] * t70;
	complex<double> t134 = -t13 * t98 / 0.6e1 - t19 * t98 / 0.6e1 + t73 / 0.12e2 + t34 * t98 / 0.6e1 - t75 / 0.12e2 + t78 / 0.12e2 - t10 * t108 / 0.6e1 + t41 * ee[0] * t70 / 0.6e1 - t16 * t108 / 0.6e1 - t38 * t116 / 0.6e1 - t22 * t116 / 0.6e1 + t44 * ee[2] * t70 / 0.6e1 + t83 / 0.12e2 + t50 * t108 / 0.6e1 - t85 / 0.12e2 + t55 * ee[1] * t70 / 0.6e1 - t80 / 0.12e2 + t47 * t116 / 0.6e1;
	//
	coef[0] = t37 + t68;
	coef[1] = t86 / 0.6e1;
	coef[2] = t71 * t22 / 0.3e1 - t71 * cc[2] * ee[1] / 0.3e1 + t76 * t16 / 0.3e1 - t76 * cc[0] * ee[2] / 0.3e1 + t81 * t19 / 0.3e1 - t81 * cc[1] * ee[0] / 0.3e1;
	coef[3] = t97 / 0.2e1;
	coef[4] = -t97 / 0.4e1;
	coef[5] = -t86 / 0.12e2;
	coef[6] = t134;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f1
// ***********************************************************************

void coefficients_f3_f1 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = cc[0] * bb[1];
	complex<double> t2 = sqrt(0.3e1);
	complex<double> t3 = dd[2] * t2;
	complex<double> t4 = t1 * t3;
	complex<double> t6 = cc[0] * bb[2];
	complex<double> t7 = ee[1] * t2;
	complex<double> t10 = cc[2] * bb[1];
	complex<double> t11 = dd[0] * t2;
	complex<double> t12 = t10 * t11;
	complex<double> t14 = dd[1] * t2;
	complex<double> t15 = t6 * t14;
	complex<double> t17 = ee[2] * t2;
	complex<double> t20 = cc[1] * bb[2];
	complex<double> t21 = t20 * t11;
	complex<double> t23 = ee[0] * t2;
	complex<double> t26 = cc[1] * bb[0];
	complex<double> t27 = t26 * t3;
	complex<double> t29 = cc[2] * bb[0];
	complex<double> t32 = t29 * t14;
	complex<double> t38 = t4 / 0.12e2 + t6 * t7 / 0.6e1 - t12 / 0.12e2 - t15 / 0.12e2 - t1 * t17 / 0.6e1 + t21 / 0.12e2 - t20 * t23 / 0.6e1 - t27 / 0.12e2 - t29 * t7 / 0.6e1 + t32 / 0.12e2 + t10 * t23 / 0.6e1 + t26 * t17 / 0.6e1;
	complex<double> t39 = t29 * dd[1];
	complex<double> t41 = t1 * dd[2];
	complex<double> t45 = t6 * dd[1];
	complex<double> t49 = t20 * dd[0];
	complex<double> t57 = t10 * dd[0];
	complex<double> t59 = t26 * dd[2];
	complex<double> t63 = -t39 / 0.12e2 - t41 / 0.12e2 + t1 * ee[2] / 0.6e1 + t45 / 0.12e2 - t6 * ee[1] / 0.6e1 - t49 / 0.12e2 - t10 * ee[0] / 0.6e1 + t29 * ee[1] / 0.6e1 + t20 * ee[0] / 0.6e1 + t57 / 0.12e2 + t59 / 0.12e2 - t26 * ee[2] / 0.6e1;
	//
	coef[0] = t38;
	coef[1] = t63;
	coef[2] = -t39 / 0.4e1 + t57 / 0.4e1 + t45 / 0.4e1 + t59 / 0.4e1 - t49 / 0.4e1 - t41 / 0.4e1;
	coef[3] = -t27 / 0.12e2 + t21 / 0.12e2 + t4 / 0.12e2 - t12 / 0.12e2 - t15 / 0.12e2 + t32 / 0.12e2;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f2
// ***********************************************************************

void coefficients_f3_f2 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t4 = cc[0] * bb[1];
	complex<double> t5 = t4 * dd[2];
	complex<double> t7 = cc[0] * bb[2];
	complex<double> t8 = t7 * dd[1];
	complex<double> t10 = cc[2] * bb[1];
	complex<double> t11 = t10 * dd[0];
	complex<double> t13 = bb[0] * ee[1];
	complex<double> t16 = cc[1] * bb[0];
	complex<double> t17 = t16 * dd[2];
	complex<double> t21 = bb[1] * ee[0];
	complex<double> t26 = bb[2] * ee[1];
	complex<double> t29 = bb[0] * ee[2];
	complex<double> t32 = bb[2] * ee[0];
	complex<double> t35 = -cc[2] * ee[1] * dd[0] / 0.3e1 - t5 / 0.12e2 + t8 / 0.12e2 + t11 / 0.12e2 - t13 * dd[2] / 0.6e1 + t17 / 0.12e2 + t4 * ee[2] / 0.6e1 + t21 * dd[2] / 0.6e1 - t16 * ee[2] / 0.6e1 + t26 * dd[0] / 0.6e1 + t29 * dd[1] / 0.6e1 - t32 * dd[1] / 0.6e1;
	complex<double> t36 = cc[1] * bb[2];
	complex<double> t37 = t36 * dd[0];
	complex<double> t43 = cc[2] * bb[0];
	complex<double> t44 = t43 * dd[1];
	complex<double> t65 = bb[1] * ee[2];
	complex<double> t68 = -t37 / 0.12e2 + t36 * ee[0] / 0.6e1 - t7 * ee[1] / 0.6e1 - t44 / 0.12e2 - t10 * ee[0] / 0.6e1 - cc[1] * ee[0] * dd[2] / 0.3e1 + cc[0] * ee[1] * dd[2] / 0.3e1 - cc[0] * ee[2] * dd[1] / 0.3e1 + cc[2] * ee[0] * dd[1] / 0.3e1 + cc[1] * ee[2] * dd[0] / 0.3e1 + t43 * ee[1] / 0.6e1 - t65 * dd[0] / 0.6e1;
	complex<double> t70 = sqrt(0.3e1);
	complex<double> t71 = t70 * cc[0];
	complex<double> t73 = t71 * bb[1] * dd[2];
	complex<double> t75 = t71 * bb[2] * dd[1];
	complex<double> t76 = t70 * cc[1];
	complex<double> t78 = t76 * bb[2] * dd[0];
	complex<double> t80 = t76 * bb[0] * dd[2];
	complex<double> t81 = t70 * cc[2];
	complex<double> t83 = t81 * bb[0] * dd[1];
	complex<double> t85 = t81 * bb[1] * dd[0];
	complex<double> t86 = -t73 + t75 - t78 + t80 - t83 + t85;
	complex<double> t87 = t5 - t8 + t37 - t17 + t44 - t11;
	complex<double> t101 = dd[2] * t70;
	complex<double> t104 = ee[1] * t70;
	complex<double> t110 = dd[0] * t70;
	complex<double> t115 = ee[2] * t70;
	complex<double> t122 = dd[1] * t70;
	complex<double> t125 = ee[0] * t70;
	complex<double> t137 = t13 * t101 / 0.6e1 + t7 * t104 / 0.6e1 - t43 * t104 / 0.6e1 + t73 / 0.12e2 + t65 * t110 / 0.6e1 - t85 / 0.12e2 - t80 / 0.12e2 - t4 * t115 / 0.6e1 - t75 / 0.12e2 - t21 * t101 / 0.6e1 + t83 / 0.12e2 - t29 * t122 / 0.6e1 - t36 * t125 / 0.6e1 + t32 * t122 / 0.6e1 - t26 * t110 / 0.6e1 + t10 * t125 / 0.6e1 + t78 / 0.12e2 + t16 * t115 / 0.6e1;
	//
	coef[0] = t35 + t68;
	coef[1] = t86 / 0.6e1;
	coef[2] = t87 / 0.2e1;
	coef[3] = -t71 * ee[1] * dd[2] / 0.3e1 + t71 * ee[2] * dd[1] / 0.3e1 - t76 * ee[2] * dd[0] / 0.3e1 + t76 * ee[0] * dd[2] / 0.3e1 - t81 * ee[0] * dd[1] / 0.3e1 + t81 * ee[1] * dd[0] / 0.3e1;
	coef[4] = -t87 / 0.4e1;
	coef[5] = t137;
	coef[6] = -t86 / 0.12e2;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f3
// ***********************************************************************

void coefficients_f3_f3 (const double r1[],const double r2[],const double r3[],const double r4[], const double r5[], const double ko, complex<double> coef[] )
{                         
	double bb[3], cc[3], dd[3], ee[3];
	//
	for (int i = 0; i < 3; i++)
	{
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
			ee[i] = r5[i] - r1[i];
	}
	//
	complex<double> t1 = sqrt(0.3e1);
	complex<double> t2 = t1 * cc[0];
	complex<double> t3 = bb[1] * ee[2];
	complex<double> t4 = t2 * t3;
	complex<double> t6 = t2 * bb[2] * ee[1];
	complex<double> t7 = t1 * cc[1];
	complex<double> t8 = bb[2] * ee[0];
	complex<double> t9 = t7 * t8;
	complex<double> t11 = t7 * bb[0] * ee[2];
	complex<double> t12 = t1 * cc[2];
	complex<double> t13 = bb[0] * ee[1];
	complex<double> t14 = t12 * t13;
	complex<double> t16 = t12 * bb[1] * ee[0];
	complex<double> t19 = cc[0] * dd[1] * ee[2];
	complex<double> t21 = cc[1] * ee[2] * dd[0];
	complex<double> t23 = cc[0] * ee[1] * dd[2];
	complex<double> t25 = cc[2] * dd[0] * ee[1];
	complex<double> t27 = cc[1] * dd[2] * ee[0];
	complex<double> t29 = cc[2] * ee[0] * dd[1];
	complex<double> t32 = t2 * ee[1] * dd[2];
	complex<double> t34 = t2 * dd[1] * ee[2];
	complex<double> t36 = t7 * ee[2] * dd[0];
	complex<double> t38 = t7 * dd[2] * ee[0];
	complex<double> t40 = t12 * ee[0] * dd[1];
	complex<double> t42 = t12 * dd[0] * ee[1];
	complex<double> t44 = cc[1] * bb[2];
	complex<double> t45 = t44 * ee[0];
	complex<double> t46 = cc[1] * bb[0];
	complex<double> t47 = t46 * ee[2];
	complex<double> t48 = cc[2] * bb[0];
	complex<double> t49 = t48 * ee[1];
	complex<double> t50 = cc[2] * bb[1];
	complex<double> t51 = t50 * ee[0];
	complex<double> t52 = cc[0] * bb[2];
	complex<double> t53 = t52 * ee[1];
	complex<double> t54 = cc[0] * bb[1];
	complex<double> t55 = t54 * ee[2];
	complex<double> t57 = dd[2] * t1;
	complex<double> t59 = t54 * t57 / 0.12e2;
	complex<double> t60 = bb[2] * dd[0];
	complex<double> t63 = t60 * ee[1] * t1 / 0.12e2;
	complex<double> t64 = dd[1] * t1;
	complex<double> t66 = t8 * t64 / 0.12e2;
	complex<double> t68 = t48 * t64 / 0.12e2;
	complex<double> t69 = dd[0] * t1;
	complex<double> t71 = t44 * t69 / 0.12e2;
	complex<double> t73 = t3 * t69 / 0.12e2;
	complex<double> t75 = t52 * t64 / 0.12e2;
	complex<double> t76 = bb[1] * dd[2];
	complex<double> t79 = t76 * ee[0] * t1 / 0.12e2;
	complex<double> t81 = t46 * t57 / 0.12e2;
	complex<double> t83 = t50 * t69 / 0.12e2;
	complex<double> t85 = t13 * t57 / 0.12e2;
	complex<double> t88 = bb[0] * dd[1];
	complex<double> t91 = t88 * ee[2] * t1 / 0.12e2;
	complex<double> t96 = t59 - t63 + t66 + t68 + t71 + t73 - t75 - t79 - t81 - t83 + t85 - t40 / 0.6e1 - t36 / 0.6e1 - t91 + t38 / 0.6e1 - t32 / 0.6e1 + t42 / 0.6e1 + t34 / 0.6e1;
	complex<double> t97 = t54 * dd[2];
	complex<double> t98 = t3 * dd[0];
	complex<double> t99 = t88 * ee[2];
	complex<double> t100 = t13 * dd[2];
	complex<double> t101 = t50 * dd[0];
	complex<double> t102 = t44 * dd[0];
	complex<double> t103 = t8 * dd[1];
	complex<double> t104 = t76 * ee[0];
	complex<double> t105 = t46 * dd[2];
	complex<double> t106 = t52 * dd[1];
	complex<double> t107 = t60 * ee[1];
	complex<double> t108 = t48 * dd[1];
	complex<double> t109 = -t97 - t98 + t99 - t100 + t101 - t102 - t103 + t104 + t105 + t106 + t107 - t108;
	complex<double> t122 = t104 / 0.12e2 - t98 / 0.12e2 - t102 / 0.12e2 + t99 / 0.12e2 + t105 / 0.12e2 - t100 / 0.12e2 - t103 / 0.12e2 - t108 / 0.12e2 + t107 / 0.12e2 - t97 / 0.12e2 + t106 / 0.12e2 + t49 / 0.6e1;
	complex<double> t135 = t101 / 0.12e2 - t27 / 0.6e1 - t47 / 0.6e1 + t23 / 0.6e1 + t55 / 0.6e1 - t25 / 0.6e1 - t51 / 0.6e1 - t19 / 0.6e1 + t45 / 0.6e1 + t29 / 0.6e1 - t53 / 0.6e1 + t21 / 0.6e1;
	complex<double> t143 = t85 - t63 + t66 + t6 / 0.6e1 - t91 - t83 - t4 / 0.6e1 + t71 + t11 / 0.6e1 + t59 - t9 / 0.6e1 - t81 + t16 / 0.6e1 + t68 - t14 / 0.6e1 - t75 - t79 + t73;
	//
	coef[0] = -t4 / 0.6e1 + t6 / 0.6e1 - t9 / 0.6e1 + t11 / 0.6e1 - t14 / 0.6e1 + t16 / 0.6e1;
	coef[1] = -t19 / 0.2e1 + t21 / 0.2e1 + t23 / 0.2e1 - t25 / 0.2e1 - t27 / 0.2e1 + t29 / 0.2e1;
	coef[2] = -t32 / 0.6e1 + t34 / 0.6e1 - t36 / 0.6e1 + t38 / 0.6e1 - t40 / 0.6e1 + t42 / 0.6e1;
	coef[3] = t45 / 0.2e1 - t47 / 0.2e1 + t49 / 0.2e1 - t51 / 0.2e1 - t53 / 0.2e1 + t55 / 0.2e1;
	coef[4] = t96;
	coef[5] = t109 / 0.4e1;
	coef[6] = t122 + t135;
	coef[7] = t143;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f1_f2
// ***********************************************************************

complex<double> Omega_function_f1_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t4 * t4;
	complex<double> t6 = t3 * t5;
	complex<double> t7 = pow(GAMMA, 0.2e1);
	complex<double> t10 = t6 / t7 / GAMMA;
	complex<double> t11 = sin(theta_p);
	complex<double> t12 = sin(theta_q);
	complex<double> t14 = coef[0];
	complex<double> t17 = cos(theta_p);
	complex<double> t19 = coef[1];
	complex<double> t27 = t6 / t7 * j;
	//
	X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t17 * t12 * t19) + K[2] * (-t27 * ko * t11 * t12 * t14 - t27 * ko * t17 * t12 * t19);	
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f1_f3
// ***********************************************************************

complex<double> Omega_function_f1_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t4 * t4;
	complex<double> t6 = t3 * t5;
	complex<double> t7 = pow(GAMMA, 0.2e1);
	complex<double> t10 = t6 / t7 / GAMMA;
	complex<double> t11 = sin(theta_q);
	complex<double> t12 = cos(theta_p);
	complex<double> t14 = coef[0];
	complex<double> t17 = sin(theta_p);
	complex<double> t19 = coef[1];
	complex<double> t22 = cos(theta_q);
	complex<double> t24 = coef[2];
	complex<double> t28 = coef[3];
	complex<double> t36 = t6 / t7 * j;
	complex<double> t37 = ko * t11;
	complex<double> t44 = ko * t22;
	//
	X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t11 * t17 * t19 - t10 * t22 * t12 * t24 - t10 * t22 * t17 * t28) + K[2] * (-t36 * t37 * t12 * t14 - t36 * t37 * t17 * t19 - t36 * t44 * t12 * t24 - t36 * t44 * t17 * t28);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f1
// ***********************************************************************

complex<double> Omega_function_f2_f1 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t4 * t4;
	complex<double> t6 = t3 * t5;
	complex<double> t7 = pow(GAMMA, 0.2e1);
	complex<double> t10 = t6 / t7 / GAMMA;
	complex<double> t11 = sin(theta_p);
	complex<double> t12 = sin(theta_q);
	complex<double> t14 = coef[0];
	complex<double> t17 = cos(theta_q);
	complex<double> t19 = coef[1];
	complex<double> t27 = t6 / t7 * j;
	//
	X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t17 * t11 * t19) + K[2] * (-t27 * ko * t11 * t12 * t14 - t27 * ko * t17 * t11 * t19);
// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f2
// ***********************************************************************

complex<double> Omega_function_f2_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = pow(GAMMA, 0.2e1);
	complex<double> t8 = 0.1e1 / t6 / GAMMA;
	complex<double> t9 = sin(theta_p);
	complex<double> t11 = coef[1];
	complex<double> t14 = t2 * t2;
	complex<double> t15 = t14 * t3;
	complex<double> t16 = sin(theta_q);
	complex<double> t18 = coef[2];
	complex<double> t24 = t14 * t4;
	complex<double> t25 = t24 * t8;
	complex<double> t27 = coef[0];
	complex<double> t30 = cos(theta_p);
	complex<double> t32 = coef[3];
	complex<double> t35 = cos(theta_q);
	complex<double> t37 = coef[4];
	complex<double> t40 = 0.1e1 / t6;
	complex<double> t42 = j * ko;
	complex<double> t54 = t24 * t40 * j;
	//
	X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18) + K[1] * (-t25 * t9 * t16 * t27 - t25 * t30 * t16 * t32 - t25 * t35 * t9 * t37 - t5 * t40 * t42 * t9 * t11 - t15 * t40 * t42 * t16 * t18) + K[2] * (-t54 * ko * t9 * t16 * t27 - t54 * ko * t30 * t16 * t32 - t54 * ko * t35 * t9 * t37);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f2_f3
// ***********************************************************************

complex<double> Omega_function_f2_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = pow(GAMMA, 0.2e1);
	complex<double> t8 = 0.1e1 / t6 / GAMMA;
	complex<double> t9 = sin(theta_q);
	complex<double> t11 = coef[1];
	complex<double> t14 = t4 * t4;
	complex<double> t15 = t2 * t14;
	complex<double> t16 = sin(theta_p);
	complex<double> t18 = coef[2];
	complex<double> t21 = cos(theta_q);
	complex<double> t23 = coef[3];
	complex<double> t29 = t3 * t14;
	complex<double> t30 = t29 * t8;
	complex<double> t32 = coef[0];
	complex<double> t35 = cos(theta_p);
	complex<double> t37 = coef[4];
	complex<double> t41 = coef[5];
	complex<double> t45 = coef[6];
	complex<double> t48 = 0.1e1 / t6;
	complex<double> t49 = t5 * t48;
	complex<double> t50 = j * ko;
	complex<double> t65 = t29 * t48 * j;
	complex<double> t70 = ko * t21;
	//
	X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18 - t5 * t8 * t21 * t23) + K[1] * (-t30 * t16 * t9 * t32 - t30 * t21 * t35 * t37 - t30 * t9 * t35 * t41 - t30 * t21 * t16 * t45 - t49 * t50 * t9 * t11 - t15 * t48 * t50 * t16 * t18 - t49 * t50 * t21 * t23) + K[2] * (-t65 * ko * t16 * t9 * t32 - t65 * t70 * t35 * t37 - t65 * ko * t9 * t35 * t41 - t65 * t70 * t16 * t45);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f1
// ***********************************************************************

complex<double> Omega_function_f3_f1 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = t2 * t2;
	complex<double> t4 = cos(Psi);
	complex<double> t5 = t4 * t4;
	complex<double> t6 = t3 * t5;
	complex<double> t7 = pow(GAMMA, 0.2e1);
	complex<double> t10 = t6 / t7 / GAMMA;
	complex<double> t11 = sin(theta_q);
	complex<double> t12 = cos(theta_p);
	complex<double> t14 = coef[0];
	complex<double> t17 = sin(theta_p);
	complex<double> t19 = coef[1];
	complex<double> t22 = cos(theta_q);
	complex<double> t24 = coef[2];
	complex<double> t28 = coef[3];
	complex<double> t36 = t6 / t7 * j;
	complex<double> t37 = ko * t11;
	complex<double> t44 = ko * t22;
	//
	X = K[1] * (-t10 * t11 * t12 * t14 - t10 * t11 * t17 * t19 - t10 * t22 * t12 * t24 - t10 * t22 * t17 * t28) + K[2] * (-t36 * t37 * t12 * t14 - t36 * t37 * t17 * t19 - t36 * t44 * t12 * t24 - t36 * t44 * t17 * t28);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f2
// ***********************************************************************

complex<double> Omega_function_f3_f2 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = pow(GAMMA, 0.2e1);
	complex<double> t8 = 0.1e1 / t6 / GAMMA;
	complex<double> t9 = sin(theta_p);
	complex<double> t11 = coef[1];
	complex<double> t14 = cos(theta_p);
	complex<double> t16 = coef[2];
	complex<double> t19 = t2 * t2;
	complex<double> t20 = t19 * t3;
	complex<double> t21 = sin(theta_q);
	complex<double> t23 = coef[3];
	complex<double> t29 = t19 * t4;
	complex<double> t30 = t29 * t8;
	complex<double> t32 = coef[0];
	complex<double> t35 = cos(theta_q);
	complex<double> t37 = coef[4];
	complex<double> t41 = coef[5];
	complex<double> t45 = coef[6];
	complex<double> t48 = 0.1e1 / t6;
	complex<double> t49 = t5 * t48;
	complex<double> t50 = j * ko;
	complex<double> t65 = t29 * t48 * j;
	complex<double> t70 = ko * t35;
	//
	X = K[0] * (-t5 * t8 * t9 * t11 - t5 * t8 * t14 * t16 - t20 * t8 * t21 * t23) + K[1] * (-t30 * t9 * t21 * t32 - t30 * t35 * t14 * t37 - t30 * t21 * t14 * t41 - t30 * t35 * t9 * t45 - t49 * t50 * t9 * t11 - t49 * t50 * t14 * t16 - t20 * t48 * t50 * t21 * t23) + K[2] * (-t65 * ko * t9 * t21 * t32 - t65 * t70 * t14 * t37 - t65 * ko * t21 * t14 * t41 - t65 * t70 * t9 * t45);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> Omega_function_f3_f3
// ***********************************************************************

complex<double> Omega_function_f3_f3 (double theta_p, double theta_q, double Psi, double GAMMA, complex<double> coef[], const double ko, complex<double> K[])
{                                 
	complex<double> X;
	//
	complex<double> j   = Iunit;
	//
	complex<double> t2 = sin(Psi);
	complex<double> t3 = cos(Psi);
	complex<double> t4 = t3 * t3;
	complex<double> t5 = t2 * t4;
	complex<double> t6 = pow(GAMMA, 0.2e1);
	complex<double> t8 = 0.1e1 / t6 / GAMMA;
	complex<double> t9 = sin(theta_p);
	complex<double> t11 = coef[0];
	complex<double> t14 = t2 * t2;
	complex<double> t15 = t14 * t3;
	complex<double> t16 = cos(theta_q);
	complex<double> t18 = coef[1];
	complex<double> t21 = sin(theta_q);
	complex<double> t23 = coef[2];
	complex<double> t26 = cos(theta_p);
	complex<double> t28 = coef[3];
	complex<double> t34 = t14 * t4;
	complex<double> t35 = t34 * t8;
	complex<double> t37 = coef[4];
	complex<double> t41 = coef[5];
	complex<double> t45 = coef[6];
	complex<double> t49 = coef[7];
	complex<double> t52 = 0.1e1 / t6;
	complex<double> t53 = t5 * t52;
	complex<double> t54 = j * ko;
	complex<double> t58 = t15 * t52;
	complex<double> t72 = t34 * t52 * j;
	complex<double> t73 = ko * t16;
	complex<double> t80 = ko * t21;
	//
	X = K[0] * (-t5 * t8 * t9 * t11 - t15 * t8 * t16 * t18 - t15 * t8 * t21 * t23 - t5 * t8 * t26 * t28) + K[1] * (-t35 * t16 * t9 * t37 - t35 * t16 * t26 * t41 - t35 * t21 * t9 * t45 - t35 * t21 * t26 * t49 - t53 * t54 * t9 * t11 - t58 * t54 * t16 * t18 - t58 * t54 * t21 * t23 - t53 * t54 * t26 * t28) + K[2] * (-t72 * t73 * t9 * t37 - t72 * t73 * t26 * t41 - t72 * t80 * t9 * t45 - t72 * t80 * t26 * t49);
	// Final Output
	return X;
}


// ***********************************************************************
//			IMPLEMENTATION OF void X_function_pre
// ***********************************************************************

void K_functions (double GAMMA, double L, complex<double> K[], const double ko)
{
	complex<double> a  = Iunit * ko * GAMMA;
	//
	K[0] = (1.0 / pow(a,2.0) ) * (1.0 - exp(-a * L) - a * L * exp(-a * L) );
	K[1] = (1.0 / pow(a,3.0) ) * (2.0 - 2.0 * exp(-a * L) - 2.0 * a * L * exp(-a * L) - pow(a,2.0) * pow(L,2.0) * exp(-a * L) );
	K[2] = -(pow(L,3.0) * exp(-a * L) ) / a + (3.0 / a) * K[1];
	K[3] = -(pow(L,4.0) * exp(-a * L) ) / a + (4.0 / a) * K[2];
}


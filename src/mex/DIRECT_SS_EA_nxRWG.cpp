/**************************************************************************************************************************
           
		           DIRECT_SS_EA_nxRWG.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the edge adjacent 4-D
strongly singular integrals over planar triangular elements.

  Licensing: This code is distributed under the GNU LGPL license. 

  Modified:  21 September 2011

  Author:    Athanasios Polimeridis

  References

  A. G. Polimeridis and T. V. Yioultsis, �On the direct evaluation of weakly singular
  integrals in Galerkin mixed potential integral equation formulations,� IEEE Trans.
  Antennas Propag., vol. 56, no. 9, pp. 3011-3019, Sep. 2008.

  A. G. Polimeridis and J. R. Mosig, �Complete semi-analytical treatment of weakly
  singular integrals on planar triangles via the direct evaluation method,� Int. J.
  Numerical Methods Eng., vol. 83, pp. 1625-1650, 2010.

  A. G. Polimeridis, J. M. Tamayo, J. M. Rius and J. R. Mosig, �Fast and accurate
  computation of hyper-singular integrals in Galerkin surface integral equation
  formulations via the direct evaluation method,� IEEE Trans.
  Antennas Propag., vol. 59, no. 6, pp. 2329-2340, Jun. 2011.

  A. G. Polimeridis and J. R. Mosig, �On the direct evaluation of surface integral
  equation impedance matrix elements involving point singularities,� IEEE Antennas
  Wireless Propag. Lett., vol. 10, pp. 599-602, 2011.

  INPUT DATA
  r1,r2,r3, r4 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
  N_theta = order of the Gauss-Legendre cubature for the 1-D smooth integral over theta
  N_psi   = order of the Gauss-Legendre cubature for the 1-D smooth integral over Psi
  ko = wavenumber

  OUTPUT DATA
  I_DE[9]   = 4-D strongly singular integrals
              I_DE[0]  = I_g1_f1
              I_DE[1]  = I_g1_f2
              I_DE[2]  = I_g1_f3
              I_DE[3]  = I_g2_f1
              I_DE[4]  = I_g2_f2
              I_DE[5]  = I_g2_f3
              I_DE[6]  = I_g3_f1
              I_DE[7]  = I_g3_f2
              I_DE[8]  = I_g3_f3

**************************************************************************************************************************/
#include "DIRECT_SS_EA_nxRWG.h"
#include <math.h>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void DIRECT_SS_EA_nxRWG (const double r1[],const double r2[],const double r3[],const double r4[], const double ko, const int N_theta,  const int N_psi, const double z_theta[], const double w_theta[], const double z_psi[], const double w_psi[], complex<double> I_DE[] )
{
	// ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	double rp1[3], rp2[3], rp3[3], rp21[3], rp31[3], Ap_cross[3];
	double Ap, Jp;
	double rq1[3], rq2[3], rq3[3], rq21[3], rq31[3], Aq_cross[3];
	double Aq, Jq;
	double alpha[3], beta[3], gamma[3], r21[3], r31[3], r32[3], r41[3], r42[3];
	//
	complex<double> Ising_g1_f1, Ising_g1_f2, Ising_g1_f3, Ising_g2_f1, Ising_g2_f2, Ising_g2_f3, Ising_g3_f1, Ising_g3_f2, Ising_g3_f3;
	//
	double THETA;
	double J_theta;
	double tPsiA;
	double tPsiB;
	double PsiA;
	double PsiB;
	double psi_A;
	double psi_B;
	//
	complex<double> N[12], Nm[12];
	//
	double PSI;
	double J_psi;
	double B;
	double Bm;
	//
	double b0;
	double b1;
	double b2;
	double B1;
	double B2;
	//
	double WPSI, WTHETA;
	//
	complex<double> I_theta_g1_f1, I_theta_g1_f2, I_theta_g1_f3, I_theta_g2_f1, I_theta_g2_f2, I_theta_g2_f3, I_theta_g3_f1, I_theta_g3_f2, I_theta_g3_f3;
	complex<double> I_psi_g1_f1, I_psi_g1_f2, I_psi_g1_f3, I_psi_g2_f1, I_psi_g2_f2, I_psi_g2_f3, I_psi_g3_f1, I_psi_g3_f2, I_psi_g3_f3;
	//
	complex<double> X_g1_f1, X_g1_f2, X_g1_f3, X_g2_f1, X_g2_f2, X_g2_f3, X_g3_f1, X_g3_f2, X_g3_f3;
	//
	complex<double> I_g1_f1[6], I_g1_f2[6], I_g1_f3[6], I_g2_f1[6], I_g2_f2[6], I_g2_f3[6], I_g3_f1[6], I_g3_f2[6], I_g3_f3[6]; 
	//
	complex<double> Ig1_f1, Ig1_f2, Ig1_f3, Ig2_f1, Ig2_f2, Ig2_f3, Ig3_f1, Ig3_f2, Ig3_f3; 

	// 2. Coefficients' parameters

	complex<double> coef_g1_f1[20], coefm_g1_f1[20];
	complex<double> coef_g1_f2[27], coefm_g1_f2[27];
	complex<double> coef_g1_f3[48], coefm_g1_f3[48];
	complex<double> coef_g2_f1[27], coefm_g2_f1[27];
	complex<double> coef_g2_f2[20], coefm_g2_f2[20];
	complex<double> coef_g2_f3[48], coefm_g2_f3[48];
	complex<double> coef_g3_f1[27], coefm_g3_f1[27];
	complex<double> coef_g3_f2[27], coefm_g3_f2[27];
	complex<double> coef_g3_f3[48], coefm_g3_f3[48];


	// 3.

	double theta_A;
	double theta_B;

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
			rq2[i] =  r2[i];
			rq3[i] =  r4[i];
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
			alpha[i] = (r2[i] - r1[i]) / 2;
			beta[i]  = (2 * r3[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			gamma[i] = -(2 * r4[i] - r1[i] - r2[i]) / (2.0 * sqrt(3.0));
			r21[i]   = r2[i] - r1[i];
			r31[i]   = r3[i] - r1[i];
			r32[i]   = r3[i] - r2[i];
			r41[i]   = r4[i] - r1[i];
	        r42[i]   = r4[i] - r2[i];
	}
	// Get the coefficients
	coefficients_g1_f1 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f1,  coefm_g1_f1 );
	coefficients_g1_f2 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f2,  coefm_g1_f2 );
	coefficients_g1_f3 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g1_f3,  coefm_g1_f3 );
	coefficients_g2_f1 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f1,  coefm_g2_f1 );
	coefficients_g2_f2 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f2,  coefm_g2_f2 );
	coefficients_g2_f3 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g2_f3,  coefm_g2_f3 );
	coefficients_g3_f1 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f1,  coefm_g3_f1 );
	coefficients_g3_f2 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f2,  coefm_g3_f2 );
	coefficients_g3_f3 ( r1,  r2,  r3,  r4,  ko, Ap,  coef_g3_f3,  coefm_g3_f3 );
	// Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_g1_f1[im] = (0.0 ,0.0);
		 I_g1_f2[im] = (0.0 ,0.0);
		 I_g1_f3[im] = (0.0 ,0.0);
		 I_g2_f1[im] = (0.0 ,0.0);
		 I_g2_f2[im] = (0.0 ,0.0);
		 I_g2_f3[im] = (0.0 ,0.0);
		 I_g3_f1[im] = (0.0 ,0.0);
		 I_g3_f2[im] = (0.0 ,0.0);
		 I_g3_f3[im] = (0.0 ,0.0);
	 }
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_g1_f1 = (0.0 ,0.0);
		 I_theta_g1_f2 = (0.0 ,0.0);
		 I_theta_g1_f3 = (0.0 ,0.0);
		 I_theta_g2_f1 = (0.0 ,0.0);
		 I_theta_g2_f2 = (0.0 ,0.0);
		 I_theta_g2_f3 = (0.0 ,0.0);
		 I_theta_g3_f1 = (0.0 ,0.0);
		 I_theta_g3_f2 = (0.0 ,0.0);
		 I_theta_g3_f3 = (0.0 ,0.0);
		 //
		 THETA_limits ( m, &theta_A, &theta_B );
		 //
		 for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 {
			 THETA   = ( (theta_B - theta_A) * 0.5) * z_theta[n_theta] + (theta_B + theta_A) * 0.5;
			 J_theta = (theta_B - theta_A) * 0.5;
			 tPsiA   = sin(THETA) - sqrt(3.0) * cos(THETA);
			 tPsiB   = sin(THETA) + sqrt(3.0) * cos(THETA);
			 //
			 PsiA = atan(tPsiA);
			 PsiB = atan(tPsiB);
			 // Evaluate b0, b1, b2
			 b0 = vector_dot(beta,beta);
			 b1 = 2 * ( vector_dot(alpha,beta) * cos(THETA) + vector_dot(beta,gamma) * sin(THETA));
			 b2 = vector_dot(alpha,alpha)  * pow(cos(THETA),2) + 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     // Evaluate B1, B2 (B_plus, B_minus)
			 B1 = 2 * ( vector_dot(beta,gamma) * sin(THETA) - vector_dot(alpha,beta) * cos(THETA));
			 B2 = vector_dot(alpha,alpha) * pow(cos(THETA),2) - 2 * vector_dot(alpha,gamma) * cos(THETA) * sin(THETA) + vector_dot(gamma,gamma) * pow(sin(THETA),2) ;
		     //
			 PSI_limits ( m, PsiA, PsiB, &psi_A, &psi_B );	
			 //
			 I_psi_g1_f1 = (0.0 ,0.0);
			 I_psi_g1_f2 = (0.0 ,0.0);
			 I_psi_g1_f3 = (0.0 ,0.0);
			 I_psi_g2_f1 = (0.0 ,0.0);
			 I_psi_g2_f2 = (0.0 ,0.0);
			 I_psi_g2_f3 = (0.0 ,0.0);
			 I_psi_g3_f1 = (0.0 ,0.0);
			 I_psi_g3_f2 = (0.0 ,0.0);
			 I_psi_g3_f3 = (0.0 ,0.0);
			 //
			 for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 {
				 PSI   = ( (psi_B - psi_A) * 0.5) * z_psi[n_psi] + (psi_B + psi_A) * 0.5;
				 J_psi = (psi_B - psi_A) * 0.5;
				 //				 
				 B     = sqrt(b0 * pow(sin(PSI),2) + b1 * cos(PSI) * sin(PSI) + b2 * pow(cos(PSI),2) );
				 Bm    = sqrt(b0 * pow(sin(PSI),2) + B1 * cos(PSI) * sin(PSI) + B2 * pow(cos(PSI),2) );
				 // Pre-processing: Evaluate N[12], Nm[12]
				 X_function_pre ( THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, N, Nm, ko);
				 // Evaluate the integrands
				 X_g1_f1 = X_function_g1_f1 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g1_f1, coefm_g1_f1, N, Nm);
				 X_g1_f2 = X_function_g1_f2 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g1_f2, coefm_g1_f2, N, Nm);
				 X_g1_f3 = X_function_g1_f3 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g1_f3, coefm_g1_f3, N, Nm);
				 X_g2_f1 = X_function_g2_f1 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g2_f1, coefm_g2_f1, N, Nm);
				 X_g2_f2 = X_function_g2_f2 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g2_f2, coefm_g2_f2, N, Nm);
				 X_g2_f3 = X_function_g2_f3 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g2_f3, coefm_g2_f3, N, Nm);
				 X_g3_f1 = X_function_g3_f1 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g3_f1, coefm_g3_f1, N, Nm);
				 X_g3_f2 = X_function_g3_f2 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g3_f2, coefm_g3_f2, N, Nm);
				 X_g3_f3 = X_function_g3_f3 (THETA, PSI, tPsiA, tPsiB, PsiA, PsiB, B, Bm, coef_g3_f3, coefm_g3_f3, N, Nm);
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_g1_f1 +=  WPSI * X_g1_f1;
				 I_psi_g1_f2 +=  WPSI * X_g1_f2;
				 I_psi_g1_f3 +=  WPSI * X_g1_f3;
				 I_psi_g2_f1 +=  WPSI * X_g2_f1;
				 I_psi_g2_f2 +=  WPSI * X_g2_f2;
				 I_psi_g2_f3 +=  WPSI * X_g2_f3;
				 I_psi_g3_f1 +=  WPSI * X_g3_f1;
				 I_psi_g3_f2 +=  WPSI * X_g3_f2;
				 I_psi_g3_f3 +=  WPSI * X_g3_f3;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_g1_f1 *=  J_psi;
			 I_psi_g1_f2 *=  J_psi;
			 I_psi_g1_f3 *=  J_psi;
			 I_psi_g2_f1 *=  J_psi;
			 I_psi_g2_f2 *=  J_psi;
			 I_psi_g2_f3 *=  J_psi;
			 I_psi_g3_f1 *=  J_psi;
			 I_psi_g3_f2 *=  J_psi;
			 I_psi_g3_f3 *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_g1_f1 +=  WTHETA * I_psi_g1_f1;
			 I_theta_g1_f2 +=  WTHETA * I_psi_g1_f2;
			 I_theta_g1_f3 +=  WTHETA * I_psi_g1_f3;
			 I_theta_g2_f1 +=  WTHETA * I_psi_g2_f1;
			 I_theta_g2_f2 +=  WTHETA * I_psi_g2_f2;
			 I_theta_g2_f3 +=  WTHETA * I_psi_g2_f3;
			 I_theta_g3_f1 +=  WTHETA * I_psi_g3_f1;
			 I_theta_g3_f2 +=  WTHETA * I_psi_g3_f2;
			 I_theta_g3_f3 +=  WTHETA * I_psi_g3_f3;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_g1_f1[m-1] = J_theta * I_theta_g1_f1;
		 I_g1_f2[m-1] = J_theta * I_theta_g1_f2;
		 I_g1_f3[m-1] = J_theta * I_theta_g1_f3;
		 I_g2_f1[m-1] = J_theta * I_theta_g2_f1;
		 I_g2_f2[m-1] = J_theta * I_theta_g2_f2;
		 I_g2_f3[m-1] = J_theta * I_theta_g2_f3;
		 I_g3_f1[m-1] = J_theta * I_theta_g3_f1;
		 I_g3_f2[m-1] = J_theta * I_theta_g3_f2;
		 I_g3_f3[m-1] = J_theta * I_theta_g3_f3;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 Ig1_f1 = I_g1_f1[0] + I_g1_f1[1] + I_g1_f1[2] + I_g1_f1[3] + I_g1_f1[4] + I_g1_f1[5];
	 Ig1_f2 = I_g1_f2[0] + I_g1_f2[1] + I_g1_f2[2] + I_g1_f2[3] + I_g1_f2[4] + I_g1_f2[5];
	 Ig1_f3 = I_g1_f3[0] + I_g1_f3[1] + I_g1_f3[2] + I_g1_f3[3] + I_g1_f3[4] + I_g1_f3[5];
	 Ig2_f1 = I_g2_f1[0] + I_g2_f1[1] + I_g2_f1[2] + I_g2_f1[3] + I_g2_f1[4] + I_g2_f1[5];
	 Ig2_f2 = I_g2_f2[0] + I_g2_f2[1] + I_g2_f2[2] + I_g2_f2[3] + I_g2_f2[4] + I_g2_f2[5];
	 Ig2_f3 = I_g2_f3[0] + I_g2_f3[1] + I_g2_f3[2] + I_g2_f3[3] + I_g2_f3[4] + I_g2_f3[5];
	 Ig3_f1 = I_g3_f1[0] + I_g3_f1[1] + I_g3_f1[2] + I_g3_f1[3] + I_g3_f1[4] + I_g3_f1[5];
	 Ig3_f2 = I_g3_f2[0] + I_g3_f2[1] + I_g3_f2[2] + I_g3_f2[3] + I_g3_f2[4] + I_g3_f2[5];
	 Ig3_f3 = I_g3_f3[0] + I_g3_f3[1] + I_g3_f3[2] + I_g3_f3[3] + I_g3_f3[4] + I_g3_f3[5];
	 //
	 double PRECOMP = Jp * Jq / (4 * Ap * Aq);
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * Ig1_f1;
	 I_DE[1] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r42,r42)) * Ig1_f2;
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r21,r21)) * Ig1_f3;
	 I_DE[3] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r41,r41)) * Ig2_f1;
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r42,r42)) * Ig2_f2;
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r21,r21)) * Ig2_f3;
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * Ig3_f1;
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r42,r42)) * Ig3_f2;
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r21,r21)) * Ig3_f3;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g1_f1
// ***********************************************************************

void coefficients_g1_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[4], cm[4];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	complex<double> j   = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = pow(bb[0], 0.2e1);
	complex<double> t3 = t1 * t2;
	complex<double> t5 = cc[0] * bb[2] * dd[1];
	complex<double> t8 = cc[2] * dd[0];
	complex<double> t9 = t8 * bb[1];
	complex<double> t12 = pow(cc[1], 0.2e1);
	complex<double> t13 = t1 * t12;
	complex<double> t14 = bb[0] * dd[2];
	complex<double> t15 = t14 * bb[1];
	complex<double> t19 = cc[0] * bb[1] * dd[2];
	complex<double> t22 = pow(cc[2], 0.2e1);
	complex<double> t23 = t1 * t22;
	complex<double> t24 = bb[1] * dd[0];
	complex<double> t25 = t24 * bb[2];
	complex<double> t28 = pow(cc[0], 0.2e1);
	complex<double> t29 = t1 * t28;
	complex<double> t32 = pow(bb[2], 0.2e1);
	complex<double> t33 = t1 * t32;
	complex<double> t36 = t1 * cc[0];
	complex<double> t41 = t1 * cc[1];
	complex<double> t42 = pow(bb[1], 0.2e1);
	complex<double> t49 = t1 * t42;
	complex<double> t50 = cc[1] * bb[0];
	complex<double> t51 = t50 * dd[2];
	complex<double> t54 = bb[0] * dd[1];
	complex<double> t55 = t54 * bb[2];
	complex<double> t58 = t1 * cc[2];
	complex<double> t67 = cc[1] * dd[0];
	complex<double> t68 = t67 * bb[2];
	complex<double> t72 = cc[2] * bb[0] * dd[1];
	complex<double> t79 = t3 * t5 / 0.12e2 + t3 * t9 / 0.12e2 - t13 * t15 / 0.6e1 - t3 * t19 / 0.12e2 - t23 * t25 / 0.6e1 + t29 * t15 / 0.6e1 + t33 * t9 / 0.12e2 + t36 * cc[2] * t2 * dd[1] / 0.6e1 - t41 * cc[2] * t42 * dd[0] / 0.6e1 + t13 * t25 / 0.6e1 + t49 * t51 / 0.12e2 + t23 * t55 / 0.6e1 - t58 * cc[0] * t32 * dd[1] / 0.6e1 + t41 * cc[0] * t42 * dd[2] / 0.6e1 - t3 * t68 / 0.12e2 - t33 * t72 / 0.12e2 - t33 * t19 / 0.12e2 + t49 * t5 / 0.12e2;
	complex<double> t97 = t1 * t2 * bb[0];
	complex<double> t98 = cc[2] * dd[1];
	complex<double> t102 = t1 * t42 * bb[1];
	complex<double> t103 = cc[0] * dd[2];
	complex<double> t107 = t1 * t32 * bb[2];
	complex<double> t108 = cc[0] * dd[1];
	complex<double> t111 = cc[1] * dd[2];
	complex<double> t118 = t41 * cc[2];
	complex<double> t119 = bb[1] * dd[1];
	complex<double> t123 = t36 * cc[2];
	complex<double> t124 = bb[0] * dd[0];
	complex<double> t129 = t1 * bb[2] * cc[2];
	complex<double> t132 = bb[1] * bb[2];
	complex<double> t136 = t41 * cc[0];
	complex<double> t140 = bb[2] * bb[0];
	complex<double> t144 = -t36 * cc[1] * t2 * dd[2] / 0.6e1 + t58 * cc[1] * t32 * dd[0] / 0.6e1 - t29 * t55 / 0.6e1 + t33 * t51 / 0.12e2 - t49 * t72 / 0.12e2 - t49 * t68 / 0.12e2 - t97 * t98 / 0.12e2 - t102 * t103 / 0.12e2 + t107 * t108 / 0.12e2 + t97 * t111 / 0.12e2 + t102 * t8 / 0.12e2 - t107 * t67 / 0.12e2 + t118 * t119 * bb[0] / 0.6e1 - t123 * t124 * bb[1] / 0.6e1 - t129 * t51 / 0.6e1 + t123 * t132 * dd[2] / 0.6e1 - t136 * t132 * dd[1] / 0.6e1 + t136 * t140 * dd[0] / 0.6e1;
	complex<double> t146 = sqrt(0.3e1);
	complex<double> t147 = dd[2] * t146;
	complex<double> t171 = t102 * t8 * t146;
	complex<double> t174 = t107 * t67 * t146;
	complex<double> t177 = t97 * t98 * t146;
	complex<double> t180 = t97 * t111 * t146;
	complex<double> t183 = t107 * t108 * t146;
	complex<double> t186 = t102 * t103 * t146;
	complex<double> t188 = t146 * t1;
	complex<double> t189 = t22 * cc[2];
	complex<double> t194 = -t123 * t132 * t147 / 0.18e2 + t129 * t50 * t147 / 0.18e2 + t136 * t132 * dd[1] * t146 / 0.18e2 - t136 * t140 * dd[0] * t146 / 0.18e2 + t123 * t124 * bb[1] * t146 / 0.18e2 - t118 * t119 * bb[0] * t146 / 0.18e2 - t171 / 0.72e2 + t174 / 0.72e2 + t177 / 0.72e2 - t180 / 0.72e2 - t183 / 0.72e2 + t186 / 0.72e2 + t188 * t189 * bb[0] * dd[1] / 0.18e2;
	complex<double> t195 = t12 * cc[1];
	complex<double> t200 = t28 * cc[0];
	complex<double> t217 = t188 * t28;
	complex<double> t220 = t188 * t22;
	complex<double> t223 = t188 * t12;
	complex<double> t232 = t33 * cc[2];
	complex<double> t233 = t24 * t146;
	complex<double> t234 = t232 * t233;
	complex<double> t236 = t49 * cc[1];
	complex<double> t237 = t14 * t146;
	complex<double> t238 = t236 * t237;
	complex<double> t242 = dd[1] * bb[2] * t146;
	complex<double> t243 = t49 * cc[0] * t242;
	complex<double> t245 = t188 * t195 * bb[2] * dd[0] / 0.18e2 - t188 * t200 * bb[2] * dd[1] / 0.18e2 - t188 * t189 * bb[1] * dd[0] / 0.18e2 - t188 * t195 * bb[0] * dd[2] / 0.18e2 + t188 * t200 * bb[1] * dd[2] / 0.18e2 + t217 * t72 / 0.18e2 + t220 * t19 / 0.18e2 - t223 * t9 / 0.18e2 + t220 * t68 / 0.18e2 - t223 * t5 / 0.18e2 - t217 * t9 / 0.18e2 - t234 / 0.72e2 - t238 / 0.72e2 - t243 / 0.72e2;
	complex<double> t249 = dd[2] * bb[1] * t146;
	complex<double> t250 = t33 * cc[0] * t249;
	complex<double> t252 = t54 * t146;
	complex<double> t253 = t232 * t252;
	complex<double> t255 = t3 * cc[0];
	complex<double> t256 = t255 * t242;
	complex<double> t260 = dd[0] * bb[2] * t146;
	complex<double> t261 = t3 * cc[1] * t260;
	complex<double> t264 = t3 * cc[2] * t233;
	complex<double> t266 = t255 * t249;
	complex<double> t269 = t49 * cc[2] * t252;
	complex<double> t272 = t33 * cc[1] * t237;
	complex<double> t291 = t250 / 0.72e2 + t253 / 0.72e2 - t256 / 0.72e2 + t261 / 0.72e2 - t264 / 0.72e2 + t266 / 0.72e2 + t269 / 0.72e2 - t272 / 0.72e2 - t118 * t32 * dd[0] * t146 / 0.18e2 - t23 * bb[0] * t242 / 0.18e2 - t29 * bb[1] * t237 / 0.18e2 + t29 * bb[2] * t252 / 0.18e2 - t123 * t2 * dd[1] * t146 / 0.18e2;
	complex<double> t329 = t236 * t260;
	complex<double> t331 = t223 * t72 / 0.18e2 - t217 * t51 / 0.18e2 + t217 * t68 / 0.18e2 - t220 * t5 / 0.18e2 + t223 * t19 / 0.18e2 - t220 * t51 / 0.18e2 + t136 * t2 * dd[2] * t146 / 0.18e2 + t13 * bb[0] * t249 / 0.18e2 - t13 * bb[2] * t233 / 0.18e2 + t123 * t32 * dd[1] * t146 / 0.18e2 - t136 * t42 * dd[2] * t146 / 0.18e2 + t118 * t42 * dd[0] * t146 / 0.18e2 + t23 * bb[1] * t260 / 0.18e2 + t329 / 0.72e2;
	complex<double> t334 = -t261 - t266 + t243 - t186 + t234 + t183 + t171 + t180 - t253 - t250 + t238 + t272 - t177 - t329 + t256 + t264 - t269 - t174;
	//
	c[0] = t79 + t144;
	c[1] = t194 + t245 + t291 + t331;
	c[2] = t334 / 0.24e2;
	c[3] = -t334 / 0.24e2;
	//
	t3 =  (0.1e1 /  ko /  j);
	t5 = 2.0 * t3 * c[2];
	complex<double> t6 =   ko *  ko;
	t9 = j *  j;
	t14 = 8.0 / t6 / ko / t9 / j * c[1];
	complex<double> t17 = 1.0 / t6 / t9;
	t19 = 3.0 * t17 * c[0];
	complex<double> t21 = 2.0 * t3 * c[3];
	t23 = 8.0 * t17 * c[1];
	t25 = 4.0 * t3 * c[1];
	//
	coef[0] = -t5;
	coef[1] = -t14;
	coef[2] = t14;
	coef[3] = -t19;
	coef[4] = t21;
	coef[5] = t23;
	coef[6] = -t23;
	coef[7] = t25;
	coef[8] = -t25;
	coef[9] = c[2];
	coef[10] = -c[1];
	coef[11] = c[1];
	coef[12] = 3.0 * t3 * c[0];
	coef[13] = c[0];
	coef[14] = t19;
	coef[15] = c[3];
	coef[16] = -t21;
	coef[17] = t5;
	coef[18] = t14;
	coef[19] = -t14;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[1], 0.2e1);
	t3 = t1 * t2;
	complex<double> t4 = cc[2] * bb[0];
	t5 = t4 * dd[1];
	t8 = pow(bb[0], 0.2e1);
	t9 = t1 * t8;
	complex<double> t10 = cc[2] * bb[1];
	complex<double> t11 = t10 * dd[0];
	t15 = cc[0] * bb[2] * dd[1];
	complex<double> t18 = pow(bb[2], 0.2e1);
	t19 = t1 * t18;
	t22 = t1 * cc[1];
	complex<double> t30 = cc[1] * bb[2] * dd[0];
	complex<double> t34 = cc[1] * bb[0] * dd[2];
	complex<double> t38 = cc[0] * bb[1] * dd[2];
	t41 = pow(cc[1], 0.2e1);
	t42 = t1 * t41;
	t55 = t1 * bb[2];
	complex<double> t56 = pow(cc[2], 0.2e1);
	complex<double> t62 = cc[0] * cc[1] * dd[2];
	complex<double> t65 = t1 * bb[0];
	complex<double> t66 = pow(cc[0], 0.2e1);
	complex<double> t73 = t3 * t5 / 0.12e2 - t9 * t11 / 0.12e2 - t9 * t15 / 0.12e2 + t19 * t5 / 0.12e2 + t22 * cc[2] * t2 * dd[0] / 0.6e1 - t19 * t11 / 0.12e2 + t9 * t30 / 0.12e2 - t19 * t34 / 0.12e2 + t9 * t38 / 0.12e2 + t42 * bb[0] * bb[1] * dd[2] / 0.6e1 + t3 * t30 / 0.12e2 - t3 * t34 / 0.12e2 + t19 * t38 / 0.12e2 - t3 * t15 / 0.12e2 - t55 * t56 * bb[0] * dd[1] / 0.6e1 - t3 * t62 / 0.6e1 - t65 * t66 * bb[1] * dd[2] / 0.6e1 + t9 * t62 / 0.6e1;
	complex<double> t82 = t1 * bb[1];
	complex<double> t91 = t1 * t66;
	complex<double> t92 = bb[2] * bb[0];
	t102 = t1 * t18 * bb[2];
	t103 = cc[1] * dd[0];
	t107 = t1 * t2 * bb[1];
	t108 = cc[2] * dd[0];
	complex<double> t112 = t1 * t8 * bb[0];
	complex<double> t113 = cc[1] * dd[2];
	complex<double> t116 = cc[2] * dd[1];
	t119 = cc[0] * dd[2];
	complex<double> t122 = cc[0] * dd[1];
	complex<double> t125 = t1 * cc[2];
	complex<double> t126 = t125 * cc[0];
	complex<double> t127 = bb[2] * bb[1];
	complex<double> t131 = t82 * cc[1];
	complex<double> t134 = t22 * cc[0];
	complex<double> t138 = t65 * cc[0];
	t144 = t125 * cc[1];
	complex<double> t148 = -t19 * cc[1] * cc[2] * dd[0] / 0.6e1 + t55 * t56 * bb[1] * dd[0] / 0.6e1 - t82 * t41 * bb[2] * dd[0] / 0.6e1 + t19 * cc[0] * cc[2] * dd[1] / 0.6e1 + t91 * t92 * dd[1] / 0.6e1 - t1 * cc[0] * cc[2] * t8 * dd[1] / 0.6e1 + t102 * t103 / 0.12e2 - t107 * t108 / 0.12e2 - t112 * t113 / 0.12e2 + t112 * t116 / 0.12e2 + t107 * t119 / 0.12e2 - t102 * t122 / 0.12e2 - t126 * t127 * dd[2] / 0.6e1 - t131 * t5 / 0.6e1 + t134 * t127 * dd[1] / 0.6e1 + t138 * t11 / 0.6e1 - t134 * t92 * dd[0] / 0.6e1 + t144 * t92 * dd[2] / 0.6e1;
	complex<double> t150 = t19 * cc[2];
	complex<double> t152 = sqrt(0.3e1);
	complex<double> t153 = bb[1] * dd[0] * t152;
	complex<double> t154 = t150 * t153;
	complex<double> t155 = t9 * cc[0];
	complex<double> t157 = bb[2] * dd[1] * t152;
	complex<double> t158 = t155 * t157;
	complex<double> t160 = bb[1] * dd[2] * t152;
	complex<double> t161 = t155 * t160;
	complex<double> t162 = t3 * cc[1];
	complex<double> t164 = bb[2] * dd[0] * t152;
	complex<double> t165 = t162 * t164;
	complex<double> t167 = bb[0] * dd[2] * t152;
	complex<double> t168 = t162 * t167;
	complex<double> t169 = t19 * cc[0];
	complex<double> t170 = t169 * t160;
	t171 = t3 * cc[0];
	complex<double> t172 = t171 * t157;
	complex<double> t173 = t9 * cc[1];
	t174 = t173 * t164;
	t177 = bb[0] * dd[1] * t152;
	complex<double> t178 = t3 * cc[2] * t177;
	t180 = t9 * cc[2] * t153;
	complex<double> t182 = t102 * t103 * t152;
	complex<double> t184 = t102 * t122 * t152;
	complex<double> t185 = t119 * t152;
	t186 = t107 * t185;
	complex<double> t187 = t150 * t177;
	t188 = t108 * t152;
	t189 = t107 * t188;
	complex<double> t190 = t113 * t152;
	complex<double> t191 = t112 * t190;
	complex<double> t192 = t19 * cc[1];
	complex<double> t193 = t192 * t167;
	t194 = t116 * t152;
	t195 = t112 * t194;
	complex<double> t196 = t154 + t158 - t161 - t165 + t168 - t170 + t172 - t174 - t178 + t180 - t182 + t184 - t186 - t187 + t189 + t191 + t193 - t195;
	complex<double> t203 = t152 * t1;
	complex<double> t204 = t66 * cc[0];
	complex<double> t209 = t41 * cc[1];
	complex<double> t218 = t56 * cc[2];
	complex<double> t231 = t203 * t66;
	t234 = t182 / 0.72e2 - t184 / 0.72e2 + t186 / 0.72e2 - t189 / 0.72e2 - t191 / 0.72e2 + t195 / 0.72e2 + t203 * t204 * bb[1] * dd[2] / 0.18e2 + t203 * t209 * bb[2] * dd[0] / 0.18e2 - t203 * t204 * bb[2] * dd[1] / 0.18e2 + t203 * t218 * bb[0] * dd[1] / 0.18e2 - t203 * t209 * bb[0] * dd[2] / 0.18e2 - t203 * t218 * bb[1] * dd[0] / 0.18e2 + t231 * t30 / 0.18e2;
	complex<double> t235 = t203 * t41;
	complex<double> t239 = t203 * t56;
	t253 = t235 * t5 - t235 * t11 - t231 * t11 + t239 * t38 - t231 * t34 + t235 * t38 + t239 * t30 + t231 * t5 - t239 * t15 - t239 * t34 - t235 * t15 + t169 * t194 + t91 * bb[2] * t177 - t82 * t41 * t164;
	complex<double> t258 = t55 * t56;
	complex<double> t284 = t42 * bb[0] * t160 / 0.18e2 + t258 * t153 / 0.18e2 - t192 * t188 / 0.18e2 + t173 * t185 / 0.18e2 - t65 * t66 * t160 / 0.18e2 - t171 * t190 / 0.18e2 - t126 * t8 * dd[1] * t152 / 0.18e2 + t144 * t2 * dd[0] * t152 / 0.18e2 - t258 * t177 / 0.18e2 - t154 / 0.72e2 - t158 / 0.72e2 + t161 / 0.72e2 + t165 / 0.72e2;
	complex<double> t293 = dd[0] * t152;
	complex<double> t297 = dd[1] * t152;
	complex<double> t304 = dd[2] * t152;
	complex<double> t314 = -t168 / 0.72e2 + t170 / 0.72e2 - t172 / 0.72e2 + t174 / 0.72e2 + t178 / 0.72e2 - t180 / 0.72e2 + t187 / 0.72e2 - t193 / 0.72e2 + t138 * t10 * t293 / 0.18e2 + t134 * t127 * t297 / 0.18e2 - t134 * t92 * t293 / 0.18e2 + t144 * t92 * t304 / 0.18e2 - t131 * t4 * t297 / 0.18e2 - t126 * t127 * t304 / 0.18e2;
	//
	cm[0] = t73 + t148;
	cm[1] = t196 / 0.24e2;
	cm[2] = -t196 / 0.24e2;
	cm[3] = t234 + t253 / 0.18e2 + t284 + t314;
	//
	t3 =  (0.1e1 /  ko /  j);
	t5 = 2.0 * t3 * cm[1];
	t6 =  ko * ko;
	t9 =  j *  j;
	t14 = 8.0 / t6 / ko / t9 / j * cm[3];
	t19 = 1.0 / t6 / t9;
	t21 = 3.0 * t19 * cm[0];
	t23 = 2.0 * t3 * cm[2];
	t25 = 8.0 * t19 * cm[3];
	complex<double> t27 = 4.0 * t3 * cm[3];
	//
	coefm[0] = -t5;
	coefm[1] = -t14;
	coefm[2] = t14;
	coefm[3] = 3.0 * t3 * cm[0];
	coefm[4] = cm[0];
	coefm[5] = -cm[3];
	coefm[6] = cm[3];
	coefm[7] = cm[1];
	coefm[8] = -t21;
	coefm[9] = t23;
	coefm[10] = t25;
	coefm[11] = -t25;
	coefm[12] = t27;
	coefm[13] = -t27;
	coefm[14] = cm[2];
	coefm[15] = t21;
	coefm[16] = -t23;
	coefm[17] = t5;
	coefm[18] = t14;
	coefm[19] = -t14;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g1_f2
// ***********************************************************************

void coefficients_g1_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[6], cm[6];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	complex<double> j   = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = t1 * bb[2];
	complex<double> t3 = t2 * cc[2];
	complex<double> t4 = cc[1] * bb[0];
	complex<double> t5 = t4 * dd[2];
	complex<double> t8 = t1 * bb[1];
	complex<double> t9 = t8 * cc[1];
	complex<double> t10 = cc[2] * bb[0];
	complex<double> t11 = t10 * dd[1];
	complex<double> t15 = t1 * cc[2] * cc[0];
	complex<double> t16 = bb[1] * bb[2];
	complex<double> t21 = t1 * cc[1] * cc[0];
	complex<double> t25 = t1 * bb[0];
	complex<double> t26 = t25 * cc[0];
	complex<double> t27 = cc[2] * bb[1];
	complex<double> t28 = t27 * dd[0];
	complex<double> t31 = bb[2] * bb[0];
	complex<double> t35 = pow(cc[2], 0.2e1);
	complex<double> t40 = pow(bb[1], 0.2e1);
	complex<double> t41 = t1 * t40;
	complex<double> t43 = cc[0] * cc[1] * dd[2];
	complex<double> t46 = pow(bb[2], 0.2e1);
	complex<double> t47 = t1 * t46;
	complex<double> t51 = cc[0] * bb[1] * dd[2];
	complex<double> t55 = cc[0] * bb[2] * dd[1];
	complex<double> t60 = pow(bb[0], 0.2e1);
	complex<double> t61 = t1 * t60;
	complex<double> t65 = cc[1] * bb[2] * dd[0];
	complex<double> t76 = -t3 * t5 / 0.6e1 + t9 * t11 / 0.6e1 + t15 * t16 * dd[2] / 0.6e1 - t21 * t16 * dd[1] / 0.6e1 - t26 * t28 / 0.6e1 + t21 * t31 * dd[0] / 0.6e1 + t2 * t35 * bb[0] * dd[1] / 0.6e1 + t41 * t43 / 0.6e1 - t47 * t11 / 0.12e2 - t47 * t51 / 0.12e2 + t41 * t55 / 0.12e2 + t41 * t5 / 0.12e2 - t61 * t51 / 0.12e2 - t61 * t65 / 0.12e2 + t47 * t5 / 0.12e2 + t47 * t28 / 0.12e2 - t41 * t65 / 0.12e2 - t41 * t11 / 0.12e2;
	complex<double> t81 = t40 * bb[1];
	complex<double> t82 = t1 * t81;
	complex<double> t83 = cc[0] * dd[2];
	complex<double> t86 = t60 * bb[0];
	complex<double> t87 = t1 * t86;
	complex<double> t88 = cc[1] * dd[2];
	complex<double> t91 = t46 * bb[2];
	complex<double> t92 = t1 * t91;
	complex<double> t93 = cc[1] * dd[0];
	complex<double> t96 = cc[2] * dd[0];
	complex<double> t99 = cc[2] * dd[1];
	complex<double> t102 = cc[0] * dd[1];
	complex<double> t106 = cc[0] * cc[2] * dd[1];
	complex<double> t109 = pow(cc[1], 0.2e1);
	complex<double> t114 = pow(cc[0], 0.2e1);
	complex<double> t122 = cc[1] * cc[2] * dd[0];
	complex<double> t141 = t61 * t28 / 0.12e2 + t61 * t55 / 0.12e2 - t82 * t83 / 0.12e2 + t87 * t88 / 0.12e2 - t92 * t93 / 0.12e2 + t82 * t96 / 0.12e2 - t87 * t99 / 0.12e2 + t92 * t102 / 0.12e2 - t47 * t106 / 0.6e1 - t8 * t109 * bb[0] * dd[2] / 0.6e1 + t25 * t114 * bb[1] * dd[2] / 0.6e1 - t61 * t43 / 0.6e1 + t47 * t122 / 0.6e1 - t2 * t35 * bb[1] * dd[0] / 0.6e1 + t8 * t109 * bb[2] * dd[0] / 0.6e1 - t41 * t122 / 0.6e1 + t61 * t106 / 0.6e1 - t25 * t114 * bb[2] * dd[1] / 0.6e1;
	complex<double> t142 = t76 + t141;
	complex<double> t143 = sqrt(0.3e1);
	complex<double> t144 = t143 * t1;
	complex<double> t145 = t144 * t46;
	complex<double> t146 = t145 * t51;
	complex<double> t147 = t144 * t40;
	complex<double> t148 = t147 * t55;
	complex<double> t149 = t147 * t5;
	complex<double> t152 = t144 * t86 * cc[1] * dd[2];
	complex<double> t153 = t147 * t11;
	complex<double> t156 = t144 * t81 * cc[2] * dd[0];
	complex<double> t159 = t144 * t86 * cc[2] * dd[1];
	complex<double> t160 = t145 * t11;
	complex<double> t163 = t144 * t91 * cc[1] * dd[0];
	complex<double> t164 = t144 * t60;
	complex<double> t165 = t164 * t55;
	complex<double> t166 = t147 * t65;
	complex<double> t167 = t145 * t28;
	complex<double> t170 = t144 * t81 * cc[0] * dd[2];
	complex<double> t171 = t164 * t51;
	complex<double> t172 = t164 * t28;
	complex<double> t175 = t144 * t91 * cc[0] * dd[1];
	complex<double> t176 = t145 * t5;
	complex<double> t177 = t164 * t65;
	complex<double> t178 = t146 - t148 - t149 - t152 + t153 - t156 + t159 + t160 + t163 - t165 + t166 - t167 + t170 + t171 - t172 - t175 - t176 + t177;
	complex<double> t179 = dd[2] * t143;
	complex<double> t183 = dd[1] * t143;
	complex<double> t193 = dd[0] * t143;
	complex<double> t207 = t3 * t4 * t179 / 0.18e2 - t9 * t10 * t183 / 0.18e2 - t15 * t16 * t179 / 0.18e2 + t21 * t16 * t183 / 0.18e2 + t26 * t27 * t193 / 0.18e2 - t21 * t31 * t193 / 0.18e2 + t146 / 0.72e2 - t148 / 0.72e2 - t149 / 0.72e2 + t153 / 0.72e2 + t160 / 0.72e2 - t165 / 0.72e2 + t166 / 0.72e2;
	complex<double> t217 = t25 * t114;
	complex<double> t226 = t2 * t35;
	complex<double> t235 = t8 * t109;
	complex<double> t252 = -t167 / 0.72e2 + t171 / 0.72e2 - t172 / 0.72e2 - t176 / 0.72e2 + t177 / 0.72e2 - t61 * cc[2] * t102 * t143 / 0.18e2 + t217 * bb[2] * dd[1] * t143 / 0.18e2 - t41 * cc[0] * t88 * t143 / 0.18e2 - t226 * bb[0] * dd[1] * t143 / 0.18e2 + t47 * cc[0] * t99 * t143 / 0.18e2 + t235 * bb[0] * dd[2] * t143 / 0.18e2 - t217 * bb[1] * dd[2] * t143 / 0.18e2 + t61 * cc[1] * t83 * t143 / 0.18e2 - t47 * cc[1] * t96 * t143 / 0.18e2;
	complex<double> t263 = t144 * t114;
	complex<double> t265 = t144 * t35;
	complex<double> t267 = t144 * t109;
	complex<double> t276 = t226 * bb[1] * dd[0] * t143 - t235 * bb[2] * dd[0] * t143 + t41 * cc[2] * t93 * t143 - t263 * t5 - t265 * t5 + t267 * t11 - t263 * t28 - t267 * t28 + t263 * t65 - t267 * t55 + t265 * t65 + t265 * t51 - t265 * t55;
	complex<double> t287 = t35 * cc[2];
	complex<double> t292 = t109 * cc[1];
	complex<double> t297 = t114 * cc[0];
	complex<double> t314 = t267 * t51 / 0.18e2 + t263 * t11 / 0.18e2 - t152 / 0.72e2 - t156 / 0.72e2 + t159 / 0.72e2 + t163 / 0.72e2 + t170 / 0.72e2 - t175 / 0.72e2 + t144 * t287 * bb[0] * dd[1] / 0.18e2 - t144 * t292 * bb[0] * dd[2] / 0.18e2 + t144 * t297 * bb[1] * dd[2] / 0.18e2 - t144 * t287 * bb[1] * dd[0] / 0.18e2 + t144 * t292 * bb[2] * dd[0] / 0.18e2 - t144 * t297 * bb[2] * dd[1] / 0.18e2;
	//
	c[0] = t142;
	c[1] = t178 / 0.12e2;
	c[2] = t178 / 0.24e2;
	c[3] = t178 / 0.24e2;
	c[4] = t142;
	c[5] = t207 + t252 + t276 / 0.18e2 + t314;
	//
	t1 =  ko *  ko;
	t4 =  j *  j;
	t9 = 8.0 / t1 / ko / t4 / j * c[5];
	complex<double> t12 = 1.0 / ko / j;
	complex<double> t14 = 2.0 * t12 * c[2];
	t16 = 2.0 * t12 * c[1];
	complex<double> t19 = 1.0 / t1 / t4;
	t21 = 3.0 * t19 * c[0];
	complex<double> t23 = 2.0 * t12 * c[3];
	t25 = 3.0 * t19 * c[4];
	complex<double> t29 = 4.0 * t12 * c[5];
	t31 = 8.0 * t19 * c[5];
	//
	coef[0] = -t9;
	coef[1] = -t14;
	coef[2] = -t16;
	coef[3] = -t21;
	coef[4] = t9;
	coef[5] = -t9;
	coef[6] = t23;
	coef[7] = t9;
	coef[8] = -t23;
	coef[9] = c[2];
	coef[10] = c[1];
	coef[11] = c[0];
	coef[12] = t25;
	coef[13] = 3.0 * t12 * c[0];
	coef[14] = -c[5];
	coef[15] = c[3];
	coef[16] = c[5];
	coef[17] = -t25;
	coef[18] = t29;
	coef[19] = t14;
	coef[20] = -t31;
	coef[21] = -t29;
	coef[22] = t16;
	coef[23] = t21;
	coef[24] = t31;
	coef[25] = c[4];
	coef[26] = 3.0 * t12 * c[4];
	//
	t1 = sqrt(0.3e1);
	t2 = 0.1e1 / AreaT;
	t3 = t1 * t2;
	t4 = pow(bb[2], 0.2e1);
	t5 = t3 * t4;
	complex<double> t6 = cc[2] * bb[1];
	complex<double> t7 = t6 * dd[0];
	t8 = t5 * t7;
	t9 = pow(bb[1], 0.2e1);
	t10 = t3 * t9;
	t12 = cc[1] * bb[2] * dd[0];
	complex<double> t13 = t10 * t12;
	t14 = cc[1] * bb[0];
	t15 = t14 * dd[2];
	t16 = t10 * t15;
	complex<double> t17 = pow(bb[0], 0.2e1);
	complex<double> t18 = t3 * t17;
	complex<double> t20 = cc[0] * bb[2] * dd[1];
	t21 = t18 * t20;
	complex<double> t22 = t9 * bb[1];
	t25 = t3 * t22 * cc[2] * dd[0];
	t26 = t18 * t12;
	t27 = t5 * t15;
	t28 = t17 * bb[0];
	t31 = t3 * t28 * cc[1] * dd[2];
	complex<double> t32 = cc[2] * bb[0];
	complex<double> t33 = t32 * dd[1];
	complex<double> t34 = t5 * t33;
	complex<double> t36 = cc[0] * bb[1] * dd[2];
	complex<double> t37 = t18 * t36;
	t40 = t3 * t22 * cc[0] * dd[2];
	t41 = t10 * t20;
	complex<double> t44 = t3 * t28 * cc[2] * dd[1];
	complex<double> t45 = t4 * bb[2];
	complex<double> t48 = t3 * t45 * cc[0] * dd[1];
	t51 = t3 * t45 * cc[1] * dd[0];
	complex<double> t52 = t10 * t33;
	complex<double> t53 = t18 * t7;
	complex<double> t54 = t5 * t36;
	t55 = -t8 + t13 - t16 - t21 - t25 + t26 - t27 - t31 + t34 + t37 + t40 - t41 + t44 - t48 + t51 + t52 - t53 + t54;
	complex<double> t56 = t2 * bb[2];
	complex<double> t57 = t56 * cc[2];
	t60 = t2 * bb[1];
	t61 = t60 * cc[1];
	t65 = t2 * cc[2] * cc[0];
	complex<double> t66 = bb[2] * bb[1];
	complex<double> t71 = t2 * cc[1] * cc[0];
	complex<double> t75 = t2 * bb[0];
	t76 = t75 * cc[0];
	complex<double> t79 = bb[0] * bb[2];
	t83 = t2 * t22;
	complex<double> t84 = cc[2] * dd[0];
	t87 = t2 * t28;
	t88 = cc[2] * dd[1];
	t91 = t2 * t45;
	t92 = cc[1] * dd[0];
	complex<double> t95 = cc[1] * dd[2];
	complex<double> t98 = cc[0] * dd[2];
	complex<double> t101 = cc[0] * dd[1];
	complex<double> t104 = t2 * t4;
	complex<double> t107 = t2 * t9;
	t114 = t2 * t17;
	complex<double> t119 = -t57 * t15 / 0.6e1 + t61 * t33 / 0.6e1 + t65 * t66 * dd[2] / 0.6e1 - t71 * t66 * dd[1] / 0.6e1 - t76 * t7 / 0.6e1 + t71 * t79 * dd[0] / 0.6e1 + t83 * t84 / 0.12e2 - t87 * t88 / 0.12e2 - t91 * t92 / 0.12e2 + t87 * t95 / 0.12e2 - t83 * t98 / 0.12e2 + t91 * t101 / 0.12e2 - t104 * t33 / 0.12e2 + t107 * t15 / 0.12e2 - t104 * t36 / 0.12e2 + t107 * t20 / 0.12e2 - t114 * t12 / 0.12e2 + t104 * t15 / 0.12e2;
	complex<double> t132 = pow(cc[1], 0.2e1);
	complex<double> t138 = cc[0] * cc[2] * dd[1];
	t141 = pow(cc[2], 0.2e1);
	t147 = cc[0] * cc[1] * dd[2];
	complex<double> t150 = pow(cc[0], 0.2e1);
	complex<double> t158 = cc[1] * cc[2] * dd[0];
	t177 = t104 * t7 / 0.12e2 - t107 * t12 / 0.12e2 - t107 * t33 / 0.12e2 - t114 * t36 / 0.12e2 + t114 * t7 / 0.12e2 + t114 * t20 / 0.12e2 - t60 * t132 * bb[0] * dd[2] / 0.6e1 - t104 * t138 / 0.6e1 + t56 * t141 * bb[0] * dd[1] / 0.6e1 + t107 * t147 / 0.6e1 + t75 * t150 * bb[1] * dd[2] / 0.6e1 - t114 * t147 / 0.6e1 + t104 * t158 / 0.6e1 - t56 * t141 * bb[1] * dd[0] / 0.6e1 + t60 * t132 * bb[2] * dd[0] / 0.6e1 - t107 * t158 / 0.6e1 + t114 * t138 / 0.6e1 - t75 * t150 * bb[2] * dd[1] / 0.6e1;
	t178 = t119 + t177;
	complex<double> t191 = dd[1] * t1;
	complex<double> t195 = -t8 / 0.72e2 + t13 / 0.72e2 - t16 / 0.72e2 - t21 / 0.72e2 + t26 / 0.72e2 - t27 / 0.72e2 + t34 / 0.72e2 + t37 / 0.72e2 - t41 / 0.72e2 + t52 / 0.72e2 - t53 / 0.72e2 + t54 / 0.72e2 + t71 * t66 * t191 / 0.18e2;
	complex<double> t196 = dd[0] * t1;
	complex<double> t203 = dd[2] * t1;
	complex<double> t219 = t150 * cc[0];
	complex<double> t224 = t141 * cc[2];
	complex<double> t229 = t132 * cc[1];
	complex<double> t234 = t76 * t6 * t196 / 0.18e2 - t71 * t79 * t196 / 0.18e2 - t65 * t66 * t203 / 0.18e2 - t61 * t32 * t191 / 0.18e2 + t57 * t14 * t203 / 0.18e2 - t25 / 0.72e2 - t31 / 0.72e2 + t40 / 0.72e2 + t44 / 0.72e2 - t48 / 0.72e2 + t51 / 0.72e2 - t3 * t219 * bb[2] * dd[1] / 0.18e2 + t3 * t224 * bb[0] * dd[1] / 0.18e2 - t3 * t229 * bb[0] * dd[2] / 0.18e2;
	complex<double> t248 = t56 * t141;
	t252 = t60 * t132;
	t265 = t75 * t150;
	complex<double> t278 = t3 * t219 * bb[1] * dd[2] - t3 * t224 * bb[1] * dd[0] + t3 * t229 * bb[2] * dd[0] - t107 * cc[0] * t95 * t1 - t248 * bb[0] * dd[1] * t1 + t252 * bb[0] * dd[2] * t1 + t104 * cc[0] * t88 * t1 + t107 * cc[2] * t92 * t1 - t114 * cc[2] * t101 * t1 + t265 * bb[2] * dd[1] * t1 - t265 * bb[1] * dd[2] * t1 + t114 * cc[1] * t98 * t1 - t104 * cc[1] * t84 * t1;
	complex<double> t285 = t3 * t141;
	t287 = t3 * t150;
	complex<double> t294 = t3 * t132;
	complex<double> t300 = t248 * bb[1] * dd[0] * t1 - t252 * bb[2] * dd[0] * t1 - t285 * t20 - t287 * t15 - t285 * t15 + t285 * t12 + t287 * t12 + t287 * t33 + t285 * t36 - t294 * t20 - t294 * t7 - t287 * t7 + t294 * t33 + t294 * t36;
	//
	cm[0] = t55 / 0.24e2;
	cm[1] = -t55 / 0.12e2;
	cm[2] = t178;
	cm[3] = -t178;
	cm[4] = t55 / 0.24e2;
	cm[5] = t195 + t234 + t278 / 0.18e2 + t300 / 0.18e2;
	//
	t1 =  ko *  ko;
	t3 =  j * j;
	t5 = 0.1e1 /  t1 /  t3;
	t7 = 3.0 * t5 * cm[3];
	t10 = 0.1e1 /  ko /  j;
	t14 = 8.0 * t5 * cm[5];
	t16 = 2.0 * t10 * cm[0];
	t18 = 4.0 * t10 * cm[5];
	t20 = 3.0 * t5 * cm[2];
	t22 = 2.0 * t10 * cm[1];
	complex<double> t24 = 2.0 * t10 * cm[4];
	t31 = 8.0 / t1 / ko / t3 / j * cm[5];
	//
	coefm[0] = cm[1];
	coefm[1] = t7;
	coefm[2] = 3.0 * t10 * cm[2];
	coefm[3] = cm[0];
	coefm[4] = cm[2];
	coefm[5] = t14;
	coefm[6] = t16;
	coefm[7] = t18;
	coefm[8] = t20;
	coefm[9] = t22;
	coefm[10] = -t14;
	coefm[11] = -t7;
	coefm[12] = -t18;
	coefm[13] = cm[4];
	coefm[14] = cm[5];
	coefm[15] = -cm[5];
	coefm[16] = t24;
	coefm[17] = t31;
	coefm[18] = -t20;
	coefm[19] = -t22;
	coefm[20] = -t16;
	coefm[21] = -t31;
	coefm[22] = -t24;
	coefm[23] = -t31;
	coefm[24] = t31;
	coefm[25] = cm[3];
	coefm[26] = 3.0 * t10 * cm[3];
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g1_f3
// ***********************************************************************

void coefficients_g1_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = t1 * bb[2];
	complex<double> t3 = t2 * cc[2];
	complex<double> t4 = cc[1] * bb[0];
	complex<double> t5 = t4 * dd[2];
	complex<double> t6 = t3 * t5;
	complex<double> t8 = t1 * bb[1];
	complex<double> t9 = t8 * cc[1];
	complex<double> t10 = cc[2] * bb[0];
	complex<double> t11 = t10 * dd[1];
	complex<double> t12 = t9 * t11;
	complex<double> t15 = t1 * cc[2] * cc[0];
	complex<double> t16 = bb[1] * bb[2];
	complex<double> t18 = t15 * t16 * dd[2];
	complex<double> t21 = t1 * cc[1] * cc[0];
	complex<double> t23 = t21 * t16 * dd[1];
	complex<double> t25 = t1 * bb[0];
	complex<double> t26 = t25 * cc[0];
	complex<double> t27 = cc[2] * bb[1];
	complex<double> t28 = t27 * dd[0];
	complex<double> t29 = t26 * t28;
	complex<double> t31 = bb[2] * bb[0];
	complex<double> t33 = t21 * t31 * dd[0];
	complex<double> t35 = pow(bb[2], 0.2e1);
	complex<double> t37 = t1 * t35 * bb[2];
	complex<double> t38 = cc[0] * dd[1];
	complex<double> t39 = t37 * t38;
	complex<double> t40 = t39 / 0.24e2;
	complex<double> t41 = pow(bb[1], 0.2e1);
	complex<double> t43 = t1 * t41 * bb[1];
	complex<double> t44 = cc[0] * dd[2];
	complex<double> t45 = t43 * t44;
	complex<double> t46 = t45 / 0.24e2;
	complex<double> t47 = pow(bb[0], 0.2e1);
	complex<double> t49 = t1 * t47 * bb[0];
	complex<double> t50 = cc[1] * dd[2];
	complex<double> t51 = t49 * t50;
	complex<double> t52 = t51 / 0.24e2;
	complex<double> t53 = cc[1] * dd[0];
	complex<double> t54 = t37 * t53;
	complex<double> t55 = t54 / 0.24e2;
	complex<double> t56 = cc[2] * dd[0];
	complex<double> t57 = t43 * t56;
	complex<double> t58 = t57 / 0.24e2;
	complex<double> t59 = cc[2] * dd[1];
	complex<double> t60 = t49 * t59;
	complex<double> t61 = t60 / 0.24e2;
	complex<double> t62 = t1 * t35;
	complex<double> t64 = cc[0] * bb[1] * dd[2];
	complex<double> t65 = t62 * t64;
	complex<double> t66 = t65 / 0.24e2;
	complex<double> t68 = cc[0] * cc[2] * dd[1];
	complex<double> t69 = t62 * t68;
	complex<double> t71 = t1 * t41;
	complex<double> t73 = cc[0] * bb[2] * dd[1];
	complex<double> t74 = t71 * t73;
	complex<double> t75 = t74 / 0.24e2;
	complex<double> t77 = cc[0] * cc[1] * dd[2];
	complex<double> t78 = t71 * t77;
	complex<double> t80 = t71 * t5;
	complex<double> t81 = t80 / 0.24e2;
	complex<double> t82 = pow(cc[1], 0.2e1);
	complex<double> t85 = t8 * t82 * bb[0] * dd[2];
	complex<double> t87 = -t6 / 0.12e2 + t12 / 0.12e2 + t18 / 0.12e2 - t23 / 0.12e2 - t29 / 0.12e2 + t33 / 0.12e2 + t40 - t46 + t52 - t55 + t58 - t61 - t66 - t69 / 0.12e2 + t75 + t78 / 0.12e2 + t81 - t85 / 0.12e2;
	complex<double> t88 = t1 * t47;
	complex<double> t89 = t88 * t64;
	complex<double> t90 = t89 / 0.24e2;
	complex<double> t91 = pow(cc[0], 0.2e1);
	complex<double> t94 = t25 * t91 * bb[1] * dd[2];
	complex<double> t97 = cc[1] * bb[2] * dd[0];
	complex<double> t98 = t88 * t97;
	complex<double> t99 = t98 / 0.24e2;
	complex<double> t100 = t88 * t77;
	complex<double> t102 = t62 * t5;
	complex<double> t103 = t102 / 0.24e2;
	complex<double> t105 = cc[1] * cc[2] * dd[0];
	complex<double> t106 = t62 * t105;
	complex<double> t108 = t62 * t28;
	complex<double> t109 = t108 / 0.24e2;
	complex<double> t110 = pow(cc[2], 0.2e1);
	complex<double> t113 = t2 * t110 * bb[1] * dd[0];
	complex<double> t115 = t71 * t97;
	complex<double> t116 = t115 / 0.24e2;
	complex<double> t119 = t8 * t82 * bb[2] * dd[0];
	complex<double> t121 = t71 * t11;
	complex<double> t122 = t121 / 0.24e2;
	complex<double> t123 = t71 * t105;
	complex<double> t125 = t88 * t28;
	complex<double> t126 = t125 / 0.24e2;
	complex<double> t127 = t88 * t68;
	complex<double> t129 = t88 * t73;
	complex<double> t130 = t129 / 0.24e2;
	complex<double> t133 = t25 * t91 * bb[2] * dd[1];
	complex<double> t137 = t2 * t110 * bb[0] * dd[1];
	complex<double> t139 = t62 * t11;
	complex<double> t140 = t139 / 0.24e2;
	complex<double> t141 = -t90 + t94 / 0.12e2 - t99 - t100 / 0.12e2 + t103 + t106 / 0.12e2 + t109 - t113 / 0.12e2 - t116 + t119 / 0.12e2 - t122 - t123 / 0.12e2 + t126 + t127 / 0.12e2 + t130 - t133 / 0.12e2 + t137 / 0.12e2 - t140;
	complex<double> t143 = sqrt(0.3e1);
	complex<double> t144 = dd[2] * t143;
	complex<double> t146 = t15 * t16 * t144;
	complex<double> t147 = dd[1] * t143;
	complex<double> t149 = t9 * t10 * t147;
	complex<double> t151 = t21 * t16 * t147;
	complex<double> t153 = t3 * t4 * t144;
	complex<double> t154 = dd[0] * t143;
	complex<double> t156 = t26 * t27 * t154;
	complex<double> t158 = t21 * t31 * t154;
	complex<double> t159 = t71 * cc[0];
	complex<double> t160 = t50 * t143;
	complex<double> t161 = t159 * t160;
	complex<double> t162 = t8 * t82;
	complex<double> t163 = bb[0] * dd[2];
	complex<double> t164 = t163 * t143;
	complex<double> t165 = t162 * t164;
	complex<double> t166 = t62 * cc[0];
	complex<double> t167 = t59 * t143;
	complex<double> t168 = t166 * t167;
	complex<double> t169 = t2 * t110;
	complex<double> t170 = bb[0] * dd[1];
	complex<double> t171 = t170 * t143;
	complex<double> t172 = t169 * t171;
	complex<double> t173 = t88 * cc[1];
	complex<double> t174 = t44 * t143;
	complex<double> t175 = t173 * t174;
	complex<double> t176 = t62 * cc[1];
	complex<double> t177 = t56 * t143;
	complex<double> t178 = t176 * t177;
	complex<double> t179 = bb[1] * dd[0];
	complex<double> t180 = t179 * t143;
	complex<double> t181 = t169 * t180;
	complex<double> t182 = -t146 - t149 + t151 + t153 + t156 - t158 - t161 + t165 + t168 - t172 + t175 - t178 + t181;
	complex<double> t183 = bb[2] * dd[0];
	complex<double> t184 = t183 * t143;
	complex<double> t185 = t162 * t184;
	complex<double> t187 = t88 * cc[2];
	complex<double> t188 = t38 * t143;
	complex<double> t189 = t187 * t188;
	complex<double> t191 = t25 * t91;
	complex<double> t192 = bb[2] * dd[1];
	complex<double> t193 = t192 * t143;
	complex<double> t194 = t191 * t193;
	complex<double> t196 = bb[1] * dd[2];
	complex<double> t197 = t196 * t143;
	complex<double> t198 = t191 * t197;
	complex<double> t200 = t71 * cc[2];
	complex<double> t201 = t53 * t143;
	complex<double> t202 = t200 * t201;
	complex<double> t204 = t71 * cc[1];
	complex<double> t205 = t204 * t164;
	complex<double> t207 = t88 * cc[0];
	complex<double> t208 = t207 * t193;
	complex<double> t210 = t159 * t193;
	complex<double> t212 = t207 * t197;
	complex<double> t214 = t62 * cc[2];
	complex<double> t215 = t214 * t171;
	complex<double> t217 = t187 * t180;
	complex<double> t219 = t204 * t184;
	complex<double> t221 = t166 * t197;
	complex<double> t223 = t173 * t184;
	complex<double> t225 = -t185 / 0.18e2 - t189 / 0.18e2 + t194 / 0.18e2 - t198 / 0.18e2 + t202 / 0.18e2 - t205 / 0.72e2 - t208 / 0.72e2 - t210 / 0.72e2 + t212 / 0.72e2 + t215 / 0.72e2 - t217 / 0.72e2 + t219 / 0.72e2 + t221 / 0.72e2 + t223 / 0.72e2;
	complex<double> t227 = t200 * t171;
	complex<double> t229 = t176 * t164;
	complex<double> t231 = t214 * t180;
	complex<double> t233 = t143 * t1;
	complex<double> t234 = t233 * t82;
	complex<double> t239 = t233 * t91;
	complex<double> t246 = t233 * t110;
	complex<double> t257 = t227 / 0.72e2 - t229 / 0.72e2 - t231 / 0.72e2 - t234 * t73 / 0.18e2 + t234 * t64 / 0.18e2 + t239 * t11 / 0.18e2 - t234 * t28 / 0.18e2 - t239 * t5 / 0.18e2 + t246 * t64 / 0.18e2 + t246 * t97 / 0.18e2 - t239 * t28 / 0.18e2 - t246 * t5 / 0.18e2 - t246 * t73 / 0.18e2;
	complex<double> t262 = t49 * t160;
	complex<double> t264 = t49 * t167;
	complex<double> t266 = t43 * t177;
	complex<double> t268 = t37 * t201;
	complex<double> t270 = t37 * t188;
	complex<double> t272 = t43 * t174;
	complex<double> t274 = t91 * cc[0];
	complex<double> t279 = t82 * cc[1];
	complex<double> t284 = t110 * cc[2];
	complex<double> t301 = t239 * t97 / 0.18e2 + t234 * t11 / 0.18e2 - t262 / 0.72e2 + t264 / 0.72e2 - t266 / 0.72e2 + t268 / 0.72e2 - t270 / 0.72e2 + t272 / 0.72e2 + t233 * t274 * bb[1] * dd[2] / 0.18e2 + t233 * t279 * bb[2] * dd[0] / 0.18e2 - t233 * t284 * bb[1] * dd[0] / 0.18e2 - t233 * t274 * bb[2] * dd[1] / 0.18e2 + t233 * t284 * bb[0] * dd[1] / 0.18e2 - t233 * t279 * bb[0] * dd[2] / 0.18e2;
	complex<double> t304 = -t217 - t208 - t266 + t215 - t205 - t270 - t210 + t221 + t268 + t212 - t231 - t229 - t262 + t219 + t264 + t272 + t223 + t227;
	complex<double> t305 = t80 - t54 + t51 - t60 + t57 - t121 + t39 - t45 - t115 - t139 + t125 - t98 + t102 - t65 - t89 + t129 + t108 + t74;
	complex<double> t306 = -t146 - t149 + t151 + t153 + t156 - t158 - t161 + t165 + t168 - t172 + t175 - t178 + t181 - t185 - t189 + t194 - t198 + t202;
	complex<double> t307 = t306 / 0.12e2 + t304 / 0.24e2;
	complex<double> t308 = t6 / 0.6e1;
	complex<double> t309 = t12 / 0.6e1;
	complex<double> t310 = t18 / 0.6e1;
	complex<double> t311 = t23 / 0.6e1;
	complex<double> t312 = t29 / 0.6e1;
	complex<double> t313 = t33 / 0.6e1;
	complex<double> t314 = t1 * t284;
	complex<double> t317 = -t308 + t309 + t310 - t311 - t312 + t313 + t40 - t46 + t52 - t55 + t58 - t61 - t314 * t170 / 0.6e1;
	complex<double> t318 = t1 * t279;
	complex<double> t321 = t1 * t274;
	complex<double> t330 = t69 / 0.6e1;
	complex<double> t331 = t78 / 0.6e1;
	complex<double> t332 = t85 / 0.6e1;
	complex<double> t333 = t94 / 0.6e1;
	complex<double> t334 = t318 * t163 / 0.6e1 - t321 * t196 / 0.6e1 + t314 * t179 / 0.6e1 - t318 * t183 / 0.6e1 + t321 * t192 / 0.6e1 - t66 - t330 + t75 + t331 + t81 - t332 - t90 + t333 - t99;
	complex<double> t336 = t100 / 0.6e1;
	complex<double> t337 = t106 / 0.6e1;
	complex<double> t338 = t113 / 0.6e1;
	complex<double> t339 = t119 / 0.6e1;
	complex<double> t340 = t123 / 0.6e1;
	complex<double> t341 = t127 / 0.6e1;
	complex<double> t342 = t133 / 0.6e1;
	complex<double> t343 = -t336 + t103 + t337 + t109 - t338 - t116 + t339 - t122 - t340 + t126 + t341 + t130 - t342;
	complex<double> t344 = t137 / 0.6e1;
	complex<double> t345 = t1 * t110;
	complex<double> t350 = t1 * t82;
	complex<double> t355 = t1 * t91;
	complex<double> t372 = t344 - t140 + t345 * t5 / 0.6e1 + t345 * t73 / 0.6e1 - t350 * t64 / 0.6e1 - t350 * t11 / 0.6e1 + t355 * t28 / 0.6e1 + t355 * t5 / 0.6e1 - t345 * t97 / 0.6e1 - t345 * t64 / 0.6e1 + t350 * t73 / 0.6e1 + t350 * t28 / 0.6e1 - t355 * t11 / 0.6e1 - t355 * t97 / 0.6e1;
	complex<double> t384 = -t308 + t309 + t310 - t311 - t312 + t313 + t39 / 0.12e2 - t45 / 0.12e2 + t51 / 0.12e2 - t54 / 0.12e2 + t57 / 0.12e2 - t60 / 0.12e2 - t65 / 0.12e2 - t330 + t74 / 0.12e2 + t331 + t80 / 0.12e2 - t332;
	complex<double> t394 = -t89 / 0.12e2 + t333 - t98 / 0.12e2 - t336 + t102 / 0.12e2 + t337 + t108 / 0.12e2 - t338 - t115 / 0.12e2 + t339 - t121 / 0.12e2 - t340 + t125 / 0.12e2 + t341 + t129 / 0.12e2 - t342 + t344 - t139 / 0.12e2;
	//
	c[0] = t87 + t141;
	c[1] = t182 / 0.18e2 + t225 + t257 + t301;
	c[2] = t304 / 0.24e2;
	c[3] = t304 / 0.24e2;
	c[4] = t305 / 0.8e1;
	c[5] = t307;
	c[6] = t307;
	c[7] = t307;
	c[8] = t317 + t334 + t343 + t372;
	c[9] = t384 + t394;
	c[10] = t305 / 0.8e1;
	//
	t3 = 0.1e1 /  ko /  j;
	t6 =  ko *  ko;
	t8 = j *  j;
	t10 = 0.1e1 /  t6 /  t8;
	t12 = 3.0 * t10 * c[5];
	complex<double> t14 = 4.0 * t3 * c[1];
	t16 = 2.0 * t3 * c[6];
	t18 = 2.0 * t3 * c[4];
	complex<double> t20 = 8.0 * t10 * c[1];
	complex<double> t22 = 3.0 * t10 * c[0];
	complex<double> t24 = 3.0 * t10 * c[9];
	t26 = 2.0 * t3 * c[3];
	t28 = 2.0 * t3 * c[2];
	complex<double> t30 = 3.0 * t10 * c[8];
	complex<double> t32 = 2.0 * t3 * c[10];
	t39 = 8.0 / t6 / ko / t8 / j * c[1];
	t41 = 3.0 * t3 * c[8];
	t43 = 2.0 * t3 * c[7];
	//
	coef[0] = 3.0 * t3 * c[9];
	coef[1] = c[9];
	coef[2] = t12;
	coef[3] = t14;
	coef[4] = t16;
	coef[5] = t18;
	coef[6] = -t14;
	coef[7] = t20;
	coef[8] = t22;
	coef[9] = -t24;
	coef[10] = c[7];
	coef[11] = c[10];
	coef[12] = t26;
	coef[13] = t28;
	coef[14] = -t20;
	coef[15] = t30;
	coef[16] = -t30;
	coef[17] = -t32;
	coef[18] = -t18;
	coef[19] = -t22;
	coef[20] = -c[8];
	coef[21] = c[8];
	coef[22] = t30;
	coef[23] = -t26;
	coef[24] = -t28;
	coef[25] = t39;
	coef[26] = -t41;
	coef[27] = t32;
	coef[28] = -t39;
	coef[29] = t43;
	coef[30] = t41;
	coef[31] = -t16;
	coef[32] = -t12;
	coef[33] = -t39;
	coef[34] = c[6];
	coef[35] = c[2];
	coef[36] = 3.0 * t3 * c[5];
	coef[37] = 3.0 * t3 * c[0];
	coef[38] = t24;
	coef[39] = c[5];
	coef[40] = c[0];
	coef[41] = c[4];
	coef[42] = c[3];
	coef[43] = -c[1];
	coef[44] = c[1];
	coef[45] = -t43;
	coef[46] = -t30;
	coef[47] = t39;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[1], 0.2e1);
	t3 = t1 * t2;
	t5 = cc[2] * cc[1] * dd[0];
	t6 = t3 * t5;
	complex<double> t7 = t6 / 0.6e1;
	t8 = pow(bb[0], 0.2e1);
	t9 = t1 * t8;
	t11 = cc[2] * cc[0] * dd[1];
	t12 = t9 * t11;
	complex<double> t13 = t12 / 0.6e1;
	t14 = t1 * bb[0];
	t15 = pow(cc[0], 0.2e1);
	t18 = t14 * t15 * bb[2] * dd[1];
	complex<double> t19 = t18 / 0.6e1;
	t21 = cc[2] * bb[0] * dd[1];
	t22 = t3 * t21;
	t25 = cc[2] * bb[1] * dd[0];
	t26 = t9 * t25;
	t29 = cc[0] * bb[2] * dd[1];
	t30 = t9 * t29;
	t32 = cc[1] * bb[0];
	t33 = t32 * dd[2];
	complex<double> t34 = t3 * t33;
	complex<double> t36 = pow(bb[2], 0.2e1);
	t37 = t1 * t36;
	t38 = t37 * t25;
	t40 = t1 * bb[2];
	t41 = pow(cc[2], 0.2e1);
	t44 = t40 * t41 * bb[1] * dd[0];
	t45 = t44 / 0.6e1;
	t47 = cc[1] * cc[0] * dd[2];
	complex<double> t48 = t9 * t47;
	t49 = t48 / 0.6e1;
	t52 = t14 * t15 * bb[1] * dd[2];
	t53 = t52 / 0.6e1;
	t54 = t37 * t33;
	t56 = t37 * t5;
	t57 = t56 / 0.6e1;
	t59 = cc[1] * bb[2] * dd[0];
	t60 = t3 * t59;
	t62 = t1 * bb[1];
	complex<double> t63 = pow(cc[1], 0.2e1);
	t66 = t62 * t63 * bb[2] * dd[0];
	complex<double> t67 = t66 / 0.6e1;
	t68 = t9 * t59;
	t71 = cc[0] * bb[1] * dd[2];
	complex<double> t72 = t9 * t71;
	complex<double> t76 = t62 * t63 * bb[0] * dd[2];
	t77 = t76 / 0.6e1;
	t78 = t7 - t13 + t19 + t22 / 0.12e2 - t26 / 0.12e2 - t30 / 0.12e2 - t34 / 0.12e2 - t38 / 0.12e2 + t45 + t49 - t53 - t54 / 0.12e2 - t57 + t60 / 0.12e2 - t67 + t68 / 0.12e2 + t72 / 0.12e2 + t77;
	t81 = t40 * t41 * bb[0] * dd[1];
	t82 = t81 / 0.6e1;
	complex<double> t83 = t37 * t11;
	complex<double> t84 = t83 / 0.6e1;
	t85 = t37 * t21;
	t87 = t37 * t71;
	t89 = t3 * t47;
	t90 = t89 / 0.6e1;
	t91 = t3 * t29;
	t94 = t1 * t2 * bb[1];
	complex<double> t95 = cc[0] * dd[2];
	complex<double> t96 = t94 * t95;
	t99 = t1 * t8 * bb[0];
	t100 = cc[1] * dd[2];
	complex<double> t101 = t99 * t100;
	complex<double> t104 = t1 * t36 * bb[2];
	t105 = cc[1] * dd[0];
	t106 = t104 * t105;
	t108 = cc[2] * dd[0];
	t109 = t94 * t108;
	complex<double> t111 = cc[2] * dd[1];
	complex<double> t112 = t99 * t111;
	complex<double> t114 = cc[0] * dd[1];
	t115 = t104 * t114;
	complex<double> t117 = t1 * cc[1];
	complex<double> t118 = t117 * cc[2];
	t119 = bb[1] * bb[0];
	t121 = t118 * t119 * dd[1];
	t122 = t121 / 0.6e1;
	complex<double> t124 = t1 * cc[0] * cc[2];
	t126 = t124 * t119 * dd[0];
	t127 = t126 / 0.6e1;
	complex<double> t128 = t40 * cc[2];
	t129 = t128 * t33;
	t130 = t129 / 0.6e1;
	complex<double> t131 = bb[2] * bb[1];
	t133 = t124 * t131 * dd[2];
	complex<double> t134 = t133 / 0.6e1;
	complex<double> t135 = t117 * cc[0];
	t137 = t135 * t131 * dd[1];
	complex<double> t138 = t137 / 0.6e1;
	t139 = bb[0] * bb[2];
	t141 = t135 * t139 * dd[0];
	complex<double> t142 = t141 / 0.6e1;
	t143 = -t82 + t84 + t85 / 0.12e2 + t87 / 0.12e2 - t90 - t91 / 0.12e2 + t96 / 0.12e2 - t101 / 0.12e2 + t106 / 0.12e2 - t109 / 0.12e2 + t112 / 0.12e2 - t115 / 0.12e2 - t122 + t127 + t130 - t134 + t138 - t142;
	complex<double> t145 = -t26 - t38 - t115 + t85 - t30 + t112 + t106 - t109 - t91 + t68 - t101 + t96 + t22 - t54 + t87 + t72 + t60 - t34;
	t146 = sqrt(0.3e1);
	t147 = dd[0] * t146;
	t149 = t124 * t119 * t147;
	complex<double> t150 = dd[2] * t146;
	complex<double> t152 = t124 * t131 * t150;
	t153 = dd[1] * t146;
	complex<double> t155 = t135 * t131 * t153;
	complex<double> t157 = t135 * t139 * t147;
	t159 = t118 * t119 * t153;
	t161 = t128 * t32 * t150;
	t162 = t14 * t15;
	t163 = bb[2] * dd[1];
	t164 = t163 * t146;
	t165 = t162 * t164;
	t166 = t62 * t63;
	t167 = bb[0] * dd[2];
	t168 = t167 * t146;
	t169 = t166 * t168;
	t170 = t9 * cc[1];
	t171 = t95 * t146;
	t172 = t170 * t171;
	t173 = t37 * cc[1];
	t174 = t108 * t146;
	t175 = t173 * t174;
	t176 = bb[1] * dd[2];
	t177 = t176 * t146;
	t178 = t162 * t177;
	t179 = t9 * cc[2];
	t180 = t114 * t146;
	t181 = t179 * t180;
	t182 = t40 * t41;
	t183 = bb[1] * dd[0];
	t184 = t183 * t146;
	t185 = t182 * t184;
	complex<double> t186 = bb[2] * dd[0];
	t187 = t186 * t146;
	t188 = t166 * t187;
	t189 = bb[0] * dd[1];
	complex<double> t190 = t189 * t146;
	t191 = t182 * t190;
	t192 = t3 * cc[2];
	t193 = t105 * t146;
	t194 = t192 * t193;
	complex<double> t195 = t3 * cc[0];
	t196 = t100 * t146;
	t197 = t195 * t196;
	t198 = t37 * cc[0];
	complex<double> t199 = t111 * t146;
	t200 = t198 * t199;
	t201 = t149 - t152 + t155 - t157 - t159 + t161 + t165 + t169 + t172 - t175 - t178 - t181 + t185 - t188 - t191 + t194 - t197 + t200;
	t202 = t99 * t196;
	complex<double> t203 = t104 * t180;
	t204 = t192 * t190;
	t205 = t94 * t174;
	complex<double> t206 = t3 * cc[1];
	t207 = t206 * t168;
	t208 = t9 * cc[0];
	complex<double> t209 = t208 * t177;
	t210 = t173 * t168;
	complex<double> t211 = t195 * t164;
	t212 = t99 * t199;
	complex<double> t213 = t170 * t187;
	t214 = t104 * t193;
	t215 = t198 * t177;
	complex<double> t216 = t94 * t171;
	t217 = t37 * cc[2];
	complex<double> t218 = t217 * t184;
	t219 = t217 * t190;
	complex<double> t220 = t206 * t187;
	t221 = t208 * t164;
	complex<double> t222 = t179 * t184;
	t223 = -t202 - t203 + t204 - t205 - t207 + t209 - t210 - t211 + t212 + t213 + t214 + t215 + t216 - t218 + t219 + t220 - t221 - t222;
	complex<double> t224 = t201 / 0.12e2 + t223 / 0.24e2;
	t225 = t22 / 0.24e2;
	complex<double> t226 = t26 / 0.24e2;
	t227 = t30 / 0.24e2;
	complex<double> t228 = t34 / 0.24e2;
	t229 = t38 / 0.24e2;
	complex<double> t230 = t54 / 0.24e2;
	t231 = -t7 + t13 - t19 - t225 + t226 + t227 + t228 + t229 - t45 - t49 + t53 + t230 + t57;
	complex<double> t232 = t60 / 0.24e2;
	t233 = t68 / 0.24e2;
	t234 = t72 / 0.24e2;
	complex<double> t235 = t85 / 0.24e2;
	complex<double> t236 = t87 / 0.24e2;
	complex<double> t237 = t91 / 0.24e2;
	complex<double> t238 = t1 * t63;
	complex<double> t243 = t1 * t15;
	t246 = -t232 + t67 - t233 - t234 - t77 + t82 - t84 - t235 - t236 + t90 + t237 + t238 * t29 / 0.6e1 + t238 * t25 / 0.6e1 - t243 * t21 / 0.6e1;
	complex<double> t252 = t1 * t41;
	complex<double> t267 = t96 / 0.24e2;
	t268 = t101 / 0.24e2;
	complex<double> t269 = t106 / 0.24e2;
	t270 = t109 / 0.24e2;
	complex<double> t271 = -t243 * t59 / 0.6e1 + t243 * t25 / 0.6e1 - t252 * t59 / 0.6e1 + t243 * t33 / 0.6e1 - t252 * t71 / 0.6e1 + t252 * t33 / 0.6e1 + t252 * t29 / 0.6e1 - t238 * t71 / 0.6e1 - t238 * t21 / 0.6e1 - t267 + t268 - t269 + t270;
	t272 = t112 / 0.24e2;
	complex<double> t273 = t115 / 0.24e2;
	t274 = t41 * cc[2];
	complex<double> t275 = t1 * t274;
	complex<double> t278 = t63 * cc[1];
	t279 = t1 * t278;
	complex<double> t282 = t15 * cc[0];
	complex<double> t283 = t1 * t282;
	complex<double> t292 = -t272 + t273 - t275 * t189 / 0.6e1 + t279 * t167 / 0.6e1 - t283 * t176 / 0.6e1 + t275 * t183 / 0.6e1 - t279 * t186 / 0.6e1 + t283 * t163 / 0.6e1 + t122 - t127 - t130 + t134 - t138 + t142;
	complex<double> t295 = t149 - t152 + t155 - t157 - t159 + t161 + t165 + t169 + t172 - t175 - t178 - t181 + t185;
	t310 = -t188 / 0.18e2 - t191 / 0.18e2 + t194 / 0.18e2 - t197 / 0.18e2 + t200 / 0.18e2 + t204 / 0.72e2 - t210 / 0.72e2 + t209 / 0.72e2 - t211 / 0.72e2 + t213 / 0.72e2 + t215 / 0.72e2 - t207 / 0.72e2 + t220 / 0.72e2 - t221 / 0.72e2;
	complex<double> t315 = t146 * t1;
	complex<double> t316 = t315 * t15;
	complex<double> t319 = t315 * t41;
	complex<double> t322 = t315 * t63;
	t339 = t219 / 0.72e2 - t222 / 0.72e2 - t218 / 0.72e2 + t316 * t59 / 0.18e2 + t319 * t59 / 0.18e2 - t322 * t29 / 0.18e2 - t319 * t29 / 0.18e2 + t319 * t71 / 0.18e2 - t319 * t33 / 0.18e2 - t316 * t33 / 0.18e2 - t316 * t25 / 0.18e2 + t322 * t21 / 0.18e2 + t316 * t21 / 0.18e2;
	complex<double> t374 = -t322 * t25 / 0.18e2 + t322 * t71 / 0.18e2 + t212 / 0.72e2 + t214 / 0.72e2 + t216 / 0.72e2 - t203 / 0.72e2 - t205 / 0.72e2 - t202 / 0.72e2 + t315 * t282 * bb[1] * dd[2] / 0.18e2 - t315 * t274 * bb[1] * dd[0] / 0.18e2 + t315 * t274 * bb[0] * dd[1] / 0.18e2 - t315 * t278 * bb[0] * dd[2] / 0.18e2 - t315 * t282 * bb[2] * dd[1] / 0.18e2 + t315 * t278 * bb[2] * dd[0] / 0.18e2;
	complex<double> t386 = -t6 / 0.12e2 + t12 / 0.12e2 - t18 / 0.12e2 - t225 + t226 + t227 + t228 + t229 - t44 / 0.12e2 - t48 / 0.12e2 + t52 / 0.12e2 + t230 + t56 / 0.12e2 - t232 + t66 / 0.12e2 - t233 - t234 - t76 / 0.12e2;
	complex<double> t396 = t81 / 0.12e2 - t83 / 0.12e2 - t235 - t236 + t89 / 0.12e2 + t237 - t267 + t268 - t269 + t270 - t272 + t273 + t121 / 0.12e2 - t126 / 0.12e2 - t129 / 0.12e2 + t133 / 0.12e2 - t137 / 0.12e2 + t141 / 0.12e2;
	//
	cm[0] = t78 + t143;
	cm[1] = t145 / 0.8e1;
	cm[2] = t224;
	cm[3] = t231 + t246 + t271 + t292;
	cm[4] = t223 / 0.24e2;
	cm[5] = -t224;
	cm[6] = -t223 / 0.24e2;
	cm[7] = t295 / 0.18e2 + t310 + t339 + t374;
	cm[8] = t386 + t396;
	cm[9] = -t145 / 0.8e1;
	cm[10] = -t224;
	//
	t3 = 0.1e1 /  ko /  j;
	t6 =  ko *  ko;
	t8 =  j *  j;
	t10 = 0.1e1 /  t6 /  t8;
	t12 = 3.0 * t10 * cm[3];
	t14 = 2.0 * t3 * cm[5];
	t16 = 2.0 * t3 * cm[6];
	t18 = 4.0 * t3 * cm[7];
	t20 = 2.0 * t3 * cm[4];
	t22 = 3.0 * t10 * cm[10];
	t24 = 3.0 * t10 * cm[0];
	t26 = 8.0 * t10 * cm[7];
	t28 = 3.0 * t10 * cm[8];
	t30 = 2.0 * t3 * cm[9];
	t37 = 8.0 / t6 / ko / t8 / j * cm[7];
	t39 = 2.0 * t3 * cm[1];
	t41 = 3.0 * t3 * cm[3];
	t43 = 2.0 * t3 * cm[2];
	//
	coefm[0] = cm[0];
	coefm[1] = 3.0 * t3 * cm[0];
	coefm[2] = t12;
	coefm[3] = t14;
	coefm[4] = t16;
	coefm[5] = -t18;
	coefm[6] = cm[2];
	coefm[7] = cm[1];
	coefm[8] = -cm[7];
	coefm[9] = cm[7];
	coefm[10] = t20;
	coefm[11] = t22;
	coefm[12] = t18;
	coefm[13] = -t24;
	coefm[14] = -t26;
	coefm[15] = t28;
	coefm[16] = t30;
	coefm[17] = t26;
	coefm[18] = -t12;
	coefm[19] = -t22;
	coefm[20] = -t20;
	coefm[21] = t37;
	coefm[22] = t39;
	coefm[23] = -t28;
	coefm[24] = -t41;
	coefm[25] = -t37;
	coefm[26] = -t14;
	coefm[27] = -cm[3];
	coefm[28] = cm[3];
	coefm[29] = t41;
	coefm[30] = t43;
	coefm[31] = -t30;
	coefm[32] = -t16;
	coefm[33] = -t39;
	coefm[34] = -t12;
	coefm[35] = t12;
	coefm[36] = t37;
	coefm[37] = -t37;
	coefm[38] = -t43;
	coefm[39] = 3.0 * t3 * cm[10];
	coefm[40] = 3.0 * t3 * cm[8];
	coefm[41] = t24;
	coefm[42] = cm[9];
	coefm[43] = cm[5];
	coefm[44] = cm[8];
	coefm[45] = cm[6];
	coefm[46] = cm[4];
	coefm[47] = cm[10];
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g2_f1
// ***********************************************************************

void coefficients_g2_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[6], cm[6];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = pow(bb[1], 0.2e1);
	complex<double> t3 = t1 * t2;
	complex<double> t5 = bb[0] * dd[1];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t7 = t5 * t6;
	complex<double> t8 = t3 * cc[2] * t7;
	complex<double> t9 = t3 * cc[1];
	complex<double> t11 = bb[2] * dd[0] * t6;
	complex<double> t12 = t9 * t11;
	complex<double> t13 = pow(bb[2], 0.2e1);
	complex<double> t15 = t1 * t13 * bb[2];
	complex<double> t16 = cc[0] * dd[1];
	complex<double> t18 = t15 * t16 * t6;
	complex<double> t19 = cc[1] * dd[0];
	complex<double> t21 = t15 * t19 * t6;
	complex<double> t23 = t1 * t2 * bb[1];
	complex<double> t24 = cc[0] * dd[2];
	complex<double> t26 = t23 * t24 * t6;
	complex<double> t27 = t1 * t13;
	complex<double> t28 = t27 * cc[2];
	complex<double> t29 = bb[1] * dd[0];
	complex<double> t30 = t29 * t6;
	complex<double> t31 = t28 * t30;
	complex<double> t32 = pow(bb[0], 0.2e1);
	complex<double> t33 = t1 * t32;
	complex<double> t34 = t33 * cc[0];
	complex<double> t36 = bb[1] * dd[2] * t6;
	complex<double> t37 = t34 * t36;
	complex<double> t38 = bb[0] * dd[2];
	complex<double> t39 = t38 * t6;
	complex<double> t40 = t9 * t39;
	complex<double> t42 = t1 * t32 * bb[0];
	complex<double> t43 = cc[1] * dd[2];
	complex<double> t45 = t42 * t43 * t6;
	complex<double> t47 = bb[2] * dd[1] * t6;
	complex<double> t48 = t34 * t47;
	complex<double> t50 = t33 * cc[1] * t11;
	complex<double> t52 = t27 * cc[0] * t36;
	complex<double> t53 = t28 * t7;
	complex<double> t55 = t27 * cc[1] * t39;
	complex<double> t56 = cc[2] * dd[1];
	complex<double> t58 = t42 * t56 * t6;
	complex<double> t60 = t33 * cc[2] * t30;
	complex<double> t61 = cc[2] * dd[0];
	complex<double> t63 = t23 * t61 * t6;
	complex<double> t65 = t3 * cc[0] * t47;
	complex<double> t66 = -t8 - t12 + t18 - t21 - t26 + t31 - t37 + t40 + t45 + t48 - t50 - t52 - t53 + t55 - t58 + t60 + t63 + t65;
	complex<double> t79 = cc[2] * bb[0];
	complex<double> t80 = t79 * dd[1];
	complex<double> t84 = cc[0] * bb[1] * dd[2];
	complex<double> t88 = cc[0] * bb[2] * dd[1];
	complex<double> t91 = cc[1] * bb[0];
	complex<double> t92 = t91 * dd[2];
	complex<double> t95 = pow(cc[2], 0.2e1);
	complex<double> t96 = t1 * t95;
	complex<double> t97 = t5 * bb[2];
	complex<double> t100 = t1 * cc[2];
	complex<double> t105 = t1 * cc[1];
	complex<double> t110 = pow(cc[1], 0.2e1);
	complex<double> t111 = t1 * t110;
	complex<double> t112 = t38 * bb[1];
	complex<double> t118 = cc[1] * bb[2] * dd[0];
	complex<double> t123 = cc[2] * bb[1];
	complex<double> t124 = t123 * dd[0];
	complex<double> t127 = -t15 * t16 / 0.12e2 + t23 * t24 / 0.12e2 - t42 * t43 / 0.12e2 + t15 * t19 / 0.12e2 - t23 * t61 / 0.12e2 + t42 * t56 / 0.12e2 + t27 * t80 / 0.12e2 + t27 * t84 / 0.12e2 - t3 * t88 / 0.12e2 - t3 * t92 / 0.12e2 - t96 * t97 / 0.6e1 + t100 * cc[0] * t13 * dd[1] / 0.6e1 - t105 * cc[0] * t2 * dd[2] / 0.6e1 + t111 * t112 / 0.6e1 + t33 * t84 / 0.12e2 + t33 * t118 / 0.12e2 - t27 * t92 / 0.12e2 - t27 * t124 / 0.12e2;
	complex<double> t128 = pow(cc[0], 0.2e1);
	complex<double> t129 = t1 * t128;
	complex<double> t132 = t1 * cc[0];
	complex<double> t141 = t29 * bb[2];
	complex<double> t165 = t1 * bb[2] * cc[2];
	complex<double> t169 = t1 * bb[1] * cc[1];
	complex<double> t172 = t100 * cc[0];
	complex<double> t173 = bb[1] * bb[2];
	complex<double> t177 = t105 * cc[0];
	complex<double> t182 = t1 * bb[0] * cc[0];
	complex<double> t185 = bb[2] * bb[0];
	complex<double> t189 = -t129 * t112 / 0.6e1 + t132 * cc[1] * t32 * dd[2] / 0.6e1 - t100 * cc[1] * t13 * dd[0] / 0.6e1 + t96 * t141 / 0.6e1 + t3 * t118 / 0.12e2 + t3 * t80 / 0.12e2 - t33 * t124 / 0.12e2 - t33 * t88 / 0.12e2 - t111 * t141 / 0.6e1 + t105 * cc[2] * t2 * dd[0] / 0.6e1 - t132 * cc[2] * t32 * dd[1] / 0.6e1 + t129 * t97 / 0.6e1 + t165 * t92 / 0.6e1 - t169 * t80 / 0.6e1 - t172 * t173 * dd[2] / 0.6e1 + t177 * t173 * dd[1] / 0.6e1 + t182 * t124 / 0.6e1 - t177 * t185 * dd[0] / 0.6e1;
	complex<double> t190 = t127 + t189;
	complex<double> t191 = t8 + t12 - t31 + t37 - t40 - t48 + t50 + t52 + t53 - t55 - t60 - t65 - t18;
	complex<double> t197 = t6 * t1;
	complex<double> t198 = t110 * cc[1];
	complex<double> t207 = t95 * cc[2];
	complex<double> t216 = t128 * cc[0];
	complex<double> t225 = dd[1] * t6;
	complex<double> t229 = dd[2] * t6;
	complex<double> t233 = dd[0] * t6;
	complex<double> t237 = t21 / 0.72e2 + t26 / 0.72e2 - t45 / 0.72e2 + t58 / 0.72e2 - t63 / 0.72e2 + t197 * t198 * bb[2] * dd[0] / 0.18e2 - t197 * t198 * bb[0] * dd[2] / 0.18e2 + t197 * t207 * bb[0] * dd[1] / 0.18e2 - t197 * t207 * bb[1] * dd[0] / 0.18e2 - t197 * t216 * bb[2] * dd[1] / 0.18e2 + t197 * t216 * bb[1] * dd[2] / 0.18e2 + t177 * t173 * t225 / 0.18e2 + t165 * t91 * t229 / 0.18e2 + t182 * t123 * t233 / 0.18e2;
	complex<double> t260 = t100 * cc[1];
	complex<double> t271 = -t177 * t185 * t233 - t169 * t79 * t225 - t172 * t173 * t229 - t177 * t2 * dd[2] * t6 + t111 * bb[0] * t36 - t96 * bb[0] * t47 + t172 * t13 * dd[1] * t6 - t129 * bb[1] * t39 + t177 * t32 * dd[2] * t6 - t260 * t13 * dd[0] * t6 + t96 * bb[1] * t11 - t111 * bb[2] * t30 + t260 * t2 * dd[0] * t6;
	complex<double> t277 = t197 * t128;
	complex<double> t279 = t197 * t110;
	complex<double> t281 = t197 * t95;
	complex<double> t292 = -t172 * t32 * dd[1] * t6 + t129 * bb[2] * t7 - t277 * t92 + t279 * t80 + t281 * t118 - t281 * t92 - t281 * t88 + t277 * t118 + t279 * t84 - t279 * t88 - t279 * t124 + t281 * t84 + t277 * t80 - t277 * t124;
	//
	c[0] = t66 / 0.12e2;
	c[1] = -t66 / 0.24e2;
	c[2] = -t66 / 0.24e2;
	c[3] = t190;
	c[4] = -t190;
	c[5] = t191 / 0.72e2 + t237 + t271 / 0.18e2 + t292 / 0.18e2;
	//
	t1 = ko *  ko;
	complex<double> t4 =  j *  j;
	t9 = 8.0 / t1 / ko / t4 / j * c[5];
	t12 = 1.0 / ko / j;
	complex<double> t17 = 1.0 / t1 / t4;
	t19 = 3.0 * t17 * c[4];
	t23 = 2.0 * t12 * c[1];
	complex<double> t25 = 2.0 * t12 * c[0];
	t27 = 3.0 * t17 * c[3];
	t29 = 2.0 * t12 * c[2];
	t31 = 4.0 * t12 * c[5];
	t33 = 8.0 * t17 * c[5];
	//
	coef[0] = -t9;
	coef[1] = c[4];
	coef[2] = 3.0 * t12 * c[4];
	coef[3] = c[0];
	coef[4] = c[3];
	coef[5] = c[1];
	coef[6] = t19;
	coef[7] = 3.0 * t12 * c[3];
	coef[8] = -t23;
	coef[9] = -t25;
	coef[10] = -t27;
	coef[11] = t9;
	coef[12] = -t9;
	coef[13] = t29;
	coef[14] = -t29;
	coef[15] = t9;
	coef[16] = c[5];
	coef[17] = c[2];
	coef[18] = -c[5];
	coef[19] = t31;
	coef[20] = t27;
	coef[21] = -t31;
	coef[22] = t33;
	coef[23] = t23;
	coef[24] = -t33;
	coef[25] = t25;
	coef[26] = -t19;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[1], 0.2e1);
	t3 = t1 * t2;
	t4 = t3 * cc[1];
	t6 = sqrt(0.3e1);
	t7 = bb[0] * dd[2] * t6;
	t8 = t4 * t7;
	t11 = bb[2] * dd[1] * t6;
	t12 = t3 * cc[0] * t11;
	t13 = pow(bb[2], 0.2e1);
	complex<double> t14 = t1 * t13;
	t17 = bb[1] * dd[2] * t6;
	t18 = t14 * cc[0] * t17;
	t19 = pow(bb[0], 0.2e1);
	t21 = t1 * t19 * bb[0];
	complex<double> t22 = cc[2] * dd[1];
	t24 = t21 * t22 * t6;
	t26 = t1 * t2 * bb[1];
	t27 = cc[2] * dd[0];
	t29 = t26 * t27 * t6;
	t31 = t1 * t13 * bb[2];
	t32 = cc[0] * dd[1];
	t34 = t31 * t32 * t6;
	complex<double> t35 = t14 * cc[2];
	t37 = bb[1] * dd[0] * t6;
	t38 = t35 * t37;
	t39 = t1 * t19;
	t42 = bb[2] * dd[0] * t6;
	t43 = t39 * cc[1] * t42;
	t45 = t14 * cc[1] * t7;
	complex<double> t46 = t39 * cc[0];
	t47 = t46 * t17;
	complex<double> t49 = bb[0] * dd[1] * t6;
	t50 = t35 * t49;
	t52 = t3 * cc[2] * t49;
	complex<double> t54 = t39 * cc[2] * t37;
	t55 = t46 * t11;
	t56 = t4 * t42;
	complex<double> t57 = cc[0] * dd[2];
	complex<double> t59 = t26 * t57 * t6;
	t60 = cc[1] * dd[2];
	complex<double> t62 = t21 * t60 * t6;
	t63 = cc[1] * dd[0];
	t65 = t31 * t63 * t6;
	t66 = -t8 - t12 + t18 + t24 - t29 - t34 - t38 + t43 - t45 + t47 + t50 + t52 - t54 - t55 + t56 + t59 - t62 + t65;
	complex<double> t74 = t1 * cc[0];
	complex<double> t75 = t74 * cc[1];
	complex<double> t76 = bb[0] * bb[2];
	complex<double> t77 = dd[0] * t6;
	complex<double> t82 = t1 * bb[1] * cc[1];
	complex<double> t83 = cc[2] * bb[0];
	t84 = dd[1] * t6;
	complex<double> t89 = t1 * bb[2] * cc[2];
	complex<double> t90 = cc[1] * bb[0];
	t91 = dd[2] * t6;
	t95 = bb[1] * bb[2];
	t100 = t1 * bb[0] * cc[0];
	complex<double> t101 = cc[2] * bb[1];
	t105 = t6 * t1;
	complex<double> t106 = pow(cc[2], 0.2e1);
	complex<double> t107 = t106 * cc[2];
	t112 = -t8 / 0.72e2 - t12 / 0.72e2 + t18 / 0.72e2 - t38 / 0.72e2 + t43 / 0.72e2 - t45 / 0.72e2 + t47 / 0.72e2 - t75 * t76 * t77 / 0.18e2 - t82 * t83 * t84 / 0.18e2 + t89 * t90 * t91 / 0.18e2 + t75 * t95 * t84 / 0.18e2 + t100 * t101 * t77 / 0.18e2 + t105 * t107 * bb[0] * dd[1] / 0.18e2;
	complex<double> t113 = pow(cc[1], 0.2e1);
	complex<double> t114 = t113 * cc[1];
	complex<double> t119 = pow(cc[0], 0.2e1);
	complex<double> t120 = t119 * cc[0];
	complex<double> t143 = t1 * cc[2];
	complex<double> t144 = t143 * cc[0];
	complex<double> t150 = t105 * t114 * bb[2] * dd[0] / 0.18e2 + t105 * t120 * bb[1] * dd[2] / 0.18e2 - t105 * t107 * bb[1] * dd[0] / 0.18e2 - t105 * t120 * bb[2] * dd[1] / 0.18e2 - t105 * t114 * bb[0] * dd[2] / 0.18e2 + t24 / 0.72e2 - t29 / 0.72e2 - t34 / 0.72e2 + t59 / 0.72e2 - t62 / 0.72e2 + t65 / 0.72e2 - t144 * t95 * t91 / 0.18e2 + t50 / 0.72e2 + t52 / 0.72e2;
	complex<double> t155 = t143 * cc[1];
	complex<double> t160 = t1 * t119;
	t172 = t1 * t113;
	complex<double> t184 = t1 * t106;
	complex<double> t194 = -t54 / 0.72e2 - t55 / 0.72e2 + t56 / 0.72e2 - t155 * t13 * dd[0] * t6 / 0.18e2 - t160 * bb[1] * t7 / 0.18e2 - t75 * t2 * dd[2] * t6 / 0.18e2 + t75 * t19 * dd[2] * t6 / 0.18e2 + t172 * bb[0] * t17 / 0.18e2 + t144 * t13 * dd[1] * t6 / 0.18e2 - t144 * t19 * dd[1] * t6 / 0.18e2 - t184 * bb[0] * t11 / 0.18e2 + t160 * bb[2] * t49 / 0.18e2 - t172 * bb[2] * t37 / 0.18e2;
	complex<double> t200 = t105 * t119;
	complex<double> t201 = t83 * dd[1];
	complex<double> t203 = t90 * dd[2];
	complex<double> t206 = bb[2] * cc[1] * dd[0];
	complex<double> t208 = t101 * dd[0];
	complex<double> t210 = t105 * t113;
	complex<double> t213 = bb[2] * cc[0] * dd[1];
	complex<double> t217 = cc[0] * bb[1] * dd[2];
	complex<double> t219 = t105 * t106;
	complex<double> t224 = t155 * t2 * dd[0] * t6 + t184 * bb[1] * t42 + t200 * t201 - t200 * t203 + t200 * t206 - t200 * t208 - t210 * t208 - t210 * t213 + t210 * t201 + t210 * t217 + t219 * t217 + t219 * t206 - t219 * t213 - t219 * t203;
	complex<double> t244 = bb[1] * bb[0] * dd[2];
	complex<double> t255 = t95 * dd[0];
	t260 = t1 * cc[1];
	complex<double> t265 = t76 * dd[1];
	complex<double> t274 = t143 * cc[0] * t13 * dd[1] / 0.6e1 - t3 * t203 / 0.12e2 - t3 * t213 / 0.12e2 + t14 * t217 / 0.12e2 + t74 * cc[1] * t19 * dd[2] / 0.6e1 - t14 * t208 / 0.12e2 - t160 * t244 / 0.6e1 + t39 * t206 / 0.12e2 - t14 * t203 / 0.12e2 + t39 * t217 / 0.12e2 + t14 * t201 / 0.12e2 + t184 * t255 / 0.6e1 + t172 * t244 / 0.6e1 - t260 * cc[0] * t2 * dd[2] / 0.6e1 - t184 * t265 / 0.6e1 + t3 * t201 / 0.12e2 - t39 * t208 / 0.12e2 - t39 * t213 / 0.12e2;
	complex<double> t320 = -t143 * cc[1] * t13 * dd[0] / 0.6e1 - t172 * t255 / 0.6e1 + t260 * cc[2] * t2 * dd[0] / 0.6e1 - t74 * cc[2] * t19 * dd[1] / 0.6e1 + t160 * t265 / 0.6e1 + t3 * t206 / 0.12e2 - t31 * t32 / 0.12e2 - t26 * t27 / 0.12e2 + t21 * t22 / 0.12e2 + t31 * t63 / 0.12e2 - t21 * t60 / 0.12e2 + t26 * t57 / 0.12e2 + t89 * t203 / 0.6e1 - t82 * t201 / 0.6e1 - t144 * t95 * dd[2] / 0.6e1 + t75 * t95 * dd[1] / 0.6e1 + t100 * t208 / 0.6e1 - t75 * t76 * dd[0] / 0.6e1;
	complex<double> t321 = t274 + t320;
	//
	cm[0] = t66 / 0.24e2;
	cm[1] = t66 / 0.12e2;
	cm[2] = t66 / 0.24e2;
	cm[3] = t112 + t150 + t194 + t224 / 0.18e2;
	cm[4] = t321;
	cm[5] = t321;
	//
	t3 = 0.1e1 /  ko /  j;
	t5 = 2.0 * t3 * cm[1];
	t7 = 2.0 * t3 * cm[0];
	t8 =  ko *  ko;
	complex<double> t10 =  j *  j;
	t12 = 0.1e1 /  t8 /  t10;
	t14 = 8.0 * t12 * cm[3];
	t16 = 3.0 * t12 * cm[5];
	t18 = 3.0 * t12 * cm[4];
	complex<double> t20 = 4.0 * t3 * cm[3];
	t24 = 2.0 * t3 * cm[2];
	t31 = 8.0 / t8 / ko / t10 / j * cm[3];
	//
	coefm[0] = t5;
	coefm[1] = t7;
	coefm[2] = t14;
	coefm[3] = -t16;
	coefm[4] = t18;
	coefm[5] = -t20;
	coefm[6] = t20;
	coefm[7] = -t14;
	coefm[8] = cm[3];
	coefm[9] = cm[2];
	coefm[10] = -cm[3];
	coefm[11] = t16;
	coefm[12] = 3.0 * t3 * cm[4];
	coefm[13] = cm[1];
	coefm[14] = cm[0];
	coefm[15] = cm[4];
	coefm[16] = -t24;
	coefm[17] = t31;
	coefm[18] = 3.0 * t3 * cm[5];
	coefm[19] = cm[5];
	coefm[20] = -t31;
	coefm[21] = t24;
	coefm[22] = t31;
	coefm[23] = -t31;
	coefm[24] = -t18;
	coefm[25] = -t5;
	coefm[26] = -t7;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g2_f2
// ***********************************************************************

void coefficients_g2_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[4], cm[4];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = pow(bb[0], 0.2e1);
	complex<double> t3 = t1 * t2;
	complex<double> t4 = t3 * cc[0];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t7 = bb[1] * dd[2] * t6;
	complex<double> t8 = t4 * t7;
	complex<double> t9 = pow(bb[2], 0.2e1);
	complex<double> t10 = t1 * t9;
	complex<double> t11 = t10 * cc[2];
	complex<double> t13 = bb[0] * dd[1] * t6;
	complex<double> t14 = t11 * t13;
	complex<double> t15 = t10 * cc[0];
	complex<double> t16 = t15 * t7;
	complex<double> t17 = pow(bb[1], 0.2e1);
	complex<double> t18 = t1 * t17;
	complex<double> t19 = t18 * cc[0];
	complex<double> t21 = bb[2] * dd[1] * t6;
	complex<double> t22 = t19 * t21;
	complex<double> t23 = t18 * cc[1];
	complex<double> t25 = bb[0] * dd[2] * t6;
	complex<double> t26 = t23 * t25;
	complex<double> t27 = t3 * cc[1];
	complex<double> t29 = bb[2] * dd[0] * t6;
	complex<double> t30 = t27 * t29;
	complex<double> t31 = t10 * cc[1];
	complex<double> t32 = t31 * t25;
	complex<double> t34 = bb[1] * dd[0] * t6;
	complex<double> t35 = t11 * t34;
	complex<double> t36 = t3 * cc[2];
	complex<double> t37 = t36 * t34;
	complex<double> t38 = t18 * cc[2];
	complex<double> t39 = t38 * t13;
	complex<double> t40 = t4 * t21;
	complex<double> t42 = t1 * t17 * bb[1];
	complex<double> t43 = cc[0] * dd[2];
	complex<double> t44 = t43 * t6;
	complex<double> t45 = t42 * t44;
	complex<double> t46 = t23 * t29;
	complex<double> t48 = t1 * t9 * bb[2];
	complex<double> t49 = cc[0] * dd[1];
	complex<double> t50 = t49 * t6;
	complex<double> t51 = t48 * t50;
	complex<double> t53 = t1 * t2 * bb[0];
	complex<double> t54 = cc[2] * dd[1];
	complex<double> t55 = t54 * t6;
	complex<double> t56 = t53 * t55;
	complex<double> t57 = cc[2] * dd[0];
	complex<double> t58 = t57 * t6;
	complex<double> t59 = t42 * t58;
	complex<double> t60 = cc[1] * dd[0];
	complex<double> t61 = t60 * t6;
	complex<double> t62 = t48 * t61;
	complex<double> t63 = cc[1] * dd[2];
	complex<double> t64 = t63 * t6;
	complex<double> t65 = t53 * t64;
	complex<double> t66 = t8 + t14 + t16 - t22 - t26 + t30 - t32 - t35 - t37 + t39 - t40 + t45 + t46 - t51 + t56 - t59 + t62 - t65;
	complex<double> t68 = t1 * cc[1] * cc[0];
	complex<double> t69 = bb[1] * bb[2];
	complex<double> t70 = dd[1] * t6;
	complex<double> t74 = t1 * bb[1];
	complex<double> t75 = t74 * cc[1];
	complex<double> t76 = cc[2] * bb[0];
	complex<double> t81 = t1 * cc[2] * cc[0];
	complex<double> t82 = dd[2] * t6;
	complex<double> t86 = bb[2] * bb[0];
	complex<double> t87 = dd[0] * t6;
	complex<double> t91 = t1 * bb[2];
	complex<double> t92 = t91 * cc[2];
	complex<double> t93 = cc[1] * bb[0];
	complex<double> t97 = t1 * bb[0];
	complex<double> t98 = t97 * cc[0];
	complex<double> t99 = cc[2] * bb[1];
	complex<double> t110 = t68 * t69 * t70 / 0.18e2 - t75 * t76 * t70 / 0.18e2 - t81 * t69 * t82 / 0.18e2 - t68 * t86 * t87 / 0.18e2 + t92 * t93 * t82 / 0.18e2 + t98 * t99 * t87 / 0.18e2 + t8 / 0.72e2 + t14 / 0.72e2 + t16 / 0.72e2 - t22 / 0.72e2 - t26 / 0.72e2 + t30 / 0.72e2 - t32 / 0.72e2;
	complex<double> t120 = pow(cc[2], 0.2e1);
	complex<double> t121 = t91 * t120;
	complex<double> t124 = pow(cc[0], 0.2e1);
	complex<double> t125 = t97 * t124;
	complex<double> t132 = pow(cc[1], 0.2e1);
	complex<double> t133 = t74 * t132;
	complex<double> t140 = -t35 / 0.72e2 - t37 / 0.72e2 + t39 / 0.72e2 - t40 / 0.72e2 + t46 / 0.72e2 + t15 * t55 / 0.18e2 - t19 * t64 / 0.18e2 - t121 * t13 / 0.18e2 - t125 * t7 / 0.18e2 + t27 * t44 / 0.18e2 - t31 * t58 / 0.18e2 + t133 * t25 / 0.18e2 + t121 * t34 / 0.18e2 + t38 * t61 / 0.18e2;
	complex<double> t145 = t6 * t1;
	complex<double> t146 = t145 * t132;
	complex<double> t147 = t99 * dd[0];
	complex<double> t149 = t145 * t120;
	complex<double> t151 = cc[0] * bb[2] * dd[1];
	complex<double> t154 = bb[1] * cc[0] * dd[2];
	complex<double> t157 = t93 * dd[2];
	complex<double> t159 = t145 * t124;
	complex<double> t161 = bb[2] * cc[1] * dd[0];
	complex<double> t165 = t76 * dd[1];
	complex<double> t168 = -t36 * t50 + t125 * t21 - t133 * t29 - t146 * t147 - t149 * t151 + t149 * t154 + t146 * t154 - t149 * t157 + t159 * t161 + t149 * t161 - t146 * t151 + t159 * t165 - t159 * t147;
	complex<double> t179 = t132 * cc[1];
	complex<double> t184 = t124 * cc[0];
	complex<double> t189 = t120 * cc[2];
	complex<double> t206 = -t159 * t157 / 0.18e2 + t146 * t165 / 0.18e2 + t45 / 0.72e2 - t51 / 0.72e2 + t56 / 0.72e2 - t59 / 0.72e2 + t62 / 0.72e2 - t65 / 0.72e2 - t145 * t179 * bb[0] * dd[2] / 0.18e2 + t145 * t184 * bb[1] * dd[2] / 0.18e2 + t145 * t189 * bb[0] * dd[1] / 0.18e2 - t145 * t189 * bb[1] * dd[0] / 0.18e2 + t145 * t179 * bb[2] * dd[0] / 0.18e2 - t145 * t184 * bb[2] * dd[1] / 0.18e2;
	complex<double> t249 = cc[0] * cc[2] * dd[1];
	complex<double> t252 = -t92 * t157 / 0.6e1 + t75 * t165 / 0.6e1 + t81 * t69 * dd[2] / 0.6e1 - t68 * t69 * dd[1] / 0.6e1 - t98 * t147 / 0.6e1 + t68 * t86 * dd[0] / 0.6e1 - t42 * t43 / 0.12e2 + t48 * t49 / 0.12e2 + t53 * t63 / 0.12e2 - t48 * t60 / 0.12e2 + t42 * t57 / 0.12e2 - t53 * t54 / 0.12e2 - t10 * t165 / 0.12e2 - t10 * t154 / 0.12e2 + t18 * t151 / 0.12e2 + t18 * t157 / 0.12e2 + t91 * t120 * bb[0] * dd[1] / 0.6e1 - t10 * t249 / 0.6e1;
	complex<double> t254 = cc[0] * cc[1] * dd[2];
	complex<double> t276 = cc[1] * cc[2] * dd[0];
	complex<double> t303 = t18 * t254 / 0.6e1 - t74 * t132 * bb[0] * dd[2] / 0.6e1 - t3 * t154 / 0.12e2 - t3 * t161 / 0.12e2 + t10 * t157 / 0.12e2 + t10 * t147 / 0.12e2 + t97 * t124 * bb[1] * dd[2] / 0.6e1 - t3 * t254 / 0.6e1 + t10 * t276 / 0.6e1 - t91 * t120 * bb[1] * dd[0] / 0.6e1 - t18 * t161 / 0.12e2 - t18 * t165 / 0.12e2 + t3 * t147 / 0.12e2 + t3 * t151 / 0.12e2 + t74 * t132 * bb[2] * dd[0] / 0.6e1 - t18 * t276 / 0.6e1 + t3 * t249 / 0.6e1 - t97 * t124 * bb[2] * dd[1] / 0.6e1;
	//
	c[0] = t66 / 0.24e2;
	c[1] = t110 + t140 + t168 / 0.18e2 + t206;
	c[2] = -t66 / 0.24e2;
	c[3] = t252 + t303;
	//
	t1 =  ko * ko;
	t3 =  j *  j;
	complex<double> t5 = 0.1e1 / t1 / t3;
	t7 = 3.0 * t5 * c[3];
	t14 = 8.0 / t1 / ko / t3 / j * c[1];
	t17 = 1.0 / ko / j;
	t19 = 2.0 * t17 * c[0];
	t21 = 2.0 * t17 * c[2];
	t23 = 8.0 * t5 * c[1];
	t25 = 4.0 * t17 * c[1];
	//
	coef[0] = c[0];
	coef[1] = t7;
	coef[2] = -t14;
	coef[3] = -t19;
	coef[4] = t21;
	coef[5] = t14;
	coef[6] = -t14;
	coef[7] = c[2];
	coef[8] = c[1];
	coef[9] = -c[1];
	coef[10] = -t7;
	coef[11] = t19;
	coef[12] = t23;
	coef[13] = -t23;
	coef[14] = t25;
	coef[15] = -t25;
	coef[16] = -t21;
	coef[17] = t14;
	coef[18] = 3.0 * t17 * c[3];
	coef[19] = c[3];
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[0], 0.2e1);
	t3 = t1 * t2;
	t4 = t3 * cc[2];
	t5 = cc[0] * dd[1];
	t6 = sqrt(0.3e1);
	t7 = t5 * t6;
	t10 = pow(bb[1], 0.2e1);
	t11 = t1 * t10;
	complex<double> t12 = t11 * cc[0];
	t13 = cc[1] * dd[2];
	t14 = t13 * t6;
	t17 = t1 * bb[0];
	t18 = pow(cc[0], 0.2e1);
	t19 = t17 * t18;
	t21 = bb[1] * dd[2] * t6;
	t25 = bb[2] * dd[1] * t6;
	complex<double> t28 = t1 * bb[1];
	t29 = pow(cc[1], 0.2e1);
	t30 = t28 * t29;
	t32 = bb[0] * dd[2] * t6;
	t35 = pow(bb[2], 0.2e1);
	t36 = t1 * t35;
	t37 = t36 * cc[0];
	t38 = cc[2] * dd[1];
	t39 = t38 * t6;
	t42 = t1 * bb[2];
	t43 = pow(cc[2], 0.2e1);
	t44 = t42 * t43;
	t46 = bb[0] * dd[1] * t6;
	t49 = t11 * cc[2];
	t50 = cc[1] * dd[0];
	t51 = t50 * t6;
	t54 = t36 * cc[1];
	t55 = cc[2] * dd[0];
	t56 = t55 * t6;
	t59 = t3 * cc[1];
	t60 = cc[0] * dd[2];
	t61 = t60 * t6;
	t65 = bb[1] * dd[0] * t6;
	t69 = bb[2] * dd[0] * t6;
	complex<double> t72 = t12 * t25;
	t74 = -t4 * t7 / 0.18e2 - t12 * t14 / 0.18e2 - t19 * t21 / 0.18e2 + t19 * t25 / 0.18e2 + t30 * t32 / 0.18e2 + t37 * t39 / 0.18e2 - t44 * t46 / 0.18e2 + t49 * t51 / 0.18e2 - t54 * t56 / 0.18e2 + t59 * t61 / 0.18e2 + t44 * t65 / 0.18e2 - t30 * t69 / 0.18e2 - t72 / 0.72e2;
	t75 = t36 * cc[2];
	t76 = t75 * t46;
	complex<double> t78 = t4 * t65;
	complex<double> t80 = t37 * t21;
	t82 = t59 * t69;
	complex<double> t84 = t75 * t65;
	t86 = t3 * cc[0];
	t87 = t86 * t25;
	complex<double> t89 = t49 * t46;
	t91 = t86 * t21;
	t93 = t11 * cc[1];
	complex<double> t94 = t93 * t32;
	complex<double> t96 = t54 * t32;
	t98 = t93 * t69;
	complex<double> t100 = t6 * t1;
	complex<double> t101 = t100 * t29;
	complex<double> t103 = cc[0] * bb[1] * dd[2];
	complex<double> t106 = t100 * t18;
	complex<double> t107 = cc[2] * bb[1];
	complex<double> t108 = t107 * dd[0];
	complex<double> t111 = cc[2] * bb[0];
	complex<double> t112 = t111 * dd[1];
	complex<double> t115 = t76 / 0.72e2 - t78 / 0.72e2 + t80 / 0.72e2 + t82 / 0.72e2 - t84 / 0.72e2 - t87 / 0.72e2 + t89 / 0.72e2 + t91 / 0.72e2 - t94 / 0.72e2 - t96 / 0.72e2 + t98 / 0.72e2 + t101 * t103 / 0.18e2 - t106 * t108 / 0.18e2 + t106 * t112 / 0.18e2;
	complex<double> t117 = t100 * t43;
	complex<double> t119 = cc[0] * bb[2] * dd[1];
	t124 = cc[1] * bb[2] * dd[0];
	complex<double> t126 = cc[1] * bb[0];
	complex<double> t127 = t126 * dd[2];
	complex<double> t134 = t1 * cc[1] * cc[0];
	complex<double> t135 = bb[1] * bb[2];
	complex<double> t136 = dd[1] * t6;
	t140 = t1 * cc[2] * cc[0];
	complex<double> t141 = dd[2] * t6;
	complex<double> t144 = bb[0] * bb[2];
	t145 = dd[0] * t6;
	complex<double> t148 = t28 * cc[1];
	t151 = -t117 * t119 + t101 * t112 - t101 * t119 + t117 * t124 - t106 * t127 - t117 * t127 + t117 * t103 - t101 * t108 + t106 * t124 + t134 * t135 * t136 - t140 * t135 * t141 - t134 * t144 * t145 - t148 * t111 * t136;
	complex<double> t152 = t42 * cc[2];
	complex<double> t156 = t17 * cc[0];
	t161 = t1 * t10 * bb[1];
	complex<double> t162 = t161 * t61;
	t165 = t1 * t2 * bb[0];
	complex<double> t166 = t165 * t14;
	complex<double> t169 = t1 * t35 * bb[2];
	complex<double> t170 = t169 * t51;
	complex<double> t172 = t165 * t39;
	complex<double> t174 = t169 * t7;
	complex<double> t176 = t161 * t56;
	complex<double> t178 = t43 * cc[2];
	complex<double> t187 = t18 * cc[0];
	complex<double> t192 = t29 * cc[1];
	complex<double> t205 = t152 * t126 * t141 / 0.18e2 + t156 * t107 * t145 / 0.18e2 + t162 / 0.72e2 - t166 / 0.72e2 + t170 / 0.72e2 + t172 / 0.72e2 - t174 / 0.72e2 - t176 / 0.72e2 - t100 * t178 * bb[1] * dd[0] / 0.18e2 + t100 * t178 * bb[0] * dd[1] / 0.18e2 + t100 * t187 * bb[1] * dd[2] / 0.18e2 - t100 * t192 * bb[0] * dd[2] / 0.18e2 + t100 * t192 * bb[2] * dd[0] / 0.18e2 - t100 * t187 * bb[2] * dd[1] / 0.18e2;
	complex<double> t251 = t161 * t60 / 0.12e2 - t165 * t13 / 0.12e2 + t169 * t50 / 0.12e2 - t161 * t55 / 0.12e2 + t165 * t38 / 0.12e2 - t169 * t5 / 0.12e2 + t152 * t127 / 0.6e1 - t148 * t112 / 0.6e1 - t140 * t135 * dd[2] / 0.6e1 + t134 * t135 * dd[1] / 0.6e1 + t156 * t108 / 0.6e1 - t134 * t144 * dd[0] / 0.6e1 + t36 * t103 / 0.12e2 - t11 * t119 / 0.12e2 - t11 * t127 / 0.12e2 + t36 * t112 / 0.12e2 - t42 * t43 * bb[0] * dd[1] / 0.6e1 + t28 * t29 * bb[0] * dd[2] / 0.6e1;
	complex<double> t259 = cc[0] * cc[2] * dd[1];
	complex<double> t263 = cc[0] * cc[1] * dd[2];
	complex<double> t275 = cc[1] * cc[2] * dd[0];
	complex<double> t302 = -t36 * t127 / 0.12e2 + t3 * t103 / 0.12e2 + t3 * t124 / 0.12e2 + t36 * t259 / 0.6e1 - t11 * t263 / 0.6e1 - t3 * t108 / 0.12e2 - t3 * t119 / 0.12e2 - t36 * t108 / 0.12e2 + t3 * t263 / 0.6e1 - t36 * t275 / 0.6e1 + t42 * t43 * bb[1] * dd[0] / 0.6e1 + t11 * t124 / 0.12e2 + t11 * t112 / 0.12e2 + t17 * t18 * bb[2] * dd[1] / 0.6e1 - t28 * t29 * bb[2] * dd[0] / 0.6e1 + t11 * t275 / 0.6e1 - t3 * t259 / 0.6e1 - t17 * t18 * bb[1] * dd[2] / 0.6e1;
	complex<double> t304 = -t174 - t94 + t170 + t89 - t176 + t91 - t72 - t166 + t76 - t84 - t96 + t172 - t87 + t98 + t80 + t162 + t82 - t78;
	//
	cm[0] = t74 + t115 + t151 / 0.18e2 + t205;
	cm[1] = t251 + t302;
	cm[2] = t304 / 0.24e2;
	cm[3] = -t304 / 0.24e2;
	//
	t3 = 0.1e1 /  ko /  j;
	t5 = 2.0 * t3 * cm[3];
	t6 = ko *  ko;
	t9 = j *  j;
	t14 = 8.0 / t6 / ko / t9 / j * cm[0];
	t17 = 1.0 / t6 / t9;
	t19 = 3.0 * t17 * cm[1];
	t21 = 2.0 * t3 * cm[2];
	t23 = 8.0 * t17 * cm[0];
	t25 = 4.0 * t3 * cm[0];
	//
	coefm[0] = -t5;
	coefm[1] = t14;
	coefm[2] = -t19;
	coefm[3] = t21;
	coefm[4] = t23;
	coefm[5] = -t23;
	coefm[6] = t25;
	coefm[7] = -t25;
	coefm[8] = cm[3];
	coefm[9] = cm[0];
	coefm[10] = -cm[0];
	coefm[11] = -t14;
	coefm[12] = 3.0 * t3 * cm[1];
	coefm[13] = cm[1];
	coefm[14] = -t21;
	coefm[15] = t5;
	coefm[16] = t14;
	coefm[17] = -t14;
	coefm[18] = cm[2];
	coefm[19] = t19;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g2_f3
// ***********************************************************************

void coefficients_g2_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double>t1 = 0.1e1 / AreaT;
	complex<double>t2 = t1 * bb[0];
	complex<double>t3 = pow(cc[0], 0.2e1);
	complex<double>t4 = t2 * t3;
	complex<double>t5 = bb[1] * dd[2];
	complex<double>t6 = sqrt(0.3e1);
	complex<double>t7 = t5 * t6;
	complex<double>t8 = t4 * t7;
	complex<double>t10 = bb[2] * dd[1];
	complex<double>t11 = t10 * t6;
	complex<double>t12 = t4 * t11;
	complex<double>t14 = t1 * bb[1];
	complex<double>t15 = pow(cc[1], 0.2e1);
	complex<double>t16 = t14 * t15;
	complex<double>t17 = bb[2] * dd[0];
	complex<double>t18 = t17 * t6;
	complex<double>t19 = t16 * t18;
	complex<double>t21 = pow(bb[0], 0.2e1);
	complex<double>t22 = t1 * t21;
	complex<double>t23 = t22 * cc[1];
	complex<double>t24 = cc[0] * dd[2];
	complex<double>t25 = t24 * t6;
	complex<double>t26 = t23 * t25;
	complex<double>t28 = pow(bb[2], 0.2e1);
	complex<double>t29 = t1 * t28;
	complex<double>t30 = t29 * cc[1];
	complex<double>t31 = cc[2] * dd[0];
	complex<double>t32 = t31 * t6;
	complex<double>t33 = t30 * t32;
	complex<double>t35 = t22 * cc[2];
	complex<double>t36 = cc[0] * dd[1];
	complex<double>t37 = t36 * t6;
	complex<double>t38 = t35 * t37;
	complex<double>t40 = t29 * cc[0];
	complex<double>t41 = cc[2] * dd[1];
	complex<double>t42 = t41 * t6;
	complex<double>t43 = t40 * t42;
	complex<double>t45 = t1 * bb[2];
	complex<double>t46 = pow(cc[2], 0.2e1);
	complex<double>t47 = t45 * t46;
	complex<double>t48 = bb[0] * dd[1];
	complex<double>t49 = t48 * t6;
	complex<double>t50 = t47 * t49;
	complex<double>t52 = pow(bb[1], 0.2e1);
	complex<double>t53 = t1 * t52;
	complex<double>t54 = t53 * cc[2];
	complex<double>t55 = cc[1] * dd[0];
	complex<double>t56 = t55 * t6;
	complex<double>t57 = t54 * t56;
	complex<double>t59 = bb[1] * dd[0];
	complex<double>t60 = t59 * t6;
	complex<double>t61 = t47 * t60;
	complex<double>t63 = t53 * cc[0];
	complex<double>t64 = cc[1] * dd[2];
	complex<double>t65 = t64 * t6;
	complex<double>t66 = t63 * t65;
	complex<double>t68 = bb[0] * dd[2];
	complex<double>t69 = t68 * t6;
	complex<double>t70 = t16 * t69;
	complex<double>t72 = t22 * cc[0];
	complex<double>t73 = t72 * t7;
	complex<double>t75 = -t8 / 0.18e2 + t12 / 0.18e2 - t19 / 0.18e2 + t26 / 0.18e2 - t33 / 0.18e2 - t38 / 0.18e2 + t43 / 0.18e2 - t50 / 0.18e2 + t57 / 0.18e2 + t61 / 0.18e2 - t66 / 0.18e2 + t70 / 0.18e2 + t73 / 0.72e2;
	complex<double>t76 = t40 * t7;
	complex<double>t78 = t30 * t69;
	complex<double>t80 = t53 * cc[1];
	complex<double>t81 = t80 * t69;
	complex<double>t83 = t80 * t18;
	complex<double>t85 = t23 * t18;
	complex<double>t87 = t54 * t49;
	complex<double>t89 = t29 * cc[2];
	complex<double>t90 = t89 * t49;
	complex<double>t92 = t35 * t60;
	complex<double>t94 = t63 * t11;
	complex<double>t96 = t72 * t11;
	complex<double>t98 = t89 * t60;
	complex<double>t100 = t2 * cc[0];
	complex<double>t101 = cc[2] * bb[1];
	complex<double>t102 = dd[0] * t6;
	complex<double>t104 = t100 * t101 * t102;
	complex<double>t107 = t1 * cc[0] * cc[1];
	complex<double>t108 = bb[2] * bb[0];
	complex<double>t110 = t107 * t108 * t102;
	complex<double>t112 = t45 * cc[2];
	complex<double>t113 = cc[1] * bb[0];
	complex<double>t114 = dd[2] * t6;
	complex<double>t116 = t112 * t113 * t114;
	complex<double>t118 = t76 / 0.72e2 - t78 / 0.72e2 - t81 / 0.72e2 + t83 / 0.72e2 + t85 / 0.72e2 + t87 / 0.72e2 + t90 / 0.72e2 - t92 / 0.72e2 - t94 / 0.72e2 - t96 / 0.72e2 - t98 / 0.72e2 + t104 / 0.18e2 - t110 / 0.18e2 + t116 / 0.18e2;
	complex<double>t120 = t14 * cc[1];
	complex<double>t121 = cc[2] * bb[0];
	complex<double>t122 = dd[1] * t6;
	complex<double>t124 = t120 * t121 * t122;
	complex<double>t126 = bb[1] * bb[2];
	complex<double>t128 = t107 * t126 * t122;
	complex<double>t130 = cc[0] * bb[1];
	complex<double>t132 = t112 * t130 * t114;
	complex<double>t134 = t6 * t1;
	complex<double>t135 = t46 * cc[2];
	complex<double>t144 = t15 * cc[1];
	complex<double>t153 = t3 * cc[0];
	complex<double>t163 = t1 * t21 * bb[0];
	complex<double>t164 = t163 * t65;
	complex<double>t167 = t1 * t28 * bb[2];
	complex<double>t168 = t167 * t56;
	complex<double> t171 = t1 * t52 * bb[1];
	complex<double> t172 = t171 * t25;
	complex<double> t174 = t167 * t37;
	complex<double> t176 = -t124 / 0.18e2 + t128 / 0.18e2 - t132 / 0.18e2 + t134 * t135 * bb[0] * dd[1] / 0.18e2 - t134 * t135 * bb[1] * dd[0] / 0.18e2 - t134 * t144 * bb[0] * dd[2] / 0.18e2 + t134 * t144 * bb[2] * dd[0] / 0.18e2 - t134 * t153 * bb[2] * dd[1] / 0.18e2 + t134 * t153 * bb[1] * dd[2] / 0.18e2 - t164 / 0.72e2 + t168 / 0.72e2 + t172 / 0.72e2 - t174 / 0.72e2;
	complex<double> t177 = t163 * t42;
	complex<double> t179 = t171 * t32;
	complex<double> t181 = t134 * t46;
	complex<double> t183 = cc[0] * bb[2] * dd[1];
	complex<double> t186 = t134 * t3;
	complex<double> t187 = t101 * dd[0];
	complex<double> t190 = t134 * t15;
	complex<double> t191 = t130 * dd[2];
	complex<double> t194 = t113 * dd[2];
	complex<double> t197 = t121 * dd[1];
	complex<double> t201 = cc[1] * bb[2] * dd[0];
	complex<double> t216 = t177 / 0.72e2 - t179 / 0.72e2 - t181 * t183 / 0.18e2 - t186 * t187 / 0.18e2 + t190 * t191 / 0.18e2 - t186 * t194 / 0.18e2 + t190 * t197 / 0.18e2 + t181 * t201 / 0.18e2 + t181 * t191 / 0.18e2 - t190 * t187 / 0.18e2 - t181 * t194 / 0.18e2 + t186 * t201 / 0.18e2 + t186 * t197 / 0.18e2 - t190 * t183 / 0.18e2;
	complex<double> t219 = t98 + t94 - t83 + t164 + t179 + t174 - t76 + t92 - t73 - t85 - t90 + t96 - t87 - t172 + t78 - t168 + t81 - t177;
	complex<double> t220 = t163 * t64;
	complex<double> t221 = t171 * t31;
	complex<double> t222 = t22 * t187;
	complex<double> t223 = t163 * t41;
	complex<double> t224 = t167 * t36;
	complex<double> t225 = t29 * t191;
	complex<double> t226 = t22 * t183;
	complex<double> t227 = t53 * t201;
	complex<double> t228 = t53 * t183;
	complex<double> t229 = t167 * t55;
	complex<double> t230 = t171 * t24;
	complex<double> t231 = t53 * t197;
	complex<double> t232 = t22 * t191;
	complex<double> t233 = t29 * t197;
	complex<double> t234 = t53 * t194;
	complex<double> t235 = t29 * t194;
	complex<double> t236 = t29 * t187;
	complex<double> t237 = t22 * t201;
	complex<double> t238 = t220 + t221 + t222 - t223 + t224 - t225 + t226 - t227 + t228 - t229 - t230 - t231 - t232 - t233 + t234 + t235 + t236 - t237;
	complex<double> t257 = -t8 / 0.12e2 + t12 / 0.12e2 - t19 / 0.12e2 + t26 / 0.12e2 - t33 / 0.12e2 - t38 / 0.12e2 + t43 / 0.12e2 - t50 / 0.12e2 + t57 / 0.12e2 + t61 / 0.12e2 - t66 / 0.12e2 + t70 / 0.12e2 + t73 / 0.24e2 + t76 / 0.24e2 - t78 / 0.24e2 - t81 / 0.24e2 + t83 / 0.24e2 + t85 / 0.24e2;
	complex<double> t276 = t87 / 0.24e2 + t90 / 0.24e2 - t92 / 0.24e2 - t94 / 0.24e2 - t96 / 0.24e2 - t98 / 0.24e2 + t104 / 0.12e2 - t110 / 0.12e2 + t116 / 0.12e2 - t124 / 0.12e2 + t128 / 0.12e2 - t132 / 0.12e2 - t164 / 0.24e2 + t168 / 0.24e2 + t172 / 0.24e2 - t174 / 0.24e2 + t177 / 0.24e2 - t179 / 0.24e2;
	complex<double> t277 = t257 + t276;
	complex<double> t279 = cc[0] * cc[2] * dd[1];
	complex<double> t280 = t29 * t279;
	complex<double> t282 = cc[0] * cc[1] * dd[2];
	complex<double> t283 = t53 * t282;
	complex<double> t286 = t14 * t15 * bb[0] * dd[2];
	complex<double> t289 = t2 * t3 * bb[1] * dd[2];
	complex<double> t290 = t22 * t282;
	complex<double> t292 = cc[1] * cc[2] * dd[0];
	complex<double> t293 = t29 * t292;
	complex<double> t296 = t45 * t46 * bb[1] * dd[0];
	complex<double> t299 = t14 * t15 * bb[2] * dd[0];
	complex<double> t300 = t53 * t292;
	complex<double> t301 = t22 * t279;
	complex<double> t304 = t2 * t3 * bb[2] * dd[1];
	complex<double> t307 = t45 * t46 * bb[0] * dd[1];
	complex<double> t308 = t112 * t194;
	complex<double> t309 = t120 * t197;
	complex<double> t310 = t100 * t187;
	complex<double> t311 = t112 * t191;
	complex<double> t313 = t107 * t126 * dd[1];
	complex<double> t315 = t107 * t108 * dd[0];
	complex<double> t316 = -t280 + t283 - t286 + t289 - t290 + t293 - t296 + t299 - t300 + t301 - t304 + t307 - t308 + t309 - t310 + t311 - t313 + t315;
	complex<double> t325 = t1 * t135;
	complex<double> t328 = t1 * t144;
	complex<double> t331 = t1 * t153;
	complex<double> t341 = t220 / 0.24e2 + t221 / 0.24e2 - t223 / 0.24e2 + t224 / 0.24e2 - t229 / 0.24e2 - t230 / 0.24e2 - t325 * t48 / 0.6e1 + t328 * t68 / 0.6e1 - t331 * t5 / 0.6e1 + t325 * t59 / 0.6e1 - t328 * t17 / 0.6e1 + t331 * t10 / 0.6e1 + t222 / 0.24e2;
	complex<double> t356 = -t225 / 0.24e2 + t226 / 0.24e2 - t227 / 0.24e2 + t228 / 0.24e2 - t231 / 0.24e2 - t232 / 0.24e2 - t233 / 0.24e2 + t234 / 0.24e2 + t235 / 0.24e2 + t236 / 0.24e2 - t237 / 0.24e2 - t280 / 0.6e1 + t283 / 0.6e1 - t286 / 0.6e1;
	complex<double> t358 = t1 * t46;
	complex<double> t361 = t1 * t15;
	complex<double> t364 = t289 - t290 + t293 - t296 + t299 - t300 + t301 - t304 + t307 + t358 * t194 + t358 * t183 - t361 * t191 - t361 * t197;
	complex<double> t365 = t1 * t3;
	complex<double> t374 = t365 * t187 + t365 * t194 - t358 * t201 - t358 * t191 + t361 * t183 + t361 * t187 - t365 * t197 - t365 * t201 - t308 + t309 - t310 + t311 - t313 + t315;
	//
	c[0] = t75 + t118 + t176 + t216;
	c[1] = t219 / 0.24e2;
	c[2] = -t219 / 0.24e2;
	c[3] = t238 / 0.8e1;
	c[4] = t277;
	c[5] = t277;
	c[6] = -t277;
	c[7] = -t238 / 0.8e1;
	c[8] = t238 / 0.12e2 + t316 / 0.6e1;
	c[9] = -t238 / 0.24e2 - t316 / 0.12e2;
	c[10] = t341 + t356 + t364 / 0.6e1 + t374 / 0.6e1;
	//
	t3 = 0.1e1 /  ko /  j;
	t5 = 2.0 * t3 * c[7];
	t6 =  ko *  ko;
	t8 =  j *  j;
	t10 = 0.1e1 / t6 /  t8;
	t12 = 3.0 * t10 * c[4];
	t14 = 2.0 * t3 * c[2];
	t16 = 2.0 * t3 * c[3];
	t18 = 2.0 * t3 * c[1];
	complex<double> t20 = 3.0 * t10 * c[9];
	complex<double> t27 = 8.0 / t6 / ko / t8 / j * c[0];
	t29 = 2.0 * t3 * c[6];
	t33 = 3.0 * t10 * c[10];
	t35 = 3.0 * t3 * c[10];
	t37 = 2.0 * t3 * c[5];
	complex<double> t39 = 3.0 * t10 * c[8];
	t41 = 4.0 * t3 * c[0];
	t43 = 8.0 * t10 * c[0];
	//
	coef[0] = c[10];
	coef[1] = -c[10];
	coef[2] = t5;
	coef[3] = -t12;
	coef[4] = -t14;
	coef[5] = -t16;
	coef[6] = -t18;
	coef[7] = -t20;
	coef[8] = -t27;
	coef[9] = t27;
	coef[10] = t27;
	coef[11] = -t29;
	coef[12] = c[8];
	coef[13] = 3.0 * t3 * c[8];
	coef[14] = -t5;
	coef[15] = -t27;
	coef[16] = -t33;
	coef[17] = t33;
	coef[18] = t14;
	coef[19] = c[7];
	coef[20] = c[6];
	coef[21] = -t35;
	coef[22] = -t37;
	coef[23] = t37;
	coef[24] = -t39;
	coef[25] = t18;
	coef[26] = t16;
	coef[27] = t20;
	coef[28] = t12;
	coef[29] = -t41;
	coef[30] = t43;
	coef[31] = -c[0];
	coef[32] = c[0];
	coef[33] = -t43;
	coef[34] = t41;
	coef[35] = t29;
	coef[36] = t35;
	coef[37] = t39;
	coef[38] = c[5];
	coef[39] = 3.0 * t3 * c[4];
	coef[40] = c[4];
	coef[41] = c[1];
	coef[42] = c[9];
	coef[43] = c[2];
	coef[44] = c[3];
	coef[45] = 3.0 * t3 * c[9];
	coef[46] = t33;
	coef[47] = -t33;
	//
	t1 = 0.1e1 / AreaT;
	t2 = t1 * bb[2];
	t3 = t2 * cc[2];
	t4 = cc[1] * bb[0];
	t5 = t4 * dd[2];
	t6 = t3 * t5;
	complex<double> t9 = t1 * cc[1] * cc[2];
	t10 = bb[1] * bb[0];
	t12 = t9 * t10 * dd[1];
	t15 = t1 * cc[0] * cc[2];
	t17 = t15 * t10 * dd[0];
	t19 = t1 * bb[1];
	t20 = t19 * cc[1];
	t21 = cc[0] * bb[2];
	t22 = t21 * dd[1];
	t23 = t20 * t22;
	t25 = t1 * bb[0];
	t26 = t25 * cc[0];
	t27 = cc[1] * bb[2];
	t28 = t27 * dd[0];
	t29 = t26 * t28;
	t31 = bb[2] * bb[1];
	t33 = t15 * t31 * dd[2];
	t35 = pow(bb[2], 0.2e1);
	t36 = t1 * t35;
	t38 = cc[2] * bb[0] * dd[1];
	t39 = t36 * t38;
	t40 = t39 / 0.24e2;
	t41 = pow(cc[2], 0.2e1);
	complex<double> t44 = t2 * t41 * bb[0] * dd[1];
	t47 = cc[0] * cc[2] * dd[1];
	t48 = t36 * t47;
	complex<double> t51 = cc[0] * bb[1] * dd[2];
	t52 = t36 * t51;
	t53 = t52 / 0.24e2;
	t54 = pow(bb[1], 0.2e1);
	t55 = t1 * t54;
	t57 = cc[0] * cc[1] * dd[2];
	complex<double> t58 = t55 * t57;
	t60 = t55 * t22;
	t61 = t60 / 0.24e2;
	complex<double> t62 = pow(cc[1], 0.2e1);
	t65 = t19 * t62 * bb[0] * dd[2];
	complex<double> t67 = t55 * t5;
	t68 = t67 / 0.24e2;
	t69 = t55 * t28;
	t70 = t69 / 0.24e2;
	t73 = t19 * t62 * bb[2] * dd[0];
	t75 = t55 * t38;
	t76 = t75 / 0.24e2;
	t78 = cc[2] * cc[1] * dd[0];
	complex<double> t79 = t55 * t78;
	t81 = t6 / 0.12e2 - t12 / 0.12e2 + t17 / 0.12e2 + t23 / 0.12e2 - t29 / 0.12e2 - t33 / 0.12e2 + t40 - t44 / 0.12e2 + t48 / 0.12e2 + t53 - t58 / 0.12e2 - t61 + t65 / 0.12e2 - t68 + t70 - t73 / 0.12e2 + t76 + t79 / 0.12e2;
	complex<double> t82 = pow(bb[0], 0.2e1);
	t83 = t1 * t82;
	complex<double> t84 = t83 * t47;
	complex<double> t86 = pow(cc[0], 0.2e1);
	t89 = t25 * t86 * bb[2] * dd[1];
	complex<double> t91 = t83 * t51;
	t92 = t91 / 0.24e2;
	t94 = cc[2] * bb[1] * dd[0];
	complex<double> t95 = t36 * t94;
	t96 = t95 / 0.24e2;
	complex<double> t97 = t83 * t28;
	t98 = t97 / 0.24e2;
	t101 = t25 * t86 * bb[1] * dd[2];
	complex<double> t105 = t2 * t41 * bb[1] * dd[0];
	t107 = t36 * t5;
	t108 = t107 / 0.24e2;
	complex<double> t109 = t36 * t78;
	complex<double> t111 = t83 * t22;
	t112 = t111 / 0.24e2;
	t113 = t83 * t94;
	t114 = t113 / 0.24e2;
	complex<double> t115 = t83 * t57;
	t118 = t1 * t35 * bb[2];
	complex<double> t119 = cc[0] * dd[1];
	t120 = t118 * t119;
	t121 = t120 / 0.24e2;
	complex<double> t123 = t1 * t54 * bb[1];
	t124 = cc[0] * dd[2];
	complex<double> t125 = t123 * t124;
	t126 = t125 / 0.24e2;
	t128 = t1 * t82 * bb[0];
	complex<double> t129 = cc[1] * dd[2];
	t130 = t128 * t129;
	complex<double> t131 = t130 / 0.24e2;
	t132 = cc[1] * dd[0];
	complex<double> t133 = t118 * t132;
	t134 = t133 / 0.24e2;
	t135 = cc[2] * dd[0];
	complex<double> t136 = t123 * t135;
	complex<double> t137 = t136 / 0.24e2;
	complex<double> t138 = cc[2] * dd[1];
	complex<double> t139 = t128 * t138;
	complex<double> t140 = t139 / 0.24e2;
	complex<double> t141 = -t84 / 0.12e2 + t89 / 0.12e2 + t92 - t96 + t98 - t101 / 0.12e2 + t105 / 0.12e2 - t108 - t109 / 0.12e2 - t112 - t114 + t115 / 0.12e2 - t121 + t126 - t131 + t134 - t137 + t140;
	complex<double> t143 = t6 / 0.6e1;
	t144 = t12 / 0.6e1;
	complex<double> t145 = t17 / 0.6e1;
	complex<double> t146 = t23 / 0.6e1;
	complex<double> t147 = t29 / 0.6e1;
	complex<double> t148 = t33 / 0.6e1;
	complex<double> t150 = t44 / 0.6e1;
	complex<double> t151 = t48 / 0.6e1;
	t153 = t58 / 0.6e1;
	complex<double> t155 = t65 / 0.6e1;
	complex<double> t158 = t73 / 0.6e1;
	complex<double> t160 = t79 / 0.6e1;
	complex<double> t161 = t143 - t144 + t145 + t146 - t147 - t148 + t39 / 0.12e2 - t150 + t151 + t52 / 0.12e2 - t153 - t60 / 0.12e2 + t155 - t67 / 0.12e2 + t69 / 0.12e2 - t158 + t75 / 0.12e2 + t160;
	complex<double> t162 = t84 / 0.6e1;
	t163 = t89 / 0.6e1;
	t167 = t101 / 0.6e1;
	t168 = t105 / 0.6e1;
	complex<double> t170 = t109 / 0.6e1;
	complex<double> t173 = t115 / 0.6e1;
	complex<double> t180 = -t162 + t163 + t91 / 0.12e2 - t95 / 0.12e2 + t97 / 0.12e2 - t167 + t168 - t107 / 0.12e2 - t170 - t111 / 0.12e2 - t113 / 0.12e2 + t173 - t120 / 0.12e2 + t125 / 0.12e2 - t130 / 0.12e2 + t133 / 0.12e2 - t136 / 0.12e2 + t139 / 0.12e2;
	complex<double> t182 = -t143 + t144 - t145 - t146 + t147 + t148 - t40 + t150 - t151 - t53 + t153 + t61 - t155;
	t183 = t68 - t70 + t158 - t76 - t160 + t162 - t163 - t92 + t96 - t98 + t167 - t168 + t108 + t170;
	complex<double> t185 = t1 * t41;
	t190 = t1 * t62;
	complex<double> t195 = t1 * t86;
	complex<double> t208 = t112 + t114 - t173 + t185 * t22 / 0.6e1 + t185 * t5 / 0.6e1 - t190 * t51 / 0.6e1 - t190 * t38 / 0.6e1 - t195 * t38 / 0.6e1 - t195 * t28 / 0.6e1 + t190 * t22 / 0.6e1 + t195 * t94 / 0.6e1 + t195 * t5 / 0.6e1 - t185 * t28 / 0.6e1;
	complex<double> t213 = t41 * cc[2];
	complex<double> t214 = t1 * t213;
	complex<double> t215 = bb[0] * dd[1];
	complex<double> t218 = t62 * cc[1];
	t219 = t1 * t218;
	t220 = bb[0] * dd[2];
	t223 = t86 * cc[0];
	t224 = t1 * t223;
	t225 = bb[1] * dd[2];
	t228 = bb[1] * dd[0];
	t231 = bb[2] * dd[1];
	t234 = bb[2] * dd[0];
	t237 = -t185 * t51 / 0.6e1 + t190 * t94 / 0.6e1 + t121 - t126 + t131 - t134 + t137 - t140 - t214 * t215 / 0.6e1 + t219 * t220 / 0.6e1 - t224 * t225 / 0.6e1 + t214 * t228 / 0.6e1 + t224 * t231 / 0.6e1 - t219 * t234 / 0.6e1;
	complex<double> t240 = t19 * t62;
	complex<double> t241 = sqrt(0.3e1);
	complex<double> t242 = t234 * t241;
	complex<double> t243 = t240 * t242;
	complex<double> t245 = t36 * cc[1];
	complex<double> t246 = t135 * t241;
	complex<double> t247 = t245 * t246;
	complex<double> t249 = t83 * cc[2];
	complex<double> t250 = t119 * t241;
	complex<double> t251 = t249 * t250;
	complex<double> t253 = t55 * cc[2];
	complex<double> t254 = t132 * t241;
	complex<double> t255 = t253 * t254;
	t257 = t25 * t86;
	complex<double> t258 = t225 * t241;
	complex<double> t259 = t257 * t258;
	complex<double> t261 = t220 * t241;
	complex<double> t262 = t240 * t261;
	complex<double> t264 = t55 * cc[0];
	complex<double> t265 = t129 * t241;
	complex<double> t266 = t264 * t265;
	complex<double> t268 = t2 * t41;
	complex<double> t269 = t215 * t241;
	complex<double> t270 = t268 * t269;
	complex<double> t272 = t231 * t241;
	complex<double> t273 = t257 * t272;
	complex<double> t275 = t36 * cc[0];
	t276 = t138 * t241;
	t277 = t275 * t276;
	t279 = t228 * t241;
	t280 = t268 * t279;
	t282 = t83 * cc[1];
	t283 = t124 * t241;
	complex<double> t284 = t282 * t283;
	t286 = t282 * t242;
	complex<double> t288 = t55 * cc[1];
	t289 = t288 * t261;
	complex<double> t291 = t83 * cc[0];
	t292 = t291 * t272;
	complex<double> t294 = t36 * cc[2];
	complex<double> t295 = t294 * t269;
	complex<double> t297 = t291 * t258;
	t299 = t253 * t269;
	t301 = t243 / 0.12e2 + t247 / 0.12e2 + t251 / 0.12e2 - t255 / 0.12e2 + t259 / 0.12e2 - t262 / 0.12e2 + t266 / 0.12e2 + t270 / 0.12e2 - t273 / 0.12e2 - t277 / 0.12e2 - t280 / 0.12e2 - t284 / 0.12e2 - t286 / 0.24e2 + t289 / 0.24e2 + t292 / 0.24e2 - t295 / 0.24e2 - t297 / 0.24e2 - t299 / 0.24e2;
	complex<double> t302 = t275 * t258;
	t304 = t249 * t279;
	complex<double> t306 = t264 * t272;
	t308 = t288 * t242;
	t310 = t294 * t279;
	complex<double> t312 = t245 * t261;
	complex<double> t314 = dd[2] * t241;
	t316 = t3 * t4 * t314;
	complex<double> t318 = dd[1] * t241;
	complex<double> t320 = t20 * t21 * t318;
	complex<double> t322 = dd[0] * t241;
	complex<double> t324 = t26 * t27 * t322;
	complex<double> t327 = t15 * t31 * t314;
	complex<double> t330 = t9 * t10 * t318;
	complex<double> t333 = t15 * t10 * t322;
	complex<double> t335 = t118 * t254;
	complex<double> t337 = t128 * t265;
	complex<double> t339 = t123 * t283;
	t341 = t118 * t250;
	complex<double> t343 = t128 * t276;
	complex<double> t345 = t123 * t246;
	complex<double> t347 = -t302 / 0.24e2 + t304 / 0.24e2 + t306 / 0.24e2 - t308 / 0.24e2 + t310 / 0.24e2 + t312 / 0.24e2 - t316 / 0.12e2 - t320 / 0.12e2 + t324 / 0.12e2 + t327 / 0.12e2 + t330 / 0.12e2 - t333 / 0.12e2 - t335 / 0.24e2 + t337 / 0.24e2 - t339 / 0.24e2 + t341 / 0.24e2 - t343 / 0.24e2 + t345 / 0.24e2;
	complex<double> t348 = t301 + t347;
	complex<double> t349 = -t133 + t130 - t125 + t120 - t139 + t136 - t91 - t75 + t67 + t95 - t52 + t113 + t111 + t60 - t97 - t39 + t107 - t69;
	complex<double> t363 = -t243 / 0.18e2 - t247 / 0.18e2 - t251 / 0.18e2 + t255 / 0.18e2 - t259 / 0.18e2 + t262 / 0.18e2 - t266 / 0.18e2 - t270 / 0.18e2 + t273 / 0.18e2 + t277 / 0.18e2 + t280 / 0.18e2 + t284 / 0.18e2 + t286 / 0.72e2;
	complex<double> t375 = t241 * t1;
	complex<double> t376 = t375 * t86;
	complex<double> t379 = t375 * t62;
	complex<double> t384 = -t289 / 0.72e2 - t292 / 0.72e2 + t295 / 0.72e2 + t297 / 0.72e2 + t299 / 0.72e2 + t302 / 0.72e2 - t304 / 0.72e2 - t306 / 0.72e2 + t308 / 0.72e2 - t310 / 0.72e2 - t312 / 0.72e2 + t376 * t28 / 0.18e2 + t379 * t38 / 0.18e2 + t376 * t38 / 0.18e2;
	complex<double> t387 = t375 * t41;
	complex<double> t396 = t379 * t51 - t387 * t22 - t379 * t94 - t387 * t5 + t387 * t28 + t387 * t51 - t379 * t22 - t376 * t5 - t376 * t94 + t316 + t320 - t324 - t327;
	complex<double> t429 = -t330 / 0.18e2 + t333 / 0.18e2 + t335 / 0.72e2 - t337 / 0.72e2 + t339 / 0.72e2 - t341 / 0.72e2 + t343 / 0.72e2 - t345 / 0.72e2 + t375 * t223 * bb[1] * dd[2] / 0.18e2 - t375 * t213 * bb[1] * dd[0] / 0.18e2 + t375 * t218 * bb[2] * dd[0] / 0.18e2 - t375 * t223 * bb[2] * dd[1] / 0.18e2 + t375 * t213 * bb[0] * dd[1] / 0.18e2 - t375 * t218 * bb[0] * dd[2] / 0.18e2;
	complex<double> t432 = t308 + t295 + t286 - t304 - t341 - t289 + t297 + t302 + t339 + t343 - t306 - t345 - t292 + t335 + t299 - t337 - t312 - t310;
	//
	cm[0] = t81 + t141;
	cm[1] = t161 + t180;
	cm[2] = t182 + t183 + t208 + t237;
	cm[3] = t348;
	cm[4] = t348;
	cm[5] = t349 / 0.8e1;
	cm[6] = t363 + t384 + t396 / 0.18e2 + t429;
	cm[7] = t432 / 0.24e2;
	cm[8] = t432 / 0.24e2;
	cm[9] = t349 / 0.8e1;
	cm[10] = t348;
	//
	t3 = 0.1e1 /  ko /  j;
	t6 =  ko *  ko;
	t9 = j *  j;
	t14 = 8.0 / t6 / ko / t9 / j * cm[6];
	t16 = 2.0 * t3 * cm[5];
	t19 = 1.0 / t6 / t9;
	t21 = 3.0 * t19 * cm[10];
	t23 = 2.0 * t3 * cm[9];
	t25 = 4.0 * t3 * cm[6];
	t27 = 8.0 * t19 * cm[6];
	t29 = 3.0 * t19 * cm[1];
	t31 = 2.0 * t3 * cm[8];
	t33 = 2.0 * t3 * cm[4];
	t35 = 3.0 * t19 * cm[0];
	t37 = 2.0 * t3 * cm[7];
	t39 = 3.0 * t3 * cm[2];
	t41 = 2.0 * t3 * cm[3];
	t43 = 3.0 * t19 * cm[2];
	//
	coefm[0] = 3.0 * t3 * cm[1];
	coefm[1] = -t14;
	coefm[2] = -t16;
	coefm[3] = t21;
	coefm[4] = t23;
	coefm[5] = -t25;
	coefm[6] = cm[1];
	coefm[7] = t27;
	coefm[8] = t25;
	coefm[9] = -cm[6];
	coefm[10] = cm[5];
	coefm[11] = cm[3];
	coefm[12] = cm[6];
	coefm[13] = -t29;
	coefm[14] = t31;
	coefm[15] = t33;
	coefm[16] = t35;
	coefm[17] = t37;
	coefm[18] = -t27;
	coefm[19] = -t31;
	coefm[20] = -t37;
	coefm[21] = -t14;
	coefm[22] = t39;
	coefm[23] = t14;
	coefm[24] = -t33;
	coefm[25] = -t39;
	coefm[26] = t41;
	coefm[27] = -t43;
	coefm[28] = -t21;
	coefm[29] = t16;
	coefm[30] = cm[2];
	coefm[31] = -cm[2];
	coefm[32] = -t35;
	coefm[33] = -t23;
	coefm[34] = t14;
	coefm[35] = 3.0 * t3 * cm[0];
	coefm[36] = 3.0 * t3 * cm[10];
	coefm[37] = t29;
	coefm[38] = cm[4];
	coefm[39] = cm[7];
	coefm[40] = cm[10];
	coefm[41] = cm[9];
	coefm[42] = cm[8];
	coefm[43] = cm[0];
	coefm[44] = t43;
	coefm[45] = -t41;
	coefm[46] = -t43;
	coefm[47] = t43;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g3_f1
// ***********************************************************************

void coefficients_g3_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[6], cm[6];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	complex<double> j   = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = pow(cc[1], 0.2e1);
	complex<double> t3 = t1 * t2;
	complex<double> t5 = bb[1] * dd[2];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t7 = t5 * t6;
	complex<double> t8 = t3 * bb[0] * t7;
	complex<double> t10 = pow(cc[2], 0.2e1);
	complex<double> t11 = t1 * t10;
	complex<double> t13 = bb[2] * dd[1];
	complex<double> t14 = t13 * t6;
	complex<double> t15 = t11 * bb[0] * t14;
	complex<double> t17 = t1 * cc[2];
	complex<double> t18 = t17 * cc[0];
	complex<double> t19 = pow(bb[2], 0.2e1);
	complex<double> t22 = t18 * t19 * dd[1] * t6;
	complex<double> t24 = t1 * cc[1];
	complex<double> t25 = t24 * cc[0];
	complex<double> t26 = pow(bb[1], 0.2e1);
	complex<double> t29 = t25 * t26 * dd[2] * t6;
	complex<double> t32 = bb[1] * dd[0];
	complex<double> t33 = t32 * t6;
	complex<double> t34 = t3 * bb[2] * t33;
	complex<double> t36 = t24 * cc[2];
	complex<double> t39 = t36 * t26 * dd[0] * t6;
	complex<double> t41 = pow(bb[0], 0.2e1);
	complex<double> t44 = t18 * t41 * dd[1] * t6;
	complex<double> t46 = pow(cc[0], 0.2e1);
	complex<double> t47 = t1 * t46;
	complex<double> t49 = bb[0] * dd[1];
	complex<double> t50 = t49 * t6;
	complex<double> t51 = t47 * bb[2] * t50;
	complex<double> t54 = bb[0] * dd[2];
	complex<double> t55 = t54 * t6;
	complex<double> t56 = t47 * bb[1] * t55;
	complex<double> t60 = t25 * t41 * dd[2] * t6;
	complex<double> t64 = t36 * t19 * dd[0] * t6;
	complex<double> t67 = bb[2] * dd[0];
	complex<double> t68 = t67 * t6;
	complex<double> t69 = t11 * bb[1] * t68;
	complex<double> t71 = t1 * t19;
	complex<double> t72 = t71 * cc[2];
	complex<double> t73 = t72 * t50;
	complex<double> t75 = t8 / 0.18e2 - t15 / 0.18e2 + t22 / 0.18e2 - t29 / 0.18e2 - t34 / 0.18e2 + t39 / 0.18e2 - t44 / 0.18e2 + t51 / 0.18e2 - t56 / 0.18e2 + t60 / 0.18e2 - t64 / 0.18e2 + t69 / 0.18e2 + t73 / 0.72e2;
	complex<double> t77 = t71 * cc[0] * t7;
	complex<double> t79 = t1 * t26;
	complex<double> t81 = t79 * cc[0] * t14;
	complex<double> t83 = t79 * cc[1];
	complex<double> t84 = t83 * t55;
	complex<double> t86 = t1 * t41;
	complex<double> t87 = t86 * cc[0];
	complex<double> t88 = t87 * t7;
	complex<double> t91 = t86 * cc[1] * t68;
	complex<double> t94 = t71 * cc[1] * t55;
	complex<double> t96 = t72 * t33;
	complex<double> t98 = t83 * t68;
	complex<double> t101 = t79 * cc[2] * t50;
	complex<double> t104 = t86 * cc[2] * t33;
	complex<double> t106 = t87 * t14;
	complex<double> t108 = t6 * t1;
	complex<double> t109 = t108 * t2;
	complex<double> t111 = bb[2] * cc[0] * dd[1];
	complex<double> t114 = t108 * t10;
	complex<double> t116 = bb[0] * cc[1] * dd[2];
	complex<double> t120 = cc[2] * bb[1] * dd[0];
	complex<double> t123 = t77 / 0.72e2 - t81 / 0.72e2 - t84 / 0.72e2 + t88 / 0.72e2 + t91 / 0.72e2 - t94 / 0.72e2 - t96 / 0.72e2 + t98 / 0.72e2 + t101 / 0.72e2 - t104 / 0.72e2 - t106 / 0.72e2 - t109 * t111 / 0.18e2 - t114 * t116 / 0.18e2 - t109 * t120 / 0.18e2;
	complex<double> t126 = t108 * t46;
	complex<double> t128 = cc[2] * bb[0] * dd[1];
	complex<double> t132 = bb[1] * cc[0] * dd[2];
	complex<double> t136 = cc[1] * bb[2] * dd[0];
	complex<double> t141 = bb[0] * dd[0];
	complex<double> t142 = bb[1] * t6;
	complex<double> t144 = t18 * t141 * t142;
	complex<double> t145 = bb[2] * t6;
	complex<double> t147 = t25 * t141 * t145;
	complex<double> t148 = bb[1] * dd[1];
	complex<double> t149 = bb[0] * t6;
	complex<double> t151 = t36 * t148 * t149;
	complex<double> t152 = bb[2] * dd[2];
	complex<double> t154 = t18 * t152 * t142;
	complex<double> t155 = -t114 * t111 + t126 * t128 - t126 * t120 + t114 * t132 + t109 * t128 + t114 * t136 - t126 * t116 + t126 * t136 + t109 * t132 + t144 - t147 - t151 - t154;
	complex<double> t157 = t36 * t152 * t149;
	complex<double> t160 = t25 * t148 * t145;
	complex<double> t163 = t1 * t19 * bb[2];
	complex<double> t164 = cc[0] * dd[1];
	complex<double> t166 = t163 * t164 * t6;
	complex<double> t169 = t1 * t26 * bb[1];
	complex<double> t170 = cc[0] * dd[2];
	complex<double> t172 = t169 * t170 * t6;
	complex<double> t175 = t1 * t41 * bb[0];
	complex<double> t176 = cc[1] * dd[2];
	complex<double> t178 = t175 * t176 * t6;
	complex<double> t180 = cc[1] * dd[0];
	complex<double> t182 = t163 * t180 * t6;
	complex<double> t184 = cc[2] * dd[0];
	complex<double> t186 = t169 * t184 * t6;
	complex<double> t188 = cc[2] * dd[1];
	complex<double> t190 = t175 * t188 * t6;
	complex<double> t192 = t2 * cc[1];
	complex<double> t197 = t10 * cc[2];
	complex<double> t210 = t46 * cc[0];
	complex<double> t219 = t157 / 0.18e2 + t160 / 0.18e2 - t166 / 0.72e2 + t172 / 0.72e2 - t178 / 0.72e2 + t182 / 0.72e2 - t186 / 0.72e2 + t190 / 0.72e2 - t108 * t192 * bb[0] * dd[2] / 0.18e2 + t108 * t197 * bb[0] * dd[1] / 0.18e2 + t108 * t192 * bb[2] * dd[0] / 0.18e2 - t108 * t197 * bb[1] * dd[0] / 0.18e2 + t108 * t210 * bb[1] * dd[2] / 0.18e2 - t108 * t210 * bb[2] * dd[1] / 0.18e2;
	complex<double> t222 = t98 - t104 - t106 + t88 - t81 - t186 + t73 + t190 + t172 - t84 - t166 - t94 + t77 + t101 + t182 - t178 - t96 + t91;
	complex<double> t223 = t60 + t69 + t144 - t56 + t8 - t64 + t22 - t29 - t147 - t15 - t154 + t157 + t160 - t151 - t34 + t39 - t44 + t51;
	complex<double> t224 = t1 * t197;
	complex<double> t227 = t1 * t192;
	complex<double> t230 = t1 * t210;
	complex<double> t241 = t49 * bb[2];
	complex<double> t242 = t11 * t241;
	complex<double> t248 = t17 * cc[0] * t19 * dd[1];
	complex<double> t254 = t24 * cc[0] * t26 * dd[2];
	complex<double> t258 = t54 * bb[1];
	complex<double> t259 = t3 * t258;
	complex<double> t263 = t47 * t258;
	complex<double> t267 = t1 * cc[0];
	complex<double> t270 = t267 * cc[1] * t41 * dd[2];
	complex<double> t272 = -t224 * t49 / 0.6e1 + t227 * t54 / 0.6e1 - t230 * t5 / 0.6e1 + t224 * t32 / 0.6e1 - t227 * t67 / 0.6e1 + t230 * t13 / 0.6e1 + t11 * t116 / 0.6e1 + t242 / 0.12e2 + t11 * t111 / 0.6e1 - t248 / 0.12e2 - t3 * t132 / 0.6e1 + t254 / 0.12e2 - t3 * t128 / 0.6e1 - t259 / 0.12e2 + t47 * t120 / 0.6e1 + t263 / 0.12e2 + t47 * t116 / 0.6e1 - t270 / 0.12e2;
	complex<double> t277 = t17 * cc[1] * t19 * dd[0];
	complex<double> t281 = t32 * bb[2];
	complex<double> t282 = t11 * t281;
	complex<double> t286 = t3 * t281;
	complex<double> t292 = t24 * cc[2] * t26 * dd[0];
	complex<double> t298 = t267 * cc[2] * t41 * dd[1];
	complex<double> t302 = t47 * t241;
	complex<double> t305 = t18 * t152 * bb[1];
	complex<double> t308 = t25 * t148 * bb[2];
	complex<double> t311 = t25 * t141 * bb[2];
	complex<double> t314 = t36 * t152 * bb[0];
	complex<double> t317 = t36 * t148 * bb[0];
	complex<double> t320 = t18 * t141 * bb[1];
	complex<double> t322 = -t11 * t136 / 0.6e1 + t277 / 0.12e2 - t11 * t132 / 0.6e1 - t282 / 0.12e2 + t3 * t111 / 0.6e1 + t286 / 0.12e2 + t3 * t120 / 0.6e1 - t292 / 0.12e2 - t47 * t128 / 0.6e1 + t298 / 0.12e2 - t47 * t136 / 0.6e1 - t302 / 0.12e2 + t305 / 0.12e2 - t308 / 0.12e2 + t311 / 0.12e2 - t314 / 0.12e2 + t317 / 0.12e2 - t320 / 0.12e2;
	complex<double> t348 = t163 * t164 / 0.12e2 - t169 * t170 / 0.12e2 + t175 * t176 / 0.12e2 - t163 * t180 / 0.12e2 + t169 * t184 / 0.12e2 - t175 * t188 / 0.12e2 + t242 / 0.6e1 - t248 / 0.6e1 + t254 / 0.6e1 - t259 / 0.6e1 + t263 / 0.6e1 - t270 / 0.6e1 + t277 / 0.6e1 - t282 / 0.6e1 + t286 / 0.6e1 - t292 / 0.6e1 + t298 / 0.6e1 - t302 / 0.6e1;
	complex<double> t379 = -t71 * t128 / 0.12e2 - t71 * t132 / 0.12e2 + t79 * t111 / 0.12e2 + t79 * t116 / 0.12e2 - t86 * t132 / 0.12e2 - t86 * t136 / 0.12e2 + t71 * t116 / 0.12e2 + t71 * t120 / 0.12e2 - t79 * t136 / 0.12e2 - t79 * t128 / 0.12e2 + t86 * t120 / 0.12e2 + t86 * t111 / 0.12e2 + t305 / 0.6e1 - t308 / 0.6e1 + t311 / 0.6e1 - t314 / 0.6e1 + t317 / 0.6e1 - t320 / 0.6e1;
	complex<double> t399 = -t8 / 0.12e2 + t15 / 0.12e2 - t22 / 0.12e2 + t29 / 0.12e2 + t34 / 0.12e2 - t39 / 0.12e2 + t44 / 0.12e2 - t51 / 0.12e2 + t56 / 0.12e2 - t60 / 0.12e2 + t64 / 0.12e2 - t69 / 0.12e2 - t73 / 0.24e2 - t77 / 0.24e2 + t81 / 0.24e2 + t84 / 0.24e2 - t88 / 0.24e2 - t91 / 0.24e2;
	complex<double> t418 = t94 / 0.24e2 + t96 / 0.24e2 - t98 / 0.24e2 - t101 / 0.24e2 + t104 / 0.24e2 + t106 / 0.24e2 - t144 / 0.12e2 + t147 / 0.12e2 + t151 / 0.12e2 + t154 / 0.12e2 - t157 / 0.12e2 - t160 / 0.12e2 + t166 / 0.24e2 - t172 / 0.24e2 + t178 / 0.24e2 - t182 / 0.24e2 + t186 / 0.24e2 - t190 / 0.24e2;
	//
	c[0] = t75 + t123 + t155 / 0.18e2 + t219;
	c[1] = t222 / 0.24e2;
	c[2] = t223 / 0.12e2;
	c[3] = t272 + t322;
	c[4] = t348 + t379;
	c[5] = t399 + t418;
	//
	t1 =  ko *  ko;
	t3 =  j *  j;
	t5 = 0.1e1 /  t1 /  t3;
	t7 = 3.0 * t5 * c[4];
	t10 = 0.1e1 /  ko /  j;
	t14 = 2.0 * t10 * c[1];
	complex<double> t16 = 2.0 * t10 * c[2];
	t18 = 3.0 * t5 * c[3];
	t25 = 8.0 / t1 / ko / t3 / j * c[0];
	complex<double> t27 = 2.0 * t10 * c[5];
	t29 = 4.0 * t10 * c[0];
	complex<double> t31 = 8.0 * t5 * c[0];
	//
	coef[0] = c[2];
	coef[1] = c[3];
	coef[2] = c[1];
	coef[3] = t7;
	coef[4] = 3.0 * t10 * c[3];
	coef[5] = -t14;
	coef[6] = -t16;
	coef[7] = -t18;
	coef[8] = t25;
	coef[9] = t27;
	coef[10] = -t25;
	coef[11] = -t25;
	coef[12] = c[0];
	coef[13] = c[5];
	coef[14] = -c[0];
	coef[15] = t29;
	coef[16] = -t31;
	coef[17] = -t7;
	coef[18] = -t29;
	coef[19] = t14;
	coef[20] = t18;
	coef[21] = t31;
	coef[22] = t16;
	coef[23] = c[4];
	coef[24] = 3.0 * t10 * c[4];
	coef[25] = -t27;
	coef[26] = t25;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[1], 0.2e1);
	t3 = t1 * t2;
	complex<double> t4 = t3 * cc[1];
	t5 = bb[0] * dd[2];
	t6 = sqrt(0.3e1);
	t7 = t5 * t6;
	t8 = t4 * t7;
	t10 = pow(cc[2], 0.2e1);
	t11 = t1 * t10;
	t13 = bb[2] * dd[0];
	t14 = t13 * t6;
	t15 = t11 * bb[1] * t14;
	t18 = bb[2] * dd[1];
	t19 = t18 * t6;
	complex<double> t20 = t3 * cc[0] * t19;
	t22 = pow(bb[2], 0.2e1);
	complex<double> t23 = t1 * t22;
	t24 = t23 * cc[2];
	t25 = bb[0] * dd[1];
	t26 = t25 * t6;
	t27 = t24 * t26;
	complex<double> t30 = bb[1] * dd[2];
	t31 = t30 * t6;
	t32 = t23 * cc[0] * t31;
	t34 = bb[1] * dd[0];
	complex<double> t35 = t34 * t6;
	t36 = t24 * t35;
	complex<double> t38 = pow(bb[0], 0.2e1);
	t39 = t1 * t38;
	t41 = t39 * cc[1] * t14;
	t44 = t23 * cc[1] * t7;
	t46 = t39 * cc[0];
	t47 = t46 * t31;
	t50 = t3 * cc[2] * t26;
	complex<double> t53 = t39 * cc[2] * t35;
	t55 = t46 * t19;
	complex<double> t57 = t4 * t14;
	complex<double> t59 = t1 * cc[2];
	t60 = t59 * cc[0];
	complex<double> t61 = bb[2] * bb[1];
	complex<double> t62 = dd[2] * t6;
	t64 = t60 * t61 * t62;
	complex<double> t66 = t59 * cc[1];
	t67 = bb[2] * bb[0];
	t69 = t66 * t67 * t62;
	t71 = bb[0] * bb[1];
	t72 = dd[0] * t6;
	complex<double> t74 = t60 * t71 * t72;
	complex<double> t76 = dd[1] * t6;
	complex<double> t78 = t66 * t71 * t76;
	complex<double> t80 = t1 * cc[0];
	t81 = t80 * cc[1];
	t83 = t81 * t67 * t72;
	complex<double> t85 = t8 / 0.24e2 - t15 / 0.12e2 + t20 / 0.24e2 - t27 / 0.24e2 - t32 / 0.24e2 + t36 / 0.24e2 - t41 / 0.24e2 + t44 / 0.24e2 - t47 / 0.24e2 - t50 / 0.24e2 + t53 / 0.24e2 + t55 / 0.24e2 - t57 / 0.24e2 + t64 / 0.12e2 - t69 / 0.12e2 - t74 / 0.12e2 + t78 / 0.12e2 + t83 / 0.12e2;
	t87 = t81 * t61 * t76;
	t91 = t60 * t22 * dd[1] * t6;
	t94 = t11 * bb[0] * t19;
	t98 = t81 * t2 * dd[2] * t6;
	complex<double> t100 = pow(cc[1], 0.2e1);
	t101 = t1 * t100;
	complex<double> t103 = t101 * bb[0] * t31;
	complex<double> t107 = t66 * t2 * dd[0] * t6;
	t111 = t60 * t38 * dd[1] * t6;
	complex<double> t113 = pow(cc[0], 0.2e1);
	t114 = t1 * t113;
	t116 = t114 * bb[2] * t26;
	complex<double> t119 = t101 * bb[2] * t35;
	complex<double> t122 = t114 * bb[1] * t7;
	t126 = t81 * t38 * dd[2] * t6;
	complex<double> t130 = t66 * t22 * dd[0] * t6;
	complex<double> t133 = t1 * t2 * bb[1];
	complex<double> t134 = cc[0] * dd[2];
	t136 = t133 * t134 * t6;
	complex<double> t139 = t1 * t22 * bb[2];
	complex<double> t140 = cc[0] * dd[1];
	t142 = t139 * t140 * t6;
	t144 = cc[1] * dd[0];
	complex<double> t146 = t139 * t144 * t6;
	t149 = t1 * t38 * bb[0];
	complex<double> t150 = cc[1] * dd[2];
	t152 = t149 * t150 * t6;
	t154 = cc[2] * dd[1];
	complex<double> t156 = t149 * t154 * t6;
	complex<double> t158 = cc[2] * dd[0];
	t160 = t133 * t158 * t6;
	complex<double> t162 = -t87 / 0.12e2 - t91 / 0.12e2 + t94 / 0.12e2 + t98 / 0.12e2 - t103 / 0.12e2 - t107 / 0.12e2 + t111 / 0.12e2 - t116 / 0.12e2 + t119 / 0.12e2 + t122 / 0.12e2 - t126 / 0.12e2 + t130 / 0.12e2 - t136 / 0.24e2 + t142 / 0.24e2 - t146 / 0.24e2 + t152 / 0.24e2 - t156 / 0.24e2 + t160 / 0.24e2;
	complex<double> t177 = -t8 / 0.72e2 + t15 / 0.18e2 - t20 / 0.72e2 + t27 / 0.72e2 + t32 / 0.72e2 - t36 / 0.72e2 + t41 / 0.72e2 - t44 / 0.72e2 + t47 / 0.72e2 + t50 / 0.72e2 - t53 / 0.72e2 - t55 / 0.72e2 + t57 / 0.72e2;
	t178 = t6 * t1;
	complex<double> t179 = t178 * t10;
	complex<double> t181 = cc[0] * bb[2] * dd[1];
	t184 = cc[1] * bb[2] * dd[0];
	t186 = t178 * t100;
	t188 = cc[0] * bb[1] * dd[2];
	t190 = t178 * t113;
	complex<double> t193 = bb[0] * cc[1] * dd[2];
	complex<double> t198 = cc[2] * bb[1] * dd[0];
	complex<double> t201 = cc[2] * bb[0] * dd[1];
	complex<double> t206 = -t179 * t181 + t179 * t184 + t186 * t188 + t190 * t184 - t179 * t193 - t186 * t181 - t190 * t193 - t186 * t198 + t190 * t201 + t179 * t188 - t190 * t198 + t186 * t201 - t64 + t69;
	complex<double> t208 = t74 - t78 - t83 + t87 + t91 - t94 - t98 + t103 + t107 - t111 + t116 - t119 - t122;
	complex<double> t217 = t113 * cc[0];
	t222 = t100 * cc[1];
	t227 = t10 * cc[2];
	complex<double> t244 = t126 / 0.18e2 - t130 / 0.18e2 + t136 / 0.72e2 - t142 / 0.72e2 + t146 / 0.72e2 - t152 / 0.72e2 + t156 / 0.72e2 - t160 / 0.72e2 + t178 * t217 * bb[1] * dd[2] / 0.18e2 - t178 * t222 * bb[0] * dd[2] / 0.18e2 - t178 * t227 * bb[1] * dd[0] / 0.18e2 - t178 * t217 * bb[2] * dd[1] / 0.18e2 + t178 * t222 * bb[2] * dd[0] / 0.18e2 + t178 * t227 * bb[0] * dd[1] / 0.18e2;
	complex<double> t247 = -t142 + t47 - t8 - t20 + t50 - t160 + t57 + t156 + t146 - t55 + t32 + t41 - t53 - t36 + t136 + t27 - t44 - t152;
	t248 = -t116 + t111 - t107 - t126 - t91 + t130 - t15 + t94 - t74 + t83 + t98 - t69 - t87 + t122 + t78 + t64 - t103 + t119;
	complex<double> t249 = t1 * t227;
	complex<double> t252 = t1 * t222;
	complex<double> t255 = t1 * t217;
	complex<double> t265 = t81 * t61 * dd[1];
	complex<double> t268 = t81 * t67 * dd[0];
	complex<double> t271 = t60 * t61 * dd[2];
	complex<double> t274 = t66 * t71 * dd[1];
	t277 = t60 * t71 * dd[0];
	complex<double> t280 = t66 * t67 * dd[2];
	complex<double> t284 = t59 * cc[0] * t22 * dd[1];
	complex<double> t290 = t67 * dd[1];
	complex<double> t291 = t11 * t290;
	complex<double> t295 = t1 * cc[1];
	t298 = t295 * cc[0] * t2 * dd[2];
	complex<double> t300 = -t249 * t25 / 0.6e1 + t252 * t5 / 0.6e1 - t255 * t30 / 0.6e1 + t249 * t34 / 0.6e1 + t255 * t18 / 0.6e1 - t252 * t13 / 0.6e1 - t265 / 0.12e2 + t268 / 0.12e2 + t271 / 0.12e2 + t274 / 0.12e2 - t277 / 0.12e2 - t280 / 0.12e2 - t284 / 0.12e2 - t101 * t188 / 0.6e1 + t11 * t193 / 0.6e1 + t291 / 0.12e2 + t11 * t181 / 0.6e1 + t298 / 0.12e2;
	complex<double> t303 = t71 * dd[2];
	complex<double> t304 = t101 * t303;
	complex<double> t310 = t80 * cc[1] * t38 * dd[2];
	complex<double> t316 = t59 * cc[1] * t22 * dd[0];
	complex<double> t318 = t114 * t303;
	complex<double> t328 = t80 * cc[2] * t38 * dd[1];
	complex<double> t332 = t114 * t290;
	complex<double> t334 = t61 * dd[0];
	complex<double> t335 = t101 * t334;
	complex<double> t341 = t295 * cc[2] * t2 * dd[0];
	complex<double> t345 = t11 * t334;
	complex<double> t347 = -t101 * t201 / 0.6e1 - t304 / 0.12e2 + t114 * t193 / 0.6e1 - t310 / 0.12e2 - t11 * t184 / 0.6e1 + t316 / 0.12e2 + t318 / 0.12e2 + t114 * t198 / 0.6e1 + t101 * t198 / 0.6e1 - t114 * t201 / 0.6e1 + t328 / 0.12e2 - t114 * t184 / 0.6e1 - t332 / 0.12e2 + t335 / 0.12e2 + t101 * t181 / 0.6e1 - t341 / 0.12e2 - t11 * t188 / 0.6e1 - t345 / 0.12e2;
	complex<double> t373 = -t139 * t140 / 0.12e2 + t133 * t134 / 0.12e2 - t133 * t158 / 0.12e2 + t149 * t154 / 0.12e2 - t149 * t150 / 0.12e2 + t139 * t144 / 0.12e2 + t265 / 0.6e1 - t268 / 0.6e1 - t271 / 0.6e1 - t274 / 0.6e1 + t277 / 0.6e1 + t280 / 0.6e1 + t284 / 0.6e1 - t291 / 0.6e1 - t298 / 0.6e1 + t304 / 0.6e1 + t310 / 0.6e1 - t316 / 0.6e1;
	complex<double> t404 = -t318 / 0.6e1 - t328 / 0.6e1 + t332 / 0.6e1 - t335 / 0.6e1 + t341 / 0.6e1 + t345 / 0.6e1 - t3 * t193 / 0.12e2 - t3 * t181 / 0.12e2 + t23 * t201 / 0.12e2 + t23 * t188 / 0.12e2 - t23 * t198 / 0.12e2 + t39 * t184 / 0.12e2 - t23 * t193 / 0.12e2 + t39 * t188 / 0.12e2 + t3 * t201 / 0.12e2 - t39 * t198 / 0.12e2 - t39 * t181 / 0.12e2 + t3 * t184 / 0.12e2;
	//
	cm[0] = t85 + t162;
	cm[1] = t177 + t206 / 0.18e2 + t208 / 0.18e2 + t244;
	cm[2] = t247 / 0.24e2;
	cm[3] = t248 / 0.12e2;
	cm[4] = t300 + t347;
	cm[5] = t373 + t404;
	//
	t3 = 0.1e1 /  ko /  j;
	t6 =  ko *  ko;
	t8 =  j * j;
	t10 = 0.1e1 /  t6 /  t8;
	complex<double> t12 = 3.0 * t10 * cm[5];
	t16 = 3.0 * t10 * cm[4];
	t18 = 2.0 * t3 * cm[2];
	t20 = 2.0 * t3 * cm[3];
	t27 = 8.0 / t6 / ko / t8 / j * cm[1];
	t29 = 2.0 * t3 * cm[0];
	t31 = 8.0 * t10 * cm[1];
	t33 = 4.0 * t3 * cm[1];
	//
	coefm[0] = cm[5];
	coefm[1] = 3.0 * t3 * cm[5];
	coefm[2] = t12;
	coefm[3] = 3.0 * t3 * cm[4];
	coefm[4] = cm[2];
	coefm[5] = cm[3];
	coefm[6] = cm[4];
	coefm[7] = -t16;
	coefm[8] = -t18;
	coefm[9] = -t20;
	coefm[10] = t27;
	coefm[11] = -t27;
	coefm[12] = t29;
	coefm[13] = -t27;
	coefm[14] = -t29;
	coefm[15] = t27;
	coefm[16] = -t31;
	coefm[17] = t16;
	coefm[18] = -t33;
	coefm[19] = t20;
	coefm[20] = t33;
	coefm[21] = t31;
	coefm[22] = t18;
	coefm[23] = -t12;
	coefm[24] = cm[1];
	coefm[25] = cm[0];
	coefm[26] = -cm[1];
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g3_f2
// ***********************************************************************

void coefficients_g3_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[6], cm[6];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = sqrt(0.3e1);
	complex<double> t2 = 0.1e1 / AreaT;
	complex<double> t3 = t1 * t2;
	complex<double> t4 = pow(bb[2], 0.2e1);
	complex<double> t5 = t3 * t4;
	complex<double> t6 = cc[2] * bb[0];
	complex<double> t7 = t6 * dd[1];
	complex<double> t8 = t5 * t7;
	complex<double> t9 = pow(bb[0], 0.2e1);
	complex<double> t10 = t3 * t9;
	complex<double> t12 = cc[0] * bb[2] * dd[1];
	complex<double> t13 = t10 * t12;
	complex<double> t14 = pow(bb[1], 0.2e1);
	complex<double> t15 = t3 * t14;
	complex<double> t17 = cc[1] * bb[2] * dd[0];
	complex<double> t18 = t15 * t17;
	complex<double> t19 = t15 * t12;
	complex<double> t20 = t14 * bb[1];
	complex<double> t23 = t3 * t20 * cc[0] * dd[2];
	complex<double> t25 = cc[0] * bb[1] * dd[2];
	complex<double> t26 = t5 * t25;
	complex<double> t27 = t9 * bb[0];
	complex<double> t30 = t3 * t27 * cc[1] * dd[2];
	complex<double> t33 = t3 * t27 * cc[2] * dd[1];
	complex<double> t34 = t4 * bb[2];
	complex<double> t37 = t3 * t34 * cc[1] * dd[0];
	complex<double> t38 = t10 * t17;
	complex<double> t39 = cc[1] * bb[0];
	complex<double> t40 = t39 * dd[2];
	complex<double> t41 = t15 * t40;
	complex<double> t42 = t5 * t40;
	complex<double> t45 = t3 * t34 * cc[0] * dd[1];
	complex<double> t46 = t15 * t7;
	complex<double> t49 = t3 * t20 * cc[2] * dd[0];
	complex<double> t50 = cc[2] * bb[1];
	complex<double> t51 = t50 * dd[0];
	complex<double> t52 = t10 * t51;
	complex<double> t53 = t5 * t51;
	complex<double> t54 = t10 * t25;
	complex<double> t55 = t8 - t13 + t18 - t19 + t23 + t26 - t30 + t33 + t37 + t38 - t41 - t42 - t45 + t46 - t49 - t52 - t53 + t54;
	complex<double> t56 = t2 * cc[1];
	complex<double> t57 = t56 * cc[0];
	complex<double> t58 = bb[1] * bb[2];
	complex<double> t59 = dd[1] * t1;
	complex<double> t61 = t57 * t58 * t59;
	complex<double> t64 = t2 * bb[1] * cc[1];
	complex<double> t66 = t64 * t6 * t59;
	complex<double> t68 = t2 * cc[2];
	complex<double> t69 = t68 * cc[0];
	complex<double> t70 = dd[2] * t1;
	complex<double> t72 = t69 * t58 * t70;
	complex<double> t74 = bb[2] * bb[0];
	complex<double> t75 = dd[0] * t1;
	complex<double> t77 = t57 * t74 * t75;
	complex<double> t80 = t2 * bb[2] * cc[2];
	complex<double> t82 = t80 * t39 * t70;
	complex<double> t85 = t2 * bb[0] * cc[0];
	complex<double> t87 = t85 * t50 * t75;
	complex<double> t96 = t61 / 0.18e2 - t66 / 0.18e2 - t72 / 0.18e2 - t77 / 0.18e2 + t82 / 0.18e2 + t87 / 0.18e2 + t8 / 0.72e2 - t13 / 0.72e2 + t18 / 0.72e2 - t19 / 0.72e2 + t26 / 0.72e2 + t38 / 0.72e2 - t41 / 0.72e2;
	complex<double> t102 = pow(cc[0], 0.2e1);
	complex<double> t103 = t2 * t102;
	complex<double> t105 = bb[0] * dd[1];
	complex<double> t107 = t103 * bb[2] * t105 * t1;
	complex<double> t109 = pow(cc[1], 0.2e1);
	complex<double> t110 = t2 * t109;
	complex<double> t112 = bb[1] * dd[0];
	complex<double> t114 = t110 * bb[2] * t112 * t1;
	complex<double> t116 = pow(cc[2], 0.2e1);
	complex<double> t117 = t2 * t116;
	complex<double> t119 = bb[2] * dd[1];
	complex<double> t121 = t117 * bb[0] * t119 * t1;
	complex<double> t125 = t57 * t14 * dd[2] * t1;
	complex<double> t129 = t69 * t4 * dd[1] * t1;
	complex<double> t132 = bb[0] * dd[2];
	complex<double> t134 = t103 * bb[1] * t132 * t1;
	complex<double> t138 = t57 * t9 * dd[2] * t1;
	complex<double> t140 = t68 * cc[1];
	complex<double> t143 = t140 * t4 * dd[0] * t1;
	complex<double> t146 = bb[1] * dd[2];
	complex<double> t148 = t110 * bb[0] * t146 * t1;
	complex<double> t150 = -t42 / 0.72e2 + t46 / 0.72e2 - t52 / 0.72e2 - t53 / 0.72e2 + t54 / 0.72e2 + t107 / 0.18e2 - t114 / 0.18e2 - t121 / 0.18e2 - t125 / 0.18e2 + t129 / 0.18e2 - t134 / 0.18e2 + t138 / 0.18e2 - t143 / 0.18e2 + t148 / 0.18e2;
	complex<double> t153 = bb[2] * dd[0];
	complex<double> t155 = t117 * bb[1] * t153 * t1;
	complex<double> t158 = t140 * t14 * dd[0] * t1;
	complex<double> t161 = t69 * t9 * dd[1] * t1;
	complex<double> t162 = t3 * t102;
	complex<double> t165 = t3 * t116;
	complex<double> t168 = t3 * t109;
	complex<double> t175 = t155 + t158 - t161 + t162 * t17 + t162 * t7 - t165 * t40 + t165 * t25 - t168 * t51 + t168 * t7 - t165 * t12 + t168 * t25 + t165 * t17 - t168 * t12;
	complex<double> t186 = t109 * cc[1];
	complex<double> t191 = t116 * cc[2];
	complex<double> t200 = t102 * cc[0];
	complex<double> t213 = -t162 * t51 / 0.18e2 - t162 * t40 / 0.18e2 + t37 / 0.72e2 - t49 / 0.72e2 - t45 / 0.72e2 + t23 / 0.72e2 - t30 / 0.72e2 + t33 / 0.72e2 - t3 * t186 * bb[0] * dd[2] / 0.18e2 + t3 * t191 * bb[0] * dd[1] / 0.18e2 - t3 * t191 * bb[1] * dd[0] / 0.18e2 - t3 * t200 * bb[2] * dd[1] / 0.18e2 + t3 * t200 * bb[1] * dd[2] / 0.18e2 + t3 * t186 * bb[2] * dd[0] / 0.18e2;
	complex<double> t234 = t61 / 0.12e2 - t66 / 0.12e2 - t72 / 0.12e2 - t77 / 0.12e2 + t82 / 0.12e2 + t87 / 0.12e2 + t8 / 0.24e2 - t13 / 0.24e2 + t18 / 0.24e2 - t19 / 0.24e2 + t26 / 0.24e2 + t38 / 0.24e2 - t41 / 0.24e2 - t42 / 0.24e2 + t46 / 0.24e2 - t52 / 0.24e2 - t53 / 0.24e2 + t54 / 0.24e2;
	complex<double> t253 = t107 / 0.12e2 - t114 / 0.12e2 - t121 / 0.12e2 - t125 / 0.12e2 + t129 / 0.12e2 - t134 / 0.12e2 + t138 / 0.12e2 - t143 / 0.12e2 + t148 / 0.12e2 + t155 / 0.12e2 + t158 / 0.12e2 - t161 / 0.12e2 + t37 / 0.24e2 - t49 / 0.24e2 - t45 / 0.24e2 + t23 / 0.24e2 - t30 / 0.24e2 + t33 / 0.24e2;
	complex<double> t255 = t80 * t40;
	complex<double> t257 = t64 * t7;
	complex<double> t260 = t69 * t58 * dd[2];
	complex<double> t263 = t57 * t58 * dd[1];
	complex<double> t265 = t85 * t51;
	complex<double> t268 = t57 * t74 * dd[0];
	complex<double> t270 = t2 * t34;
	complex<double> t273 = t270 * cc[0] * dd[1] / 0.12e2;
	complex<double> t274 = t2 * t20;
	complex<double> t277 = t274 * cc[0] * dd[2] / 0.12e2;
	complex<double> t278 = t2 * t191;
	complex<double> t281 = t2 * t186;
	complex<double> t284 = t2 * t27;
	complex<double> t287 = t284 * cc[1] * dd[2] / 0.12e2;
	complex<double> t290 = t270 * cc[1] * dd[0] / 0.12e2;
	complex<double> t291 = t2 * t200;
	complex<double> t294 = -t255 / 0.4e1 + t257 / 0.4e1 + t260 / 0.4e1 - t263 / 0.4e1 - t265 / 0.4e1 + t268 / 0.4e1 + t273 - t277 - t278 * t105 / 0.6e1 + t281 * t132 / 0.6e1 + t287 - t290 - t291 * t146 / 0.6e1;
	complex<double> t299 = t274 * cc[2] * dd[0] / 0.12e2;
	complex<double> t302 = t284 * cc[2] * dd[1] / 0.12e2;
	complex<double> t307 = t2 * t4;
	complex<double> t309 = t307 * t7 / 0.12e2;
	complex<double> t311 = t307 * t25 / 0.12e2;
	complex<double> t312 = t2 * t14;
	complex<double> t314 = t312 * t12 / 0.12e2;
	complex<double> t316 = t312 * t40 / 0.12e2;
	complex<double> t325 = t74 * dd[1];
	complex<double> t326 = t117 * t325;
	complex<double> t328 = t278 * t112 / 0.6e1 + t299 - t302 - t281 * t153 / 0.6e1 + t291 * t119 / 0.6e1 - t309 - t311 + t314 + t316 + t117 * t40 / 0.6e1 + t117 * t12 / 0.6e1 - t110 * t25 / 0.6e1 - t110 * t7 / 0.6e1 + t326 / 0.4e1;
	complex<double> t332 = t68 * cc[0] * t4 * dd[1];
	complex<double> t336 = t56 * cc[0] * t14 * dd[2];
	complex<double> t339 = bb[0] * bb[1] * dd[2];
	complex<double> t340 = t110 * t339;
	complex<double> t342 = t2 * t9;
	complex<double> t344 = t342 * t25 / 0.12e2;
	complex<double> t346 = t342 * t17 / 0.12e2;
	complex<double> t348 = t307 * t40 / 0.12e2;
	complex<double> t350 = t307 * t51 / 0.12e2;
	complex<double> t359 = t103 * t339;
	complex<double> t361 = t2 * cc[0];
	complex<double> t364 = t361 * cc[1] * t9 * dd[2];
	complex<double> t366 = -t332 / 0.4e1 + t336 / 0.4e1 - t340 / 0.4e1 - t344 - t346 + t348 + t350 + t103 * t51 / 0.6e1 + t103 * t40 / 0.6e1 - t117 * t17 / 0.6e1 - t117 * t25 / 0.6e1 + t359 / 0.4e1 - t364 / 0.4e1;
	complex<double> t369 = t68 * cc[1] * t4 * dd[0];
	complex<double> t371 = t58 * dd[0];
	complex<double> t372 = t117 * t371;
	complex<double> t375 = t312 * t17 / 0.12e2;
	complex<double> t377 = t312 * t7 / 0.12e2;
	complex<double> t379 = t342 * t51 / 0.12e2;
	complex<double> t381 = t342 * t12 / 0.12e2;
	complex<double> t390 = t110 * t371;
	complex<double> t394 = t56 * cc[2] * t14 * dd[0];
	complex<double> t398 = t361 * cc[2] * t9 * dd[1];
	complex<double> t400 = t103 * t325;
	complex<double> t402 = t369 / 0.4e1 - t372 / 0.4e1 - t375 - t377 + t379 + t381 + t110 * t12 / 0.6e1 + t110 * t51 / 0.6e1 - t103 * t7 / 0.6e1 - t103 * t17 / 0.6e1 + t390 / 0.4e1 - t394 / 0.4e1 + t398 / 0.4e1 - t400 / 0.4e1;
	complex<double> t413 = -t255 / 0.6e1 + t257 / 0.6e1 + t260 / 0.6e1 - t263 / 0.6e1 - t265 / 0.6e1 + t268 / 0.6e1 + t273 - t277 + t287 - t290 + t299 - t302 - t309 - t311 + t314 + t316 + t326 / 0.6e1 - t332 / 0.6e1;
	complex<double> t424 = t336 / 0.6e1 - t340 / 0.6e1 - t344 - t346 + t348 + t350 + t359 / 0.6e1 - t364 / 0.6e1 + t369 / 0.6e1 - t372 / 0.6e1 - t375 - t377 + t379 + t381 + t390 / 0.6e1 - t394 / 0.6e1 + t398 / 0.6e1 - t400 / 0.6e1;
	complex<double> t426 = t61 - t66 - t72 - t77 + t82 + t87 + t8 - t13 + t18 - t19 + t26 + t38 - t41 - t42 + t46 - t52 - t53 + t54;
	complex<double> t427 = t107 - t114 - t121 - t125 + t129 - t134 + t138 - t143 + t148 + t155 + t158 - t161 + t37 - t49 - t45 + t23 - t30 + t33;
	//
	c[0] = t55 / 0.24e2;
	c[1] = t96 + t150 + t175 / 0.18e2 + t213;
	c[2] = t234 + t253;
	c[3] = t294 + t328 + t366 + t402;
	c[4] = t413 + t424;
	c[5] = t426 / 0.12e2 + t427 / 0.12e2;
	//
	t1 = ko * ko;
	t4 =  j *  j;
	t9 = 8.0 / t1 / ko / t4 / j * c[1];
	t12 = 1.0 / ko / j;
	t17 = 1.0 / t1 / t4;
	t19 = 3.0 * t17 * c[4];
	t23 = 2.0 * t12 * c[0];
	t25 = 2.0 * t12 * c[5];
	t27 = 3.0 * t17 * c[3];
	complex<double> t29 = 2.0 * t12 * c[2];
	complex<double> t31 = 4.0 * t12 * c[1];
	t33 = 8.0 * t17 * c[1];
	//
	coef[0] = -t9;
	coef[1] = c[4];
	coef[2] = 3.0 * t12 * c[4];
	coef[3] = c[0];
	coef[4] = c[3];
	coef[5] = c[5];
	coef[6] = t19;
	coef[7] = 3.0 * t12 * c[3];
	coef[8] = -t23;
	coef[9] = -t25;
	coef[10] = -t27;
	coef[11] = t9;
	coef[12] = -t9;
	coef[13] = t29;
	coef[14] = -t29;
	coef[15] = t9;
	coef[16] = c[1];
	coef[17] = c[2];
	coef[18] = -c[1];
	coef[19] = t27;
	coef[20] = t31;
	coef[21] = t25;
	coef[22] = t23;
	coef[23] = -t31;
	coef[24] = t33;
	coef[25] = -t19;
	coef[26] = -t33;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[2], 0.2e1);
	t3 = t1 * t2;
	t5 = cc[2] * bb[0] * dd[1];
	t7 = t3 * t5 / 0.12e2;
	t8 = pow(cc[0], 0.2e1);
	t9 = t1 * t8;
	t10 = bb[2] * bb[0];
	complex<double> t11 = t10 * dd[1];
	t12 = t9 * t11;
	t14 = bb[1] * bb[0];
	t15 = t14 * dd[2];
	complex<double> t16 = t9 * t15;
	t19 = bb[2] * cc[1] * dd[0];
	complex<double> t22 = pow(cc[1], 0.2e1);
	t23 = t1 * t22;
	complex<double> t24 = bb[2] * bb[1];
	t25 = t24 * dd[0];
	t26 = t23 * t25;
	t29 = bb[2] * cc[0] * dd[1];
	complex<double> t32 = t1 * cc[1];
	t33 = pow(bb[1], 0.2e1);
	complex<double> t36 = t32 * cc[2] * t33 * dd[0];
	t38 = t23 * t15;
	t40 = pow(cc[2], 0.2e1);
	t41 = t1 * t40;
	t42 = t41 * t11;
	complex<double> t47 = cc[0] * bb[1] * dd[2];
	t53 = bb[0] * cc[1] * dd[2];
	t56 = -t7 - t12 / 0.4e1 + t16 / 0.4e1 - t9 * t19 / 0.6e1 + t26 / 0.4e1 + t23 * t29 / 0.6e1 - t36 / 0.4e1 - t38 / 0.4e1 + t42 / 0.4e1 - t23 * t5 / 0.6e1 - t23 * t47 / 0.6e1 + t41 * t29 / 0.6e1 + t41 * t53 / 0.6e1;
	t58 = t3 * t47 / 0.12e2;
	t59 = t1 * t33;
	t61 = t59 * t29 / 0.12e2;
	complex<double> t63 = t59 * t53 / 0.12e2;
	t66 = t32 * cc[0] * t33 * dd[2];
	t69 = cc[2] * bb[1] * dd[0];
	t72 = t1 * cc[0];
	complex<double> t73 = pow(bb[0], 0.2e1);
	complex<double> t76 = t72 * cc[2] * t73 * dd[1];
	t80 = t1 * cc[2];
	complex<double> t83 = t80 * cc[0] * t2 * dd[1];
	t85 = t1 * t73;
	t87 = t85 * t47 / 0.12e2;
	complex<double> t89 = t85 * t19 / 0.12e2;
	complex<double> t93 = t3 * t53 / 0.12e2;
	complex<double> t95 = t3 * t69 / 0.12e2;
	complex<double> t98 = -t58 + t61 + t63 + t66 / 0.4e1 + t23 * t69 / 0.6e1 + t76 / 0.4e1 - t9 * t5 / 0.6e1 - t83 / 0.4e1 - t87 - t89 + t9 * t53 / 0.6e1 + t93 + t95 + t9 * t69 / 0.6e1;
	complex<double> t101 = t85 * t29 / 0.12e2;
	t103 = t85 * t69 / 0.12e2;
	t105 = t59 * t19 / 0.12e2;
	t107 = t59 * t5 / 0.12e2;
	t110 = t72 * cc[1] * t73 * dd[2];
	t114 = t80 * cc[1] * t2 * dd[0];
	complex<double> t118 = t41 * t25;
	complex<double> t122 = t8 * cc[0];
	complex<double> t123 = t1 * t122;
	complex<double> t124 = bb[1] * dd[2];
	complex<double> t127 = t40 * cc[2];
	complex<double> t128 = t1 * t127;
	t129 = bb[1] * dd[0];
	t132 = t2 * bb[2];
	complex<double> t133 = t1 * t132;
	complex<double> t136 = t133 * cc[1] * dd[0] / 0.12e2;
	complex<double> t137 = bb[0] * dd[1];
	t140 = t101 + t103 - t105 - t107 - t110 / 0.4e1 + t114 / 0.4e1 - t41 * t19 / 0.6e1 - t118 / 0.4e1 - t41 * t47 / 0.6e1 - t123 * t124 / 0.6e1 + t128 * t129 / 0.6e1 - t136 - t128 * t137 / 0.6e1;
	complex<double> t141 = t22 * cc[1];
	complex<double> t142 = t1 * t141;
	t143 = bb[0] * dd[2];
	t146 = t73 * bb[0];
	complex<double> t147 = t1 * t146;
	t150 = t147 * cc[1] * dd[2] / 0.12e2;
	complex<double> t151 = bb[2] * dd[0];
	complex<double> t154 = bb[2] * dd[1];
	complex<double> t157 = t33 * bb[1];
	t158 = t1 * t157;
	t161 = t158 * cc[2] * dd[0] / 0.12e2;
	complex<double> t164 = t147 * cc[2] * dd[1] / 0.12e2;
	complex<double> t167 = t133 * cc[0] * dd[1] / 0.12e2;
	complex<double> t170 = t158 * cc[0] * dd[2] / 0.12e2;
	complex<double> t171 = t72 * cc[1];
	complex<double> t173 = t171 * t10 * dd[0];
	t175 = t80 * cc[1];
	complex<double> t177 = t175 * t10 * dd[2];
	complex<double> t180 = t175 * t14 * dd[1];
	complex<double> t182 = t72 * cc[2];
	complex<double> t184 = t182 * t14 * dd[0];
	complex<double> t187 = t182 * t24 * dd[2];
	complex<double> t190 = t171 * t24 * dd[1];
	complex<double> t192 = t142 * t143 / 0.6e1 + t150 - t142 * t151 / 0.6e1 + t123 * t154 / 0.6e1 + t161 - t164 + t167 - t170 + t173 / 0.4e1 - t177 / 0.4e1 + t180 / 0.4e1 - t184 / 0.4e1 + t187 / 0.4e1 - t190 / 0.4e1;
	complex<double> t204 = t7 + t12 / 0.6e1 - t16 / 0.6e1 - t26 / 0.6e1 + t36 / 0.6e1 + t38 / 0.6e1 - t42 / 0.6e1 + t58 - t61 - t63 - t66 / 0.6e1 - t76 / 0.6e1 + t83 / 0.6e1 + t87 + t89 - t93 - t95 - t101;
	complex<double> t214 = -t103 + t105 + t107 + t110 / 0.6e1 - t114 / 0.6e1 + t118 / 0.6e1 + t136 - t150 - t161 + t164 - t167 + t170 - t173 / 0.6e1 + t177 / 0.6e1 - t180 / 0.6e1 + t184 / 0.6e1 - t187 / 0.6e1 + t190 / 0.6e1;
	complex<double> t216 = sqrt(0.3e1);
	complex<double> t217 = t216 * t1;
	complex<double> t218 = t217 * t40;
	complex<double> t221 = t217 * t8;
	complex<double> t224 = t217 * t22;
	complex<double> t235 = t217 * t146 * cc[2] * dd[1];
	complex<double> t239 = t217 * t132 * cc[0] * dd[1];
	complex<double> t243 = t217 * t157 * cc[2] * dd[0];
	complex<double> t247 = t217 * t132 * cc[1] * dd[0];
	complex<double> t251 = t217 * t146 * cc[1] * dd[2];
	t255 = t217 * t157 * cc[0] * dd[2];
	complex<double> t261 = -t218 * t29 / 0.18e2 - t221 * t69 / 0.18e2 - t224 * t29 / 0.18e2 - t218 * t53 / 0.18e2 - t224 * t69 / 0.18e2 - t221 * t53 / 0.18e2 + t235 / 0.72e2 - t239 / 0.72e2 - t243 / 0.72e2 + t247 / 0.72e2 - t251 / 0.72e2 + t255 / 0.72e2 - t217 * t127 * bb[1] * dd[0] / 0.18e2;
	t278 = dd[0] * t216;
	complex<double> t280 = t171 * t10 * t278;
	t281 = dd[2] * t216;
	complex<double> t283 = t182 * t24 * t281;
	complex<double> t285 = t175 * t10 * t281;
	complex<double> t286 = dd[1] * t216;
	complex<double> t288 = t171 * t24 * t286;
	t290 = t175 * t14 * t286;
	complex<double> t292 = t182 * t14 * t278;
	complex<double> t295 = t9 * bb[2] * t137 * t216;
	complex<double> t298 = t41 * bb[0] * t154 * t216;
	t299 = -t217 * t141 * bb[0] * dd[2] - t217 * t122 * bb[2] * dd[1] + t217 * t122 * bb[1] * dd[2] + t217 * t127 * bb[0] * dd[1] + t217 * t141 * bb[2] * dd[0] + t218 * t19 - t280 - t283 + t285 + t288 - t290 + t292 + t295 - t298;
	complex<double> t303 = t41 * bb[1] * t151 * t216;
	t307 = t171 * t73 * dd[2] * t216;
	t311 = t182 * t2 * dd[1] * t216;
	complex<double> t315 = t171 * t33 * dd[2] * t216;
	complex<double> t319 = t182 * t73 * dd[1] * t216;
	complex<double> t323 = t23 * bb[0] * t124 * t216;
	complex<double> t327 = t23 * bb[2] * t129 * t216;
	complex<double> t331 = t175 * t33 * dd[0] * t216;
	complex<double> t335 = t9 * bb[1] * t143 * t216;
	t339 = t175 * t2 * dd[0] * t216;
	complex<double> t341 = t217 * t33;
	t342 = t341 * t29;
	t344 = t217 * t2;
	complex<double> t345 = t344 * t53;
	complex<double> t347 = t217 * t73;
	t348 = t347 * t47;
	t350 = t303 / 0.18e2 + t307 / 0.18e2 + t311 / 0.18e2 - t315 / 0.18e2 - t319 / 0.18e2 + t323 / 0.18e2 - t327 / 0.18e2 + t331 / 0.18e2 - t335 / 0.18e2 - t339 / 0.18e2 - t342 / 0.72e2 - t345 / 0.72e2 + t348 / 0.72e2;
	complex<double> t351 = t341 * t53;
	complex<double> t353 = t344 * t5;
	complex<double> t355 = t344 * t69;
	complex<double> t357 = t347 * t19;
	t359 = t344 * t47;
	t361 = t341 * t19;
	complex<double> t363 = t341 * t5;
	complex<double> t365 = t347 * t69;
	complex<double> t367 = t347 * t29;
	t379 = -t351 / 0.72e2 + t353 / 0.72e2 - t355 / 0.72e2 + t357 / 0.72e2 + t359 / 0.72e2 + t361 / 0.72e2 + t363 / 0.72e2 - t365 / 0.72e2 - t367 / 0.72e2 + t224 * t47 / 0.18e2 + t224 * t5 / 0.18e2 + t221 * t19 / 0.18e2 + t221 * t5 / 0.18e2 + t218 * t47 / 0.18e2;
	complex<double> t382 = -t355 + t357 - t351 + t359 - t243 - t251 + t255 - t365 + t363 + t361 - t345 + t353 + t247 - t239 - t367 - t342 + t348 + t235;
	complex<double> t383 = -t235 + t239 + t243 - t247 + t251 - t255 + t280 + t283 - t285 - t288 + t290 - t292 - t295 + t298 - t303 - t307 - t311 + t315;
	complex<double> t384 = t319 - t323 + t327 - t331 + t335 + t339 + t342 + t345 - t348 + t351 - t353 + t355 - t357 - t359 - t361 - t363 + t365 + t367;
	complex<double> t404 = t235 / 0.24e2 - t239 / 0.24e2 - t243 / 0.24e2 + t247 / 0.24e2 - t251 / 0.24e2 + t255 / 0.24e2 - t280 / 0.12e2 - t283 / 0.12e2 + t285 / 0.12e2 + t288 / 0.12e2 - t290 / 0.12e2 + t292 / 0.12e2 + t295 / 0.12e2 - t298 / 0.12e2 + t303 / 0.12e2 + t307 / 0.12e2 + t311 / 0.12e2 - t315 / 0.12e2;
	complex<double> t423 = -t319 / 0.12e2 + t323 / 0.12e2 - t327 / 0.12e2 + t331 / 0.12e2 - t335 / 0.12e2 - t339 / 0.12e2 - t342 / 0.24e2 - t345 / 0.24e2 + t348 / 0.24e2 - t351 / 0.24e2 + t353 / 0.24e2 - t355 / 0.24e2 + t357 / 0.24e2 + t359 / 0.24e2 + t361 / 0.24e2 + t363 / 0.24e2 - t365 / 0.24e2 - t367 / 0.24e2;
	//
	cm[0] = t56 + t98 + t140 + t192;
	cm[1] = t204 + t214;
	cm[2] = t261 + t299 / 0.18e2 + t350 + t379;
	cm[3] = t382 / 0.24e2;
	cm[4] = t383 / 0.12e2 + t384 / 0.12e2;
	cm[5] = t404 + t423;
	//
	t1 =  ko *  ko;
	t4 =  j * j;
	t9 = 8.0 / t1 / ko / t4 / j * cm[2];
	t12 = 1.0 / ko / j;
	t17 = 1.0 / t1 / t4;
	t19 = 3.0 * t17 * cm[1];
	t23 = 3.0 * t17 * cm[0];
	t25 = 2.0 * t12 * cm[4];
	t27 = 2.0 * t12 * cm[3];
	t29 = 2.0 * t12 * cm[5];
	t31 = 4.0 * t12 * cm[2];
	t33 = 8.0 * t17 * cm[2];
	//
	coefm[0] = t9;
	coefm[1] = cm[1];
	coefm[2] = 3.0 * t12 * cm[1];
	coefm[3] = cm[4];
	coefm[4] = cm[3];
	coefm[5] = cm[0];
	coefm[6] = t19;
	coefm[7] = 3.0 * t12 * cm[0];
	coefm[8] = -t23;
	coefm[9] = -t25;
	coefm[10] = -t27;
	coefm[11] = t9;
	coefm[12] = -t9;
	coefm[13] = t29;
	coefm[14] = -t29;
	coefm[15] = -t9;
	coefm[16] = cm[2];
	coefm[17] = cm[5];
	coefm[18] = -cm[2];
	coefm[19] = -t31;
	coefm[20] = t27;
	coefm[21] = -t33;
	coefm[22] = t23;
	coefm[23] = t31;
	coefm[24] = t25;
	coefm[25] = -t19;
	coefm[26] = t33;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_g3_f3
// ***********************************************************************

void coefficients_g3_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, double Ap, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[11], cm[11];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> j = Iunit;
	//
	double AreaT = Ap;
	//
	complex<double> t1 = 0.1e1 / AreaT;
	complex<double> t2 = t1 * bb[0];
	complex<double> t3 = pow(cc[0], 0.2e1);
	complex<double> t4 = t2 * t3;
	complex<double> t5 = bb[2] * dd[1];
	complex<double> t6 = sqrt(0.3e1);
	complex<double> t7 = t5 * t6;
	complex<double> t8 = t4 * t7;
	complex<double> t10 = bb[1] * dd[2];
	complex<double> t11 = t10 * t6;
	complex<double> t12 = t4 * t11;
	complex<double> t14 = pow(bb[2], 0.2e1);
	complex<double> t15 = t1 * t14;
	complex<double> t16 = t15 * cc[1];
	complex<double> t17 = cc[2] * dd[0];
	complex<double> t18 = t17 * t6;
	complex<double> t19 = t16 * t18;
	complex<double> t21 = pow(bb[0], 0.2e1);
	complex<double> t22 = t1 * t21;
	complex<double> t23 = t22 * cc[2];
	complex<double> t24 = cc[0] * dd[1];
	complex<double> t25 = t24 * t6;
	complex<double> t26 = t23 * t25;
	complex<double> t28 = t15 * cc[0];
	complex<double> t29 = cc[2] * dd[1];
	complex<double> t30 = t29 * t6;
	complex<double> t31 = t28 * t30;
	complex<double> t33 = t22 * cc[1];
	complex<double> t34 = cc[0] * dd[2];
	complex<double> t35 = t34 * t6;
	complex<double> t36 = t33 * t35;
	complex<double> t38 = t1 * bb[1];
	complex<double> t39 = pow(cc[1], 0.2e1);
	complex<double> t40 = t38 * t39;
	complex<double> t41 = bb[2] * dd[0];
	complex<double> t42 = t41 * t6;
	complex<double> t43 = t40 * t42;
	complex<double> t45 = t1 * bb[2];
	complex<double> t46 = pow(cc[2], 0.2e1);
	complex<double> t47 = t45 * t46;
	complex<double> t48 = bb[1] * dd[0];
	complex<double> t49 = t48 * t6;
	complex<double> t50 = t47 * t49;
	complex<double> t52 = pow(bb[1], 0.2e1);
	complex<double> t53 = t1 * t52;
	complex<double> t54 = t53 * cc[0];
	complex<double> t55 = cc[1] * dd[2];
	complex<double> t56 = t55 * t6;
	complex<double> t57 = t54 * t56;
	complex<double> t59 = bb[0] * dd[2];
	complex<double> t60 = t59 * t6;
	complex<double> t61 = t40 * t60;
	complex<double> t63 = bb[0] * dd[1];
	complex<double> t64 = t63 * t6;
	complex<double> t65 = t47 * t64;
	complex<double> t67 = t53 * cc[2];
	complex<double> t68 = cc[1] * dd[0];
	complex<double> t69 = t68 * t6;
	complex<double> t70 = t67 * t69;
	complex<double> t72 = t33 * t42;
	complex<double> t74 = t8 / 0.18e2 - t12 / 0.18e2 - t19 / 0.18e2 - t26 / 0.18e2 + t31 / 0.18e2 + t36 / 0.18e2 - t43 / 0.18e2 + t50 / 0.18e2 - t57 / 0.18e2 + t61 / 0.18e2 - t65 / 0.18e2 + t70 / 0.18e2 + t72 / 0.72e2;
	complex<double> t75 = t22 * cc[0];
	complex<double> t76 = t75 * t7;
	complex<double> t78 = t16 * t60;
	complex<double> t80 = t53 * cc[1];
	complex<double> t81 = t80 * t60;
	complex<double> t83 = t23 * t49;
	complex<double> t85 = t75 * t11;
	complex<double> t87 = t15 * cc[2];
	complex<double> t88 = t87 * t64;
	complex<double> t90 = t54 * t7;
	complex<double> t92 = t28 * t11;
	complex<double> t94 = t80 * t42;
	complex<double> t96 = t87 * t49;
	complex<double> t98 = t67 * t64;
	complex<double> t100 = t6 * t1;
	complex<double> t101 = t100 * t46;
	complex<double> t103 = cc[1] * bb[2] * dd[0];
	complex<double> t104 = t101 * t103;
	complex<double> t107 = bb[1] * cc[0] * dd[2];
	complex<double> t108 = t101 * t107;
	complex<double> t110 = t100 * t3;
	complex<double> t111 = cc[2] * bb[0];
	complex<double> t112 = t111 * dd[1];
	complex<double> t113 = t110 * t112;
	complex<double> t115 = -t76 / 0.72e2 - t78 / 0.72e2 - t81 / 0.72e2 - t83 / 0.72e2 + t85 / 0.72e2 + t88 / 0.72e2 - t90 / 0.72e2 + t92 / 0.72e2 + t94 / 0.72e2 - t96 / 0.72e2 + t98 / 0.72e2 + t104 / 0.18e2 + t108 / 0.18e2 + t113 / 0.18e2;
	complex<double> t118 = cc[0] * bb[2] * dd[1];
	complex<double> t121 = cc[2] * bb[1];
	complex<double> t122 = t121 * dd[0];
	complex<double> t124 = t100 * t39;
	complex<double> t127 = cc[1] * bb[0] * dd[2];
	complex<double> t134 = t1 * cc[0] * cc[1];
	complex<double> t135 = bb[2] * bb[0];
	complex<double> t136 = dd[0] * t6;
	complex<double> t138 = t134 * t135 * t136;
	complex<double> t139 = t2 * cc[0];
	complex<double> t141 = t139 * t121 * t136;
	complex<double> t142 = t38 * cc[1];
	complex<double> t143 = dd[1] * t6;
	complex<double> t145 = t142 * t111 * t143;
	complex<double> t146 = bb[1] * bb[2];
	complex<double> t148 = t134 * t146 * t143;
	complex<double> t149 = -t101 * t118 + t110 * t103 - t110 * t122 + t124 * t112 - t110 * t127 - t124 * t118 - t124 * t122 - t101 * t127 + t124 * t107 - t138 + t141 - t145 + t148;
	complex<double> t150 = t1 * cc[2];
	complex<double> t151 = t150 * cc[0];
	complex<double> t152 = dd[2] * t6;
	complex<double> t154 = t151 * t146 * t152;
	complex<double> t156 = t150 * cc[1];
	complex<double> t158 = t156 * t135 * t152;
	complex<double> t161 = t1 * t52 * bb[1];
	complex<double> t162 = t161 * t35;
	complex<double> t165 = t1 * t14 * bb[2];
	complex<double> t166 = t165 * t25;
	complex<double> t168 = t161 * t18;
	complex<double> t171 = t1 * t21 * bb[0];
	complex<double> t172 = t171 * t30;
	complex<double> t174 = t171 * t56;
	complex<double> t176 = t165 * t69;
	complex<double> t178 = t46 * cc[2];
	complex<double> t181 = t100 * t178 * bb[0] * dd[1];
	complex<double> t183 = t39 * cc[1];
	complex<double> t186 = t100 * t183 * bb[0] * dd[2];
	complex<double> t188 = t3 * cc[0];
	complex<double> t191 = t100 * t188 * bb[1] * dd[2];
	complex<double> t195 = t100 * t183 * bb[2] * dd[0];
	complex<double> t199 = t100 * t188 * bb[2] * dd[1];
	complex<double> t203 = t100 * t178 * bb[1] * dd[0];
	complex<double> t205 = -t154 / 0.18e2 + t158 / 0.18e2 + t162 / 0.72e2 - t166 / 0.72e2 - t168 / 0.72e2 + t172 / 0.72e2 - t174 / 0.72e2 + t176 / 0.72e2 + t181 / 0.18e2 - t186 / 0.18e2 + t191 / 0.18e2 + t195 / 0.18e2 - t199 / 0.18e2 - t203 / 0.18e2;
	complex<double> t220 = t72 / 0.24e2;
	complex<double> t221 = t76 / 0.24e2;
	complex<double> t222 = t78 / 0.24e2;
	complex<double> t223 = t81 / 0.24e2;
	complex<double> t224 = t83 / 0.24e2;
	complex<double> t225 = t85 / 0.24e2;
	complex<double> t226 = t8 / 0.12e2 - t12 / 0.12e2 - t19 / 0.12e2 - t26 / 0.12e2 + t31 / 0.12e2 + t36 / 0.12e2 - t43 / 0.12e2 + t50 / 0.12e2 - t57 / 0.12e2 + t61 / 0.12e2 - t65 / 0.12e2 + t70 / 0.12e2 + t220 - t221 - t222 - t223 - t224 + t225;
	complex<double> t227 = t88 / 0.24e2;
	complex<double> t228 = t90 / 0.24e2;
	complex<double> t229 = t92 / 0.24e2;
	complex<double> t230 = t94 / 0.24e2;
	complex<double> t231 = t96 / 0.24e2;
	complex<double> t232 = t98 / 0.24e2;
	complex<double> t239 = t162 / 0.24e2;
	complex<double> t240 = t166 / 0.24e2;
	complex<double> t241 = t168 / 0.24e2;
	complex<double> t242 = t172 / 0.24e2;
	complex<double> t243 = t174 / 0.24e2;
	complex<double> t244 = t176 / 0.24e2;
	complex<double> t245 = t227 - t228 + t229 + t230 - t231 + t232 - t138 / 0.12e2 + t141 / 0.12e2 - t145 / 0.12e2 + t148 / 0.12e2 - t154 / 0.12e2 + t158 / 0.12e2 + t239 - t240 - t241 + t242 - t243 + t244;
	complex<double> t246 = t226 + t245;
	complex<double> t247 = t94 - t96 - t90 + t92 - t83 - t168 + t98 + t72 + t85 + t162 - t78 - t166 - t76 - t81 + t176 + t88 - t174 + t172;
	complex<double> t248 = t22 * t118;
	complex<double> t249 = t15 * t107;
	complex<double> t250 = t15 * t127;
	complex<double> t251 = t53 * t127;
	complex<double> t252 = t165 * t68;
	complex<double> t253 = t15 * t112;
	complex<double> t254 = t161 * t34;
	complex<double> t255 = t53 * t118;
	complex<double> t256 = t171 * t55;
	complex<double> t257 = t171 * t29;
	complex<double> t258 = t22 * t122;
	complex<double> t259 = t53 * t103;
	complex<double> t260 = t22 * t107;
	complex<double> t261 = t22 * t103;
	complex<double> t262 = t165 * t24;
	complex<double> t263 = t53 * t112;
	complex<double> t264 = t161 * t17;
	complex<double> t265 = t15 * t122;
	complex<double> t266 = t248 - t249 + t250 + t251 - t252 - t253 - t254 + t255 + t256 - t257 + t258 - t259 - t260 - t261 + t262 - t263 + t264 + t265;
	complex<double> t268 = t151 * t146 * dd[2];
	complex<double> t271 = t134 * t146 * dd[1];
	complex<double> t274 = t134 * t135 * dd[0];
	complex<double> t277 = t156 * t135 * dd[2];
	complex<double> t279 = t142 * t112;
	complex<double> t281 = t139 * t122;
	complex<double> t295 = t268 / 0.4e1 - t271 / 0.4e1 + t274 / 0.4e1 - t277 / 0.4e1 + t279 / 0.4e1 - t281 / 0.4e1 - t252 / 0.8e1 - t254 / 0.8e1 + t256 / 0.8e1 - t257 / 0.8e1 + t262 / 0.8e1 + t264 / 0.8e1 + t248 / 0.8e1 - t249 / 0.8e1 + t250 / 0.8e1 + t251 / 0.8e1 - t253 / 0.8e1 + t255 / 0.8e1;
	complex<double> t303 = cc[0] * cc[2] * dd[1];
	complex<double> t304 = t15 * t303;
	complex<double> t307 = cc[0] * cc[1] * dd[2];
	complex<double> t308 = t53 * t307;
	complex<double> t312 = t38 * t39 * bb[0] * dd[2];
	complex<double> t316 = t2 * t3 * bb[1] * dd[2];
	complex<double> t318 = t22 * t307;
	complex<double> t321 = cc[1] * cc[2] * dd[0];
	complex<double> t322 = t15 * t321;
	complex<double> t326 = t45 * t46 * bb[1] * dd[0];
	complex<double> t330 = t38 * t39 * bb[2] * dd[0];
	complex<double> t332 = t53 * t321;
	complex<double> t334 = t22 * t303;
	complex<double> t338 = t2 * t3 * bb[2] * dd[1];
	complex<double> t342 = t45 * t46 * bb[0] * dd[1];
	complex<double> t344 = t258 / 0.8e1 - t259 / 0.8e1 - t260 / 0.8e1 - t261 / 0.8e1 - t263 / 0.8e1 + t265 / 0.8e1 - t304 / 0.4e1 + t308 / 0.4e1 - t312 / 0.4e1 + t316 / 0.4e1 - t318 / 0.4e1 + t322 / 0.4e1 - t326 / 0.4e1 + t330 / 0.4e1 - t332 / 0.4e1 + t334 / 0.4e1 - t338 / 0.4e1 + t342 / 0.4e1;
	complex<double> t358 = t8 / 0.6e1 - t12 / 0.6e1 - t19 / 0.6e1 - t26 / 0.6e1 + t31 / 0.6e1 + t36 / 0.6e1 - t43 / 0.6e1 + t50 / 0.6e1 - t57 / 0.6e1 + t61 / 0.6e1 - t65 / 0.6e1 + t70 / 0.6e1 + t220;
	complex<double> t362 = -t221 - t222 - t223 - t224 + t225 + t227 - t228 + t229 + t230 - t231 + t232 + t104 / 0.6e1 + t108 / 0.6e1 + t113 / 0.6e1;
	complex<double> t372 = -t154 / 0.6e1 + t158 / 0.6e1 + t239 - t240 - t241 + t242 - t243 + t244 + t181 / 0.6e1 - t186 / 0.6e1 + t191 / 0.6e1 + t195 / 0.6e1 - t199 / 0.6e1 - t203 / 0.6e1;
	complex<double> t375 = t268 / 0.6e1;
	complex<double> t376 = t271 / 0.6e1;
	complex<double> t377 = t274 / 0.6e1;
	complex<double> t378 = t277 / 0.6e1;
	complex<double> t379 = t279 / 0.6e1;
	complex<double> t380 = t281 / 0.6e1;
	complex<double> t387 = t1 * t178;
	complex<double> t390 = t375 - t376 + t377 - t378 + t379 - t380 - t252 / 0.24e2 - t254 / 0.24e2 + t256 / 0.24e2 - t257 / 0.24e2 + t262 / 0.24e2 + t264 / 0.24e2 - t387 * t63 / 0.6e1;
	complex<double> t391 = t1 * t183;
	complex<double> t394 = t1 * t188;
	complex<double> t412 = t391 * t59 / 0.6e1 - t394 * t10 / 0.6e1 + t387 * t48 / 0.6e1 - t391 * t41 / 0.6e1 + t394 * t5 / 0.6e1 + t248 / 0.24e2 - t249 / 0.24e2 + t250 / 0.24e2 + t251 / 0.24e2 - t253 / 0.24e2 + t255 / 0.24e2 + t258 / 0.24e2 - t259 / 0.24e2 - t260 / 0.24e2;
	complex<double> t417 = t304 / 0.6e1;
	complex<double> t418 = t308 / 0.6e1;
	complex<double> t419 = t312 / 0.6e1;
	complex<double> t420 = t316 / 0.6e1;
	complex<double> t421 = t318 / 0.6e1;
	complex<double> t422 = t322 / 0.6e1;
	complex<double> t423 = t326 / 0.6e1;
	complex<double> t424 = t330 / 0.6e1;
	complex<double> t425 = t332 / 0.6e1;
	complex<double> t426 = t334 / 0.6e1;
	complex<double> t427 = -t261 / 0.24e2 - t263 / 0.24e2 + t265 / 0.24e2 - t417 + t418 - t419 + t420 - t421 + t422 - t423 + t424 - t425 + t426;
	complex<double> t428 = t1 * t46;
	complex<double> t431 = t1 * t39;
	complex<double> t434 = t1 * t3;
	complex<double> t443 = -t338 + t342 + t428 * t127 + t428 * t118 - t431 * t107 - t431 * t112 + t434 * t122 + t434 * t127 - t428 * t103 - t428 * t107 + t431 * t118 + t431 * t122 - t434 * t112 - t434 * t103;
	complex<double> t445 = t390 + t412 + t427 + t443 / 0.6e1;
	complex<double> t458 = t375 - t376 + t377 - t378 + t379 - t380 - t252 / 0.12e2 - t254 / 0.12e2 + t256 / 0.12e2 - t257 / 0.12e2 + t262 / 0.12e2 + t264 / 0.12e2 + t248 / 0.12e2 - t249 / 0.12e2 + t250 / 0.12e2 + t251 / 0.12e2 - t253 / 0.12e2 + t255 / 0.12e2;
	complex<double> t467 = t258 / 0.12e2 - t259 / 0.12e2 - t260 / 0.12e2 - t261 / 0.12e2 - t263 / 0.12e2 + t265 / 0.12e2 - t417 + t418 - t419 + t420 - t421 + t422 - t423 + t424 - t425 + t426 - t338 / 0.6e1 + t342 / 0.6e1;
	//
	c[0] = t74 + t115 + t149 / 0.18e2 + t205;
	c[1] = t246;
	c[2] = t247 / 0.24e2;
	c[3] = t266 / 0.8e1;
	c[4] = t246;
	c[5] = t246;
	c[6] = t295 + t344;
	c[7] = t358 + t362 + t149 / 0.6e1 + t372;
	c[8] = t445;
	c[9] = t458 + t467;
	c[10] = t445;
	//
	t1 =  ko *  ko;
	t4 =  j *  j;
	complex<double> t9 = 8.0 / t1 / ko / t4 / j * c[0];
	t12 = 1 / ko / j;
	t16 = 3.0 * t12 * c[10];
	t18 = 2.0 * t12 * c[7];
	t21 = 1.0 / t1 / t4;
	t23 = 3.0 * t21 * c[4];
	t25 = 2.0 * t12 * c[1];
	complex<double> t27 = 2.0 * t12 * c[5];
	t29 = 2.0 * t12 * c[3];
	t31 = 2.0 * t12 * c[2];
	t33 = 3.0 * t21 * c[8];
	t35 = 2.0 * t12 * c[6];
	complex<double> t37 = 3.0 * t21 * c[10];
	t39 = 3.0 * t21 * c[9];
	t45 = 8.0 * t21 * c[0];
	t47 = 4.0 * t12 * c[0];
	//
	coef[0] = -t9;
	coef[1] = 3.0 * t12 * c[9];
	coef[2] = c[9];
	coef[3] = t9;
	coef[4] = t16;
	coef[5] = t18;
	coef[6] = t9;
	coef[7] = -t9;
	coef[8] = -c[10];
	coef[9] = -t23;
	coef[10] = -t25;
	coef[11] = -t27;
	coef[12] = -t29;
	coef[13] = -t31;
	coef[14] = -t16;
	coef[15] = -t33;
	coef[16] = t35;
	coef[17] = c[10];
	coef[18] = -t37;
	coef[19] = -t37;
	coef[20] = t37;
	coef[21] = -t35;
	coef[22] = t37;
	coef[23] = c[4];
	coef[24] = c[8];
	coef[25] = c[1];
	coef[26] = c[3];
	coef[27] = t39;
	coef[28] = 3.0 * t12 * c[4];
	coef[29] = 3.0 * t12 * c[8];
	coef[30] = c[2];
	coef[31] = c[5];
	coef[32] = -t18;
	coef[33] = c[6];
	coef[34] = c[7];
	coef[35] = -c[0];
	coef[36] = t23;
	coef[37] = t33;
	coef[38] = t25;
	coef[39] = -t39;
	coef[40] = t29;
	coef[41] = -t45;
	coef[42] = t45;
	coef[43] = t47;
	coef[44] = -t47;
	coef[45] = t31;
	coef[46] = c[0];
	coef[47] = t27;
	//
	t1 = 0.1e1 / AreaT;
	t2 = pow(bb[2], 0.2e1);
	t4 = t1 * t2 * bb[2];
	t5 = cc[0] * dd[1];
	t6 = t4 * t5;
	t8 = pow(bb[1], 0.2e1);
	t10 = t1 * t8 * bb[1];
	t11 = cc[0] * dd[2];
	t12 = t10 * t11;
	t14 = pow(bb[0], 0.2e1);
	t16 = t1 * t14 * bb[0];
	t17 = cc[1] * dd[2];
	t18 = t16 * t17;
	complex<double> t20 = cc[1] * dd[0];
	t21 = t4 * t20;
	t23 = cc[2] * dd[0];
	t24 = t10 * t23;
	t26 = cc[2] * dd[1];
	t27 = t16 * t26;
	t29 = t1 * cc[1];
	t30 = t29 * cc[2];
	t31 = bb[1] * bb[0];
	t33 = t30 * t31 * dd[1];
	t35 = t1 * cc[0];
	t36 = t35 * cc[2];
	t38 = t36 * t31 * dd[0];
	t40 = bb[2] * bb[1];
	t42 = t36 * t40 * dd[2];
	complex<double> t44 = t29 * cc[0];
	t46 = t44 * t40 * dd[1];
	t48 = bb[2] * bb[0];
	t50 = t30 * t48 * dd[2];
	t53 = t44 * t48 * dd[0];
	t55 = t1 * bb[1];
	t56 = pow(cc[1], 0.2e1);
	t59 = t55 * t56 * bb[2] * dd[0];
	t61 = t1 * t8;
	t64 = t61 * cc[2] * cc[1] * dd[0];
	complex<double> t66 = t1 * t14;
	t69 = t66 * cc[2] * cc[0] * dd[1];
	complex<double> t71 = t1 * bb[0];
	t72 = pow(cc[0], 0.2e1);
	t75 = t71 * t72 * bb[2] * dd[1];
	complex<double> t77 = t1 * t72;
	t78 = t31 * dd[2];
	complex<double> t79 = t77 * t78;
	complex<double> t82 = cc[1] * bb[2] * dd[0];
	t83 = t61 * t82;
	t85 = -t6 / 0.8e1 + t12 / 0.8e1 - t18 / 0.8e1 + t21 / 0.8e1 - t24 / 0.8e1 + t27 / 0.8e1 - t33 / 0.4e1 + t38 / 0.4e1 - t42 / 0.4e1 + t46 / 0.4e1 + t50 / 0.4e1 - t53 / 0.4e1 - t59 / 0.4e1 + t64 / 0.4e1 - t69 / 0.4e1 + t75 / 0.4e1 - t79 / 0.4e1 + t83 / 0.8e1;
	t87 = cc[2] * bb[1] * dd[0];
	t88 = t66 * t87;
	complex<double> t91 = cc[2] * bb[0] * dd[1];
	t92 = t61 * t91;
	t94 = pow(cc[2], 0.2e1);
	complex<double> t95 = t1 * t94;
	complex<double> t97 = t95 * t40 * dd[0];
	complex<double> t99 = t1 * cc[2];
	complex<double> t102 = t99 * cc[1] * t2 * dd[0];
	complex<double> t105 = cc[0] * bb[2] * dd[1];
	complex<double> t106 = t66 * t105;
	t110 = t35 * cc[1] * t14 * dd[2];
	t112 = t1 * t2;
	t113 = t112 * t91;
	complex<double> t116 = cc[0] * bb[1] * dd[2];
	complex<double> t117 = t112 * t116;
	complex<double> t119 = t61 * t105;
	t122 = cc[1] * bb[0] * dd[2];
	complex<double> t123 = t61 * t122;
	complex<double> t126 = t95 * t48 * dd[1];
	complex<double> t130 = t99 * cc[0] * t2 * dd[1];
	t134 = t29 * cc[0] * t8 * dd[2];
	t136 = t1 * t56;
	complex<double> t137 = t136 * t78;
	t139 = t66 * t116;
	t141 = t66 * t82;
	t143 = t112 * t122;
	t145 = t112 * t87;
	complex<double> t147 = -t88 / 0.8e1 + t92 / 0.8e1 + t97 / 0.4e1 - t102 / 0.4e1 - t106 / 0.8e1 + t110 / 0.4e1 + t113 / 0.8e1 + t117 / 0.8e1 - t119 / 0.8e1 - t123 / 0.8e1 - t126 / 0.4e1 + t130 / 0.4e1 - t134 / 0.4e1 + t137 / 0.4e1 + t139 / 0.8e1 + t141 / 0.8e1 - t143 / 0.8e1 - t145 / 0.8e1;
	t149 = t119 - t83 - t141 + t123 - t113 - t117 + t88 - t12 + t18 - t21 + t106 + t143 + t24 - t139 - t92 + t145 + t6 - t27;
	t150 = sqrt(0.3e1);
	t151 = dd[2] * t150;
	complex<double> t153 = t36 * t40 * t151;
	t154 = dd[0] * t150;
	t156 = t44 * t48 * t154;
	complex<double> t157 = dd[1] * t150;
	complex<double> t159 = t44 * t40 * t157;
	t161 = t30 * t48 * t151;
	complex<double> t163 = t36 * t31 * t154;
	t165 = t30 * t31 * t157;
	complex<double> t167 = bb[1] * dd[2];
	t168 = t167 * t150;
	complex<double> t169 = t136 * bb[0] * t168;
	t172 = t44 * t8 * dd[2] * t150;
	complex<double> t175 = t44 * t14 * dd[2] * t150;
	complex<double> t177 = bb[2] * dd[0];
	t178 = t177 * t150;
	complex<double> t179 = t95 * bb[1] * t178;
	complex<double> t180 = t61 * cc[2];
	t181 = t20 * t150;
	complex<double> t182 = t180 * t181;
	t183 = t66 * cc[2];
	complex<double> t184 = t5 * t150;
	complex<double> t185 = t183 * t184;
	complex<double> t187 = t55 * t56 * t178;
	complex<double> t189 = bb[0] * dd[2];
	complex<double> t190 = t189 * t150;
	t191 = t77 * bb[1] * t190;
	complex<double> t193 = bb[2] * dd[1];
	complex<double> t194 = t193 * t150;
	t195 = t71 * t72 * t194;
	complex<double> t198 = t36 * t2 * dd[1] * t150;
	complex<double> t201 = t30 * t2 * dd[0] * t150;
	t203 = t95 * bb[0] * t194;
	complex<double> t204 = t153 + t156 - t159 - t161 - t163 + t165 - t169 + t172 - t175 - t179 - t182 + t185 + t187 + t191 - t195 - t198 + t201 + t203;
	t205 = bb[0] * dd[1];
	complex<double> t206 = t205 * t150;
	complex<double> t207 = t180 * t206;
	complex<double> t208 = t66 * cc[0];
	complex<double> t209 = t208 * t194;
	complex<double> t211 = t66 * cc[1] * t178;
	complex<double> t212 = t112 * cc[2];
	complex<double> t213 = t212 * t206;
	complex<double> t214 = bb[1] * dd[0];
	complex<double> t215 = t214 * t150;
	complex<double> t216 = t183 * t215;
	complex<double> t217 = t61 * cc[1];
	complex<double> t218 = t217 * t178;
	complex<double> t219 = t217 * t190;
	t221 = t61 * cc[0] * t194;
	t223 = t112 * cc[1] * t190;
	t224 = t212 * t215;
	t226 = t112 * cc[0] * t168;
	t227 = t208 * t168;
	t229 = t16 * t26 * t150;
	t231 = t10 * t23 * t150;
	complex<double> t233 = t10 * t11 * t150;
	complex<double> t235 = t16 * t17 * t150;
	complex<double> t236 = t4 * t181;
	complex<double> t237 = t4 * t184;
	complex<double> t238 = -t207 + t209 - t211 - t213 + t216 - t218 + t219 + t221 + t223 + t224 - t226 - t227 - t229 + t231 - t233 + t235 - t236 + t237;
	t239 = t204 / 0.12e2 + t238 / 0.24e2;
	t240 = -t153 - t156 + t159 + t161 + t163 - t165 + t169 - t172 + t175 + t179 + t182 - t185 - t187;
	t255 = -t191 / 0.18e2 + t195 / 0.18e2 + t198 / 0.18e2 - t201 / 0.18e2 - t203 / 0.18e2 + t207 / 0.72e2 - t209 / 0.72e2 + t211 / 0.72e2 + t213 / 0.72e2 - t216 / 0.72e2 + t218 / 0.72e2 - t219 / 0.72e2 - t221 / 0.72e2 - t223 / 0.72e2;
	t260 = t150 * t1;
	t261 = t260 * t72;
	t262 = t261 * t122;
	t264 = t260 * t56;
	t265 = t264 * t105;
	complex<double> t267 = t264 * t116;
	complex<double> t269 = t264 * t87;
	t271 = t260 * t94;
	complex<double> t272 = t271 * t122;
	t274 = t261 * t91;
	complex<double> t276 = t261 * t87;
	complex<double> t278 = t271 * t82;
	complex<double> t280 = t271 * t105;
	complex<double> t282 = t261 * t82;
	complex<double> t284 = -t224 / 0.72e2 + t226 / 0.72e2 + t227 / 0.72e2 - t262 / 0.18e2 - t265 / 0.18e2 + t267 / 0.18e2 - t269 / 0.18e2 - t272 / 0.18e2 + t274 / 0.18e2 - t276 / 0.18e2 + t278 / 0.18e2 - t280 / 0.18e2 + t282 / 0.18e2;
	complex<double> t285 = t264 * t91;
	complex<double> t287 = t271 * t116;
	t295 = t72 * cc[0];
	complex<double> t298 = t260 * t295 * bb[2] * dd[1];
	complex<double> t300 = t94 * cc[2];
	t303 = t260 * t300 * bb[1] * dd[0];
	complex<double> t305 = t56 * cc[1];
	t308 = t260 * t305 * bb[2] * dd[0];
	t312 = t260 * t295 * bb[1] * dd[2];
	t316 = t260 * t305 * bb[0] * dd[2];
	complex<double> t320 = t260 * t300 * bb[0] * dd[1];
	t322 = t285 / 0.18e2 + t287 / 0.18e2 + t229 / 0.72e2 - t231 / 0.72e2 + t233 / 0.72e2 - t235 / 0.72e2 + t236 / 0.72e2 - t237 / 0.72e2 - t298 / 0.18e2 - t303 / 0.18e2 + t308 / 0.18e2 + t312 / 0.18e2 - t316 / 0.18e2 + t320 / 0.18e2;
	complex<double> t331 = t1 * t300;
	t334 = t1 * t305;
	complex<double> t337 = t1 * t295;
	complex<double> t346 = t33 / 0.6e1;
	complex<double> t347 = t6 / 0.24e2 - t12 / 0.24e2 + t18 / 0.24e2 - t21 / 0.24e2 + t24 / 0.24e2 - t27 / 0.24e2 - t331 * t205 / 0.6e1 + t334 * t189 / 0.6e1 - t337 * t167 / 0.6e1 + t331 * t214 / 0.6e1 - t334 * t177 / 0.6e1 + t337 * t193 / 0.6e1 + t346;
	complex<double> t348 = t38 / 0.6e1;
	complex<double> t349 = t42 / 0.6e1;
	complex<double> t350 = t46 / 0.6e1;
	complex<double> t351 = t50 / 0.6e1;
	complex<double> t352 = t53 / 0.6e1;
	complex<double> t353 = t59 / 0.6e1;
	complex<double> t354 = t64 / 0.6e1;
	complex<double> t355 = t69 / 0.6e1;
	complex<double> t356 = t75 / 0.6e1;
	complex<double> t357 = t79 / 0.6e1;
	complex<double> t361 = t97 / 0.6e1;
	t362 = -t348 + t349 - t350 - t351 + t352 + t353 - t354 + t355 - t356 + t357 - t83 / 0.24e2 + t88 / 0.24e2 - t92 / 0.24e2 - t361;
	complex<double> t364 = t102 / 0.6e1;
	complex<double> t366 = t110 / 0.6e1;
	complex<double> t371 = t126 / 0.6e1;
	t372 = t130 / 0.6e1;
	complex<double> t373 = t134 / 0.6e1;
	complex<double> t374 = t137 / 0.6e1;
	t377 = t364 + t106 / 0.24e2 - t366 - t113 / 0.24e2 - t117 / 0.24e2 + t119 / 0.24e2 + t123 / 0.24e2 + t371 - t372 + t373 - t374 - t139 / 0.24e2 - t141 / 0.24e2;
	complex<double> t404 = t143 / 0.24e2 + t145 / 0.24e2 - t77 * t82 / 0.6e1 - t77 * t91 / 0.6e1 + t136 * t87 / 0.6e1 + t77 * t87 / 0.6e1 - t95 * t116 / 0.6e1 - t95 * t82 / 0.6e1 + t136 * t105 / 0.6e1 + t95 * t122 / 0.6e1 + t95 * t105 / 0.6e1 - t136 * t116 / 0.6e1 - t136 * t91 / 0.6e1 + t77 * t122 / 0.6e1;
	complex<double> t406 = t347 + t362 + t377 + t404;
	complex<double> t414 = -t6 / 0.12e2 + t12 / 0.12e2 - t18 / 0.12e2 + t21 / 0.12e2 - t24 / 0.12e2 + t27 / 0.12e2 - t346 + t348 - t349 + t350 + t351 - t352 - t353 + t354 - t355 + t356 - t357 + t83 / 0.12e2;
	t426 = -t88 / 0.12e2 + t92 / 0.12e2 + t361 - t364 - t106 / 0.12e2 + t366 + t113 / 0.12e2 + t117 / 0.12e2 - t119 / 0.12e2 - t123 / 0.12e2 - t371 + t372 - t373 + t374 + t139 / 0.12e2 + t141 / 0.12e2 - t143 / 0.12e2 - t145 / 0.12e2;
	complex<double> t442 = -t191 / 0.6e1 + t195 / 0.6e1 + t198 / 0.6e1 - t201 / 0.6e1 - t203 / 0.6e1 + t207 / 0.24e2 - t209 / 0.24e2 + t211 / 0.24e2 + t213 / 0.24e2 - t216 / 0.24e2 + t218 / 0.24e2 - t219 / 0.24e2 - t221 / 0.24e2 - t223 / 0.24e2;
	complex<double> t457 = -t224 / 0.24e2 + t226 / 0.24e2 + t227 / 0.24e2 - t262 / 0.6e1 - t265 / 0.6e1 + t267 / 0.6e1 - t269 / 0.6e1 - t272 / 0.6e1 + t274 / 0.6e1 - t276 / 0.6e1 + t278 / 0.6e1 - t280 / 0.6e1 + t282 / 0.6e1;
	complex<double> t472 = t285 / 0.6e1 + t287 / 0.6e1 + t229 / 0.24e2 - t231 / 0.24e2 + t233 / 0.24e2 - t235 / 0.24e2 + t236 / 0.24e2 - t237 / 0.24e2 - t298 / 0.6e1 - t303 / 0.6e1 + t308 / 0.6e1 + t312 / 0.6e1 - t316 / 0.6e1 + t320 / 0.6e1;
	//
	cm[0] = t85 + t147;
	cm[1] = t149 / 0.8e1;
	cm[2] = t239;
	cm[3] = t239;
	cm[4] = t240 / 0.18e2 + t255 + t284 + t322;
	cm[5] = -t238 / 0.24e2;
	cm[6] = t406;
	cm[7] = t406;
	cm[8] = t414 + t426;
	cm[9] = t240 / 0.6e1 + t442 + t457 + t472;
	cm[10] = t239;
	//
	t3 = 0.1e1 /  ko /  j;
	t7 = 2.0 * t3 * cm[9];
	t8 =  ko *  ko;
	t10 =  j *  j;
	t12 = 0.1e1 /  t8 / t10;
	t14 = 3.0 * t12 * cm[7];
	t16 = 2.0 * t3 * cm[0];
	t18 = 3.0 * t3 * cm[6];
	t20 = 2.0 * t3 * cm[5];
	t27 = 8.0 / t8 / ko / t10 / j * cm[4];
	t29 = 2.0 * t3 * cm[10];
	t31 = 2.0 * t3 * cm[3];
	t33 = 2.0 * t3 * cm[1];
	t35 = 3.0 * t12 * cm[2];
	t37 = 3.0 * t12 * cm[8];
	t43 = 3.0 * t12 * cm[6];
	t45 = 4.0 * t3 * cm[4];
	t47 = 8.0 * t12 * cm[4];
	//
	coefm[0] = cm[8];
	coefm[1] = 3.0 * t3 * cm[8];
	coefm[2] = t7;
	coefm[3] = -t14;
	coefm[4] = t16;
	coefm[5] = t18;
	coefm[6] = -t20;
	coefm[7] = -cm[6];
	coefm[8] = cm[6];
	coefm[9] = -t27;
	coefm[10] = -t29;
	coefm[11] = -t18;
	coefm[12] = t27;
	coefm[13] = -t31;
	coefm[14] = -t33;
	coefm[15] = -t35;
	coefm[16] = -t27;
	coefm[17] = -t16;
	coefm[18] = cm[7];
	coefm[19] = cm[10];
	coefm[20] = cm[1];
	coefm[21] = cm[2];
	coefm[22] = cm[5];
	coefm[23] = t37;
	coefm[24] = 3.0 * t3 * cm[2];
	coefm[25] = 3.0 * t3 * cm[7];
	coefm[26] = cm[3];
	coefm[27] = t43;
	coefm[28] = -t7;
	coefm[29] = -t43;
	coefm[30] = t33;
	coefm[31] = t45;
	coefm[32] = t35;
	coefm[33] = t43;
	coefm[34] = -t43;
	coefm[35] = t27;
	coefm[36] = cm[0];
	coefm[37] = cm[4];
	coefm[38] = cm[9];
	coefm[39] = -cm[4];
	coefm[40] = -t45;
	coefm[41] = -t37;
	coefm[42] = t20;
	coefm[43] = t47;
	coefm[44] = t14;
	coefm[45] = t31;
	coefm[46] = t29;
	coefm[47] = -t47;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g1_f1
// ***********************************************************************

complex<double> X_function_g1_f1 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g1_f2
// ***********************************************************************

complex<double> X_function_g1_f2 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g1_f3
// ***********************************************************************

complex<double> X_function_g1_f3 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g2_f1
// ***********************************************************************

complex<double> X_function_g2_f1 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g2_f2
// ***********************************************************************

complex<double> X_function_g2_f2 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g2_f3
// ***********************************************************************

complex<double> X_function_g2_f3 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g3_f1
// ***********************************************************************

complex<double> X_function_g3_f1 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g3_f2
// ***********************************************************************

complex<double> X_function_g3_f2 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_g3_f3
// ***********************************************************************

complex<double> X_function_g3_f3 (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
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
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF void X_function_pre
// ***********************************************************************

void X_function_pre (double theta, double Psi, double tPsiA, double tPsiB, double PsiA, double PsiB, double B, double Bm, complex<double> N[], complex<double> Nm[], const double ko)
{
	double  D, D1, D2;
	//
	if (Psi >= PsiB && Psi >= PsiA)
	{
		D  = sqrt(3.0) / sin(Psi);
		X2(D, ko, B, N);
		X2 (D, ko, Bm, Nm);
	}
	else if (theta >= M_PI_2)
	{
		D  = sqrt(3.0) / (cos(Psi) * tPsiA);
		X2 (D, ko, B, N);
		X2 (D, ko, Bm, Nm);
	}
	else if (Psi >= PsiA)
	{
		D1 = 2.0 * sqrt(3.0) / (cos(Psi) * ( tPsiB + tan(Psi) )  );
		D2 = sin(Psi) / sqrt(3.0);
		X1( Psi, D1, D2, tPsiB, ko, B,  N);
		X1( Psi, D1, D2, tPsiB, ko, Bm, Nm);
	}
	else
	{
		D1 = sqrt(3.0) / (cos(Psi) * sin(theta) );
		D2 = ( cos(Psi) * tPsiA ) / sqrt(3.0);
		X1( Psi, D1, D2, tPsiB, ko, B, N);
		X1( Psi, D1, D2, tPsiB, ko, Bm, Nm);
	}

}

// ***********************************************************************
//			IMPLEMENTATION OF void void X1
// ***********************************************************************

void X1 ( double psi, complex<double> D1, complex<double> D2, double tpsiB, const double ko, complex<double> B, complex<double> N[] )
{
	complex<double> D3, aD3, expaD3, aD1, expaD1, D2a, H, T11, T21, T31, T41, T12, T22, T32, T42;
	complex<double>	N21, N31, N81, N91, N101, N111, N121, N41, N51, N22, N32, N82, N92, N102, N112, N122, N42, N52;
	//
	complex<double> a  = Iunit * ko * B;
	//
	D3     = sqrt(3.0) / (cos(psi) * tpsiB);
	aD3    = a * D3;
	expaD3 = exp(-aD3);
	aD1    = a * D1;
	expaD1 = exp(-aD1);
	D2a    = D2 / a;
	//
	H      = 1.0 - D1 * D2;
	//
	T11    = (expaD3 - expaD1) / aD3;
	T21    = (T11 - H * expaD1) / aD3;
	T31    = (2.0 * T21 - pow(H,2) * expaD1) / aD3;
	T41    = (3.0 * T31 - pow(H,3) * expaD1) / aD3;
	//
	T12    = D2a * (1.0 - expaD1);
	T22    = D2a * (1.0 - H * expaD1 - T12);
	T32    = D2a * (1.0 - pow(H,2) * expaD1 - 2.0 * T22);
	T42    = D2a * (1.0 - pow(H,3) * expaD1 - 3.0 * T32);
	//
	N21    = T11;
	N31    = D3 * (T11 + T21);
	N81    = T21;
	N91    = T31;
	N101   = D3 * (T21 + T31);
	N111   = D3 * (T31 + T41);
	N121   = D3 * (N101 + N111);
	N41    = D3 * (N31 + N101);
	N51    = D3 * (N41 + N121);
	//
	N22    = T12;
	N32    = (T12 - T22) / D2;
	N82    = T22;
	N92    = T32;
	N102   = (T22 - T32) / D2;
	N112   = (T32 - T42) / D2;
	N122   = (N102 - N112) / D2;
	N42    = (N32 - N102) / D2;
	N52    = (N42 - N122) / D2;
	// Final Output
	N[0]   = 1.0;
	N[1]   = N21+N22;
	N[2]   = N31+N32;
	N[3]   = N41+N42;
	N[4]   = N51+N52;
	N[5]   = 0.5;
	N[6]   = 1.0/3.0;
	N[7]   = N81+N82;
	N[8]   = N91+N92;
	N[9]   = N101+N102;
	N[10]  = N111+N112;
	N[11]  = N121+N122;
}

// ***********************************************************************
//			IMPLEMENTATION OF void void X2
// ***********************************************************************

void X2 ( complex<double> D, const double ko, complex<double> B, complex<double> N[] )
{
	complex<double> aD, expaD, T1, T2, T3, T4;
	//
	complex<double> a  = Iunit * ko * B;
	//
	aD    = a * D;
	expaD = exp(-aD);
	//
	T1    = (1.0 - expaD) / aD;
	T2    = (1.0 - T1) / aD;
	T3    = (1.0 - 2.0 * T2) / aD;
	T4    = (1.0 - 3.0 * T3) / aD;
	// Final Output
	N[0]   = 1.0;
	N[1]   = T1;
	N[2]   = D * (T1 - T2);
	N[3]   = 0.0;
	N[4]   = 0.0;
	N[5]   = 0.5;
	N[6]   = 1.0/3.0;
	N[7]   = T2;
	N[8]   = T3;
	N[9]   = D * (T2 - T3);
	N[10]  = D * (T3 - T4);
	N[11]  = 0.0;
	//
	N[11]  = D * (N[9] - N[10]);
	N[3]   = D * (N[2] - N[9]);
	N[4]   = D * (N[3] - N[11]);
}

// ***********************************************************************
//			IMPLEMENTATION OF void PSI_limits
// ***********************************************************************

void PSI_limits ( int argument, double PsiA, double PsiB, double *psi_A, double *psi_B )
{
	switch(argument)
	{
		case 1:
			 *psi_A = PsiB;
			 *psi_B = M_PI_2;
			 break;
		case 2:
			 *psi_A = PsiA;
			 *psi_B = M_PI_2;
			 break;
		case 3:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
		case 4:
			 *psi_A = 0;
			 *psi_B = PsiB;
			 break;
		case 5:
			 *psi_A = PsiA;
			 *psi_B = PsiB;
			 break;
		case 6:
			 *psi_A = 0;
			 *psi_B = PsiA;
			 break;
	}
}

// ***********************************************************************
//			IMPLEMENTATION OF void THETA_limits
// ***********************************************************************

void THETA_limits ( int argument, double *theta_A, double *theta_B )
{
	switch(argument)
	{
		case 1:
			 *theta_A = 0.0;
			 *theta_B = M_PI_2 ;
			 break;
		case 2:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 3:
			 *theta_A = M_PI_2;
			 *theta_B = M_PI;
			 break;
		case 4:
			 *theta_A = 0.0;
			 *theta_B = M_PI / 3;
			 break;
		case 5:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
		case 6:
			 *theta_A = M_PI / 3;
			 *theta_B = M_PI_2;
			 break;
	}
}

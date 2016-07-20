/**************************************************************************************************************************
           
		           DIRECT_SS_EA_RWG.cpp

Main body of the DIRECT EVALUATION method for the evaluation of the edge adjacent 4-D
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
  r1,r2,r3, r4 = point vectors of the triangular element's vertices
  Outer triangle P:(rp1,rp2,rp3)=(r1,r2,r3)
  Inner triangle Q:(rq1,rq2,rq3)=(r2,r1,r4)
  N_theta = order of the Gauss-Legendre cubature for the 1-D smooth integral over theta
  N_psi   = order of the Gauss-Legendre cubature for the 1-D smooth integral over Psi
  ko = wavenumber

  OUTPUT DATA
  I_DE[9]   = 4-D strongly singular integrals
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

#include "DIRECT_SS_EA_RWG.h"
#include <math.h>
using namespace std;

// ***********************************************************************
//			IMPLEMENTATION OF void DIRECT
// ***********************************************************************

void DIRECT_SS_EA_RWG (const double r1[],const double r2[],const double r3[],const double r4[], const double ko, const int N_theta, const int N_psi, const double z_theta[], const double w_theta[],  const double z_psi[], const double w_psi[], complex<double> I_DE[] )
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
	complex<double> Ising_f1_f1, Ising_f1_f3, Ising_f2_f2, Ising_f2_f3, Ising_f3_f1, Ising_f3_f2, Ising_f3_f3;
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
	complex<double> I_theta_f1_f1, I_theta_f1_f3, I_theta_f2_f2, I_theta_f2_f3, I_theta_f3_f1, I_theta_f3_f2, I_theta_f3_f3;
	complex<double> I_psi_f1_f1, I_psi_f1_f3, I_psi_f2_f2, I_psi_f2_f3, I_psi_f3_f1, I_psi_f3_f2, I_psi_f3_f3;
	//
	complex<double> X_f1_f1, X_f1_f3, X_f2_f2, X_f2_f3, X_f3_f1, X_f3_f2, X_f3_f3;
	//
	complex<double> I_f1_f1[6], I_f1_f3[6], I_f2_f2[6], I_f2_f3[6], I_f3_f1[6], I_f3_f2[6], I_f3_f3[6]; 
	//
	complex<double> If1_f1, If1_f3, If2_f2, If2_f3, If3_f1, If3_f2, If3_f3; 

	// 2. Coefficients' parameters

	complex<double> coef_f1_f1[4],  coefm_f1_f1[4];
	complex<double> coef_f1_f3[14], coefm_f1_f3[14];
	complex<double> coef_f2_f2[4],  coefm_f2_f2[4];
	complex<double> coef_f2_f3[14], coefm_f2_f3[14];
	complex<double> coef_f3_f1[10], coefm_f3_f1[10];
	complex<double> coef_f3_f2[10], coefm_f3_f2[10];
	complex<double> coef_f3_f3[13], coefm_f3_f3[13];

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
	coefficients_f1_f1 ( r1,  r2,  r3,  r4,  ko,  coef_f1_f1,  coefm_f1_f1 );
	coefficients_f1_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f1_f3,  coefm_f1_f3 );
	coefficients_f2_f2 ( r1,  r2,  r3,  r4,  ko,  coef_f2_f2,  coefm_f2_f2 );
	coefficients_f2_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f2_f3,  coefm_f2_f3 );
	coefficients_f3_f1 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f1,  coefm_f3_f1 );
	coefficients_f3_f2 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f2,  coefm_f3_f2 );
	coefficients_f3_f3 ( r1,  r2,  r3,  r4,  ko,  coef_f3_f3,  coefm_f3_f3 );
     // Initialization of I_
	 for ( int im = 0; im <  6; im++ )
	 {
		 I_f1_f1[im] = (0.0,0.0);
		 I_f1_f3[im] = (0.0,0.0);
		 I_f2_f2[im] = (0.0,0.0);
		 I_f2_f3[im] = (0.0,0.0);
		 I_f3_f1[im] = (0.0,0.0);
		 I_f3_f2[im] = (0.0,0.0);
		 I_f3_f3[im] = (0.0,0.0);
	 }
	 //
	 for ( int m = 1; m <  7; m++ )
	 {
		 I_theta_f1_f1 = (0.0,0.0);
		 I_theta_f1_f3 = (0.0,0.0);
		 I_theta_f2_f2 = (0.0,0.0);
		 I_theta_f2_f3 = (0.0,0.0);
		 I_theta_f3_f1 = (0.0,0.0);
		 I_theta_f3_f2 = (0.0,0.0);
		 I_theta_f3_f3 = (0.0,0.0);
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
			 I_psi_f1_f1 = (0.0,0.0);
			 I_psi_f1_f3 = (0.0,0.0);
			 I_psi_f2_f2 = (0.0,0.0);
			 I_psi_f2_f3 = (0.0,0.0);
			 I_psi_f3_f1 = (0.0,0.0);
			 I_psi_f3_f2 = (0.0,0.0);
			 I_psi_f3_f3 = (0.0,0.0);
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
				 X_f1_f1 = X_function_f1_f1 (THETA, PSI, B, Bm, coef_f1_f1, coefm_f1_f1, N, Nm);
				 X_f1_f3 = X_function_f1_f3 (THETA, PSI, B, Bm, coef_f1_f3, coefm_f1_f3, N, Nm);
				 X_f2_f2 = X_function_f2_f2 (THETA, PSI, B, Bm, coef_f2_f2, coefm_f2_f2, N, Nm);
				 X_f2_f3 = X_function_f2_f3 (THETA, PSI, B, Bm, coef_f2_f3, coefm_f2_f3, N, Nm);
				 X_f3_f1 = X_function_f3_f1 (THETA, PSI, B, Bm, coef_f3_f1, coefm_f3_f1, N, Nm);
				 X_f3_f2 = X_function_f3_f2 (THETA, PSI, B, Bm, coef_f3_f2, coefm_f3_f2, N, Nm);
				 X_f3_f3 = X_function_f3_f3 (THETA, PSI, B, Bm, coef_f3_f3, coefm_f3_f3, N, Nm);
				 // Integrate
				 WPSI = w_psi[n_psi];
				 //
				 I_psi_f1_f1 +=  WPSI * X_f1_f1;
				 I_psi_f1_f3 +=  WPSI * X_f1_f3;
				 I_psi_f2_f2 +=  WPSI * X_f2_f2;
				 I_psi_f2_f3 +=  WPSI * X_f2_f3;
				 I_psi_f3_f1 +=  WPSI * X_f3_f1;
				 I_psi_f3_f2 +=  WPSI * X_f3_f2;
				 I_psi_f3_f3 +=  WPSI * X_f3_f3;
			 }//end for ( int n_psi = 0 ; n_psi <  N_psi ; n_psi++ )
			 //
			 I_psi_f1_f1 *=  J_psi;
			 I_psi_f1_f3 *=  J_psi;
			 I_psi_f2_f2 *=  J_psi;
			 I_psi_f2_f3 *=  J_psi;
			 I_psi_f3_f1 *=  J_psi;
			 I_psi_f3_f2 *=  J_psi;
			 I_psi_f3_f3 *=  J_psi;
			 //
			 WTHETA = w_theta[n_theta];
			 //
			 I_theta_f1_f1 +=  WTHETA * I_psi_f1_f1;
			 I_theta_f1_f3 +=  WTHETA * I_psi_f1_f3;
			 I_theta_f2_f2 +=  WTHETA * I_psi_f2_f2;
			 I_theta_f2_f3 +=  WTHETA * I_psi_f2_f3;
			 I_theta_f3_f1 +=  WTHETA * I_psi_f3_f1;
			 I_theta_f3_f2 +=  WTHETA * I_psi_f3_f2;
			 I_theta_f3_f3 +=  WTHETA * I_psi_f3_f3;
		 } //end for ( int n_theta = 0 ; n_theta <  N_theta ; n_theta++ )
		 //
		 I_f1_f1[m-1] = J_theta * I_theta_f1_f1;
		 I_f1_f3[m-1] = J_theta * I_theta_f1_f3;
		 I_f2_f2[m-1] = J_theta * I_theta_f2_f2;
		 I_f2_f3[m-1] = J_theta * I_theta_f2_f3;
		 I_f3_f1[m-1] = J_theta * I_theta_f3_f1;
		 I_f3_f2[m-1] = J_theta * I_theta_f3_f2;
		 I_f3_f3[m-1] = J_theta * I_theta_f3_f3;
	 } //end for ( int m = 1; m <  7; m++ )
	 //
	 If1_f1 = I_f1_f1[0] + I_f1_f1[1] + I_f1_f1[2] + I_f1_f1[3] + I_f1_f1[4] + I_f1_f1[5];
	 If1_f3 = I_f1_f3[0] + I_f1_f3[1] + I_f1_f3[2] + I_f1_f3[3] + I_f1_f3[4] + I_f1_f3[5];
	 If2_f2 = I_f2_f2[0] + I_f2_f2[1] + I_f2_f2[2] + I_f2_f2[3] + I_f2_f2[4] + I_f2_f2[5];
	 If2_f3 = I_f2_f3[0] + I_f2_f3[1] + I_f2_f3[2] + I_f2_f3[3] + I_f2_f3[4] + I_f2_f3[5];
	 If3_f1 = I_f3_f1[0] + I_f3_f1[1] + I_f3_f1[2] + I_f3_f1[3] + I_f3_f1[4] + I_f3_f1[5];
	 If3_f2 = I_f3_f2[0] + I_f3_f2[1] + I_f3_f2[2] + I_f3_f2[3] + I_f3_f2[4] + I_f3_f2[5];
	 If3_f3 = I_f3_f3[0] + I_f3_f3[1] + I_f3_f3[2] + I_f3_f3[3] + I_f3_f3[4] + I_f3_f3[5];
	 //
	 double PRECOMP = Jp * Jq / (4 * Ap * Aq);
	 // FINAL OUTPUT
	 I_DE[0] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r41,r41)) * If1_f1;
	 I_DE[1] = 0.0;
	 I_DE[2] = PRECOMP * sqrt(vector_dot(r32,r32)) * sqrt(vector_dot(r21,r21)) * If1_f3;
	 I_DE[3] = 0.0;
	 I_DE[4] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r42,r42)) * If2_f2;
	 I_DE[5] = PRECOMP * sqrt(vector_dot(r31,r31)) * sqrt(vector_dot(r21,r21)) * If2_f3;
	 I_DE[6] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r41,r41)) * If3_f1;
	 I_DE[7] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r42,r42)) * If3_f2;
	 I_DE[8] = PRECOMP * sqrt(vector_dot(r21,r21)) * sqrt(vector_dot(r21,r21)) * If3_f3;
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f1
// ***********************************************************************

void coefficients_f1_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	complex<double> j   = Iunit;

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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f1_f3
// ***********************************************************************

void coefficients_f1_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
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
	t5  = cc[0] * bb[2];
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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f2
// ***********************************************************************

void coefficients_f2_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}

	complex<double> j = Iunit;
	//
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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f2_f3
// ***********************************************************************

void coefficients_f2_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f1
// ***********************************************************************

void coefficients_f3_f1 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	//
	for (int i = 0; i < 3; i++)
	{
			aa[i] = r1[i] - r1[i];
			bb[i] = r2[i] - r1[i];
			cc[i] = r3[i] - r1[i];
			dd[i] = r4[i] - r1[i];
	}
	//
	complex<double> t1  = cc[1] * dd[2];
	complex<double> t2  = sqrt(3.0);
	complex<double> t3  = bb[0] * t2;
	complex<double> t5  = cc[0] * dd[1];
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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f2
// ***********************************************************************

void coefficients_f3_f2 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
{
	double aa[3], bb[3], cc[3], dd[3];
	complex<double> c[3], cm[3];
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
}

// ***********************************************************************
//			IMPLEMENTATION OF void coefficients_f3_f3
// ***********************************************************************

void coefficients_f3_f3 (const double r1[],const double r2[],const double r3[],const double r4[],const double ko, complex<double> coef[], complex<double> coefm[] )
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
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f1_f1
// ***********************************************************************

complex<double> X_function_f1_f1 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2  = pow(B,2);
	double t3  = pow(t2,2);
	double t5  = 1.0 / t3 /B;
	double t9  = cos(Psi);
	double t10 = pow(t9,2);
	double t11 = sin(Psi);
	double t13 = sin(theta);
	double t14 = t10 * t11 * t13;
	double t20 = t11 * t13;
	double t39 = pow(Bm,2);
	double t40 = pow(t39,2);
	double t42 = 1.0 / t40 / Bm;
	//
	X = N[0] * t5 * coef[2] * t14 + N[2] / t3 * t10 * t20 * coef[0] + N[3] / t2 / B * t10 * t20 * coef[1]
	+ N[1] * t5 * t10 * t20 * coef[3] + Nm[0] * t42 * coefm[0] * t14 + Nm[2] / t40 * t10 * t20 * coefm[2]
	+ Nm[3] / t39 / Bm * t10 * t20 * coefm[3] + Nm[1] * t42 * t10 * t20 * coefm[1];
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f1_f3
// ***********************************************************************

complex<double> X_function_f1_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
    double t3 = t2 * t2;
    double t4 = 0.1e1 / t3;
    double t7 = cos(Psi);
    double t8 = sin(Psi);
    double t9 = t7 * t8;
    double t12 = 0.1e1 / t3 / B;
    double t15 = t7 * t7;
    double t17 = sin(theta);
    double t18 = t15 * t8 * t17;
    double t22 = cos(theta);
    double t24 = t15 * t22 * t8;
    double t40 = 0.1e1 / t2 / B;
    double t79 = Bm * Bm;
    double t80 = t79 * t79;
    double t81 = 0.1e1 / t80;
    double t86 = 0.1e1 / t80 / Bm;
    double t107 = 0.1e1 / t79 / Bm;
	//
	X = N[0] * (t4 * coef[0] * t9 + t12 * coef[9] * t18 + t12 * coef[8] * t24) + N[7] * t4 * t9 * coef[11] + N[5] * t4 * t9 * coef[6]
	+ N[3] * (t40 * coef[1] * t18 + t40 * coef[3] * t24) + N[9] * t40 * t9 * coef[2] + N[2] * (t4 * coef[4] * t24 + t4 * coef[5] * t18 + t40 * coef[10] * t9)
	+ N[1] * (t12 * coef[12] * t24 + t12 * coef[13] * t18 + t4 * coef[7] * t9) + Nm[0] * (t81 * coefm[0] * t9 + t86 * coefm[2] * t24 + t86 * coefm[3] * t18)
	+ Nm[5] * t81 * t9 * coefm[1] + Nm[7] * t81 * t9 * coefm[10] + Nm[3] * (t107 * coefm[8] * t24 + t107 * coefm[9] * t18)
	+ Nm[2] * (t81 * coefm[6] * t24 + t81 * coefm[7] * t18 + t107 * coefm[13] * t9) + Nm[1] * (t86 * coefm[11] * t24 + t86 * coefm[12] * t18 + t81 * coefm[4] * t9) + Nm[9] * t107 * t9 * coefm[5];
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f2_f2
// ***********************************************************************

complex<double> X_function_f2_f2 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2   = pow(B,2);
	double t3   = pow(t2,2);
	double t5   = 0.1e1 / t3 / B;
	double t9   = cos(Psi);
	double t10  = pow(t9,2);
	double t11  = sin(Psi);
	double t13  = sin(theta);
	double t14  = t10 * t11 * t13;
	double t20  = t10 * t13;
	double t39  = pow(Bm,2);
	double t40  = pow(t39,2);
	double t42  = 0.1e1 / t40 / Bm;
	//
	X = N[0] * t5 * coef[1] * t14 + N[2] / t3 * t11 * t20 * coef[2] + N[3] / t2 / B * t11 * t20 * coef[0] + N[1] * t5 * t11 * t20 * coef[3] + Nm[0] * t42 * coefm[2] * t14 + Nm[2] / t40 * t11 * t20 * coefm[1] + Nm[3] / t39 / Bm * t11 * t20 * coefm[0] + Nm[1] * t42 * t11 * t20 * coefm[3];
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f2_f3
// ***********************************************************************

complex<double> X_function_f2_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2   = pow(B,2);
	double t3   = pow(t2,2);
	double t4   = 0.1e1 / t3;
	double t7   = cos(Psi);
	double t8   = sin(Psi);
	double t9   = t7 * t8;
	double t12  = 0.1e1 / t3 / B;
	double t15  = t7 * t7;
	double t17  = sin(theta);
	double t18  = t15 * t8 * t17;
	double t22  = cos(theta);
	double t24  = t15 * t22 * t8;
	double t52  = 0.1e1 / t2 / B;
	double t79  = Bm * Bm;
	double t80  = t79 * t79;
	double t81  = 0.1e1 / t80;
	double t86  = 0.1e1 / t80 / Bm;
	double t107 = 0.1e1 / t79 / Bm;
	//
	X = N[0] * (t4 * coef[0] * t9 + t12 * coef[2] * t18 + t12 * coef[3] * t24) + N[5] * t4 * t9 * coef[1]
	  + N[7] * t4 * t9 * coef[11] + N[1] * (t12 * coef[12] * t18 + t12 * coef[13] * t24 + t4 * coef[4] * t9)
	  + N[9] * t52 * t9 * coef[5] + N[3] * (t52 * coef[6] * t18 + t52 * coef[7] * t24) 
	  + N[2] * (t4 * coef[8] * t18 + t4 * coef[9] * t24 + t52 * coef[10] * t9) 
	  + Nm[0] * (t81 * coefm[8] * t9 + t86 * coefm[1] * t24 + t86 * coefm[2] * t18) 
	  + Nm[5] * t81 * t9 * coefm[3] + Nm[7] * t81 * t9 * coefm[4] + Nm[3] * (t107 * coefm[11] * t24 + t107 * coefm[12] * t18) 
	  + Nm[9] * t107 * t9 * coefm[13] + Nm[1] * (t86 * coefm[5] * t18 + t86 * coefm[6] * t24 + t81 * coefm[0] * t9) 
	  + Nm[2] * (t81 * coefm[9] * t18 + t81 * coefm[10] * t24 + t107 * coefm[7] * t9);
	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f1
// ***********************************************************************

complex<double> X_function_f3_f1 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2  = pow(B,2);
	double t3  = pow(t2,2);
	double t4  = 1.0 / t3;
	double t7  = cos(Psi);
	double t8  = pow(t7,2);
	double t9  = sin(theta);
	double t10 = t8 * t9;
	double t13 = 1.0 / t3 / B;
	double t16 = sin(Psi);
	double t18 = t8 * t16 * t9;
	double t34 = 1.0 / t2 / B;
	double t64 = pow(Bm,2);
	double t65 = pow(t64,2);
	double t67 = 1.0 / t65 / Bm;
	double t71 = 1.0 / t65;
	double t89 = 1.0 / t64 / Bm;
	//
	X = N[0] * ( t4 * coef[0] * t10 + t13 * coef[4] * t18 )   +   N[7] * t4 * t10 * coef[8]   +   N[5] * t4 * t10 * coef[5]   +   N[9] * t34 * t10 * coef[1]
	+   N[2] * ( t34 * coef[7] * t10 + t4 * coef[3] * t18 )   +   N[3] * t16 * t34 * t10 * coef[2] 
	+   N[1] * ( t4 * coef[6] * t10 + t13 * coef[9] * t18 )   +   Nm[0] * ( t67 * coefm[1] * t18 + t71 * coefm[3] * t10 )
	+	 Nm[5] * t71 * t10 * coefm[0]   +   Nm[7] * t71 * t10 * coefm[7]   +    Nm[2] * ( t89 * coefm[9] * t10 + t71 * coefm[6] * t18 )
	+   Nm[3] * t16 * t89 * t10 * coefm[4]   +   Nm[1] * ( t71 * coefm[2] * t10 + t67 * coefm[8] * t18 )   +   Nm[9] * t89 * t10 * coefm[5] ;
	// Final Output
	return X;
}
 // ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f2
// ***********************************************************************

complex<double> X_function_f3_f2 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
{
	complex<double> X;
	//
	double t2 = B * B;
	double t3 = t2 * t2;
	double t5 = 0.1e1 / t3 / B;
	double t8 = cos(Psi);
	double t9 = t8 * t8;
	double t10 = sin(Psi);
	double t12 = sin(theta);
	double t13 = t9 * t10 * t12;
	double t15 = 0.1e1 / t3;
	double t18 = t9 * t12;
	double t29 = 0.1e1 / t2 / B;
	double t64 = Bm * Bm;
	double t65 = t64 * t64;
	double t67 = 0.1e1 / t65 / Bm;
	double t71 = 0.1e1 / t65;
	double t79 = 0.1e1 / t64 / Bm;
	//
	X = N[0] * (t5 * coef[2] * t13 + t15 * coef[3] * t18) + N[7] * t15 * t18 * coef[5] + N[2] * (t29 * coef[4] * t18 + t15 * coef[9] * t13)
	  + N[1] * (t15 * coef[1] * t18 + t5 * coef[6] * t13) + N[5] * t15 * t18 * coef[0] + N[9] * t29 * t18 * coef[8]
	  + N[3] * t10 * t29 * t18 * coef[7] + Nm[0] * (t67 * coefm[0] * t13 + t71 * coefm[2] * t18)
	  + Nm[2] * (t79 * coefm[8] * t18 + t71 * coefm[5] * t13) + Nm[5] * t71 * t18 * coefm[9]
	  + Nm[1] * (t71 * coefm[1] * t18 + t67 * coefm[7] * t13) + Nm[9] * t79 * t18 * coefm[3]
	  + Nm[7] * t71 * t18 * coefm[6] + Nm[3] * t10 * t79 * t18 * coefm[4];

	// Final Output
	return X;
}

// ***********************************************************************
//			IMPLEMENTATION OF complex<double> X_function_f3_f3
// ***********************************************************************

complex<double> X_function_f3_f3 (double theta, double Psi, double B, double Bm, complex<double> coef[], complex<double> coefm[], complex<double> N[], complex<double> Nm[])
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
	double t16 = sin(Psi);
	double t17 = t10 * t16;
	double t24 = sin(theta);
	double t25 = t8 * t24;
	double t29 = t7 * t16;
	double t35 = 0.1e1 / t2 / B;
	double t79 = Bm * Bm;
	double t80 = t79 * t79;
	double t81 = 0.1e1 / t80;
	double t86 = 0.1e1 / t80 / Bm;
	double t112 = 0.1e1 / t79 / Bm;
	//
	X = N[0] * (t4 * coef[0] * t10 + t13 * coef[3] * t17) + N[5] * (t4 * coef[1] * t25 + t4 * coef[2] * t29)
	  + N[9] * (t35 * coef[5] * t25 + t35 * coef[6] * t29) + N[3] * t9 * t16 * t35 * t8 * coef[7] + N[2] * (t35 * coef[9] * t10 + t4 * coef[8] * t17)
	  + N[1] * (t4 * coef[4] * t10 + t13 * coef[12] * t17) + N[7] * (t4 * coef[10] * t25 + t4 * coef[11] * t29)
	  + Nm[0] * (t81 * coefm[4] * t10 + t86 * coefm[7] * t17) + Nm[5] * (t81 * coefm[5] * t25 + t81 * coefm[6] * t29)
	  + Nm[7] * (t81 * coefm[1] * t25 + t81 * coefm[2] * t29) + Nm[9] * (t112 * coefm[9] * t25 + t112 * coefm[10] * t29)
	  + Nm[3] * t9 * t16 * t112 * t8 * coefm[11] + Nm[2] * (t112 * coefm[0] * t10 + t81 * coefm[12] * t17)
	  + Nm[1] * (t81 * coefm[8] * t10 + t86 * coefm[3] * t17);
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
	complex<double> a  = Iunit * ko * B;
	complex<double> am = Iunit * ko * Bm;
	//
	if (Psi >= PsiB && Psi >= PsiA)
	{
		D  = sqrt(3.0) / sin(Psi);
		X2(D, a, N);
		X2 (D, am, Nm);
	}
	else if (theta >= M_PI_2)
	{
		D  = sqrt(3.0) / (cos(Psi) * tPsiA);
		X2 (D, a, N);
		X2 (D, am, Nm);
	}
	else if (Psi >= PsiA)
	{
		D1 = 2.0 * sqrt(3.0) / (cos(Psi) * ( tPsiB + tan(Psi) )  );
		D2 = sin(Psi) / sqrt(3.0);
		X1(Psi, D1, D2, tPsiB, a, N);
		X1(Psi, D1, D2, tPsiB, am, Nm);
	}
	else
	{
		D1 = sqrt(3.0) / (cos(Psi) * sin(theta) );
		D2 = ( cos(Psi) * tPsiA ) / sqrt(3.0);
		X1(Psi, D1, D2, tPsiB, a, N);
		X1(Psi, D1, D2, tPsiB, am, Nm);
	}

}

// ***********************************************************************
//			IMPLEMENTATION OF void void X1
// ***********************************************************************

void X1 ( double psi, complex<double> D1, complex<double> D2, double tpsiB, complex<double> a, complex<double> N[] )
{
	complex<double> D3, aD3, expaD3, aD1, expaD1, D2a, H, T11, T21, T31, T41, T12, T22, T32, T42;
	complex<double>	N21, N31, N81, N91, N101, N111, N121, N41, N51, N22, N32, N82, N92, N102, N112, N122, N42, N52;
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
	N[6]   = 1.0 / 3.0;
	N[7]   = N81+N82;
	N[8]   = N91+N92;
	N[9]   = N101+N102;
	N[10]  = N111+N112;
	N[11]  = N121+N122;
}

// ***********************************************************************
//			IMPLEMENTATION OF void void X2
// ***********************************************************************

void X2 ( complex<double> D, complex<double> a, complex<double> N[] )
{
	complex<double> aD, expaD, T1, T2, T3, T4;
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
	N[6]   = 1.0 / 3.0;
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


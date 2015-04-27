




#include <stdlib.h> 
#include <mpi.h>
#include <math.h>
#include <omp.h>
#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void X_boundary_condition
// ============================================================================ //
(
int myid,

double deltaT,

double e,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
)
// ============================================================================ //
{

#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

	double L1, L2, L3, L4, L5;

	double dp_dx, du_dx, dv_dx, dw_dx, dT_dx;

	double U_p, C_p;

	double rho, U, V, W, VV, P, C, T;

	double _rho, _U, _V, _W, _VV, _P, _C, _T;

	double rho_, U_, V_, W_, VV_, P_, C_, T_;

	double rhoold, Uold, Vold, Wold, VVold, Pold, Cold, Told; 

	double beta;

	double deltaTau = deltaT/100.;

	/*
	if (myid == 0) {

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

				istart=3;

				U1_[istart-1][j][k] = U1_[istart][j][k];
				U2_[istart-1][j][k] = U2_[istart][j][k];
				U3_[istart-1][j][k] = U3_[istart][j][k];
				U4_[istart-1][j][k] = U4_[istart][j][k];
				U5_[istart-1][j][k] = U5_[istart][j][k];

			}
		}


	}
	
	
	if (myid == nproc-1) {

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {


				iend = gend[myid];

				U1_[iend+1][j][k] = U1_[iend][j][k];
				U2_[iend+1][j][k] = U2_[iend][j][k];
				U3_[iend+1][j][k] = U3_[iend][j][k];
				U4_[iend+1][j][k] = U4_[iend][j][k];
				U5_[iend+1][j][k] = U5_[iend][j][k];

			}
		}

	}
	*/

//#pragma omp barrier
	
	

	if (myid == 0) {

		istart=3;

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

				rhoold = U1q[istart-1][j][k]*J[istart-1][j][k];
				Uold = U2q[istart-1][j][k]/U1q[istart-1][j][k];
				Vold = U3q[istart-1][j][k]/U1q[istart-1][j][k];
				Wold = U4q[istart-1][j][k]/U1q[istart-1][j][k];     
				VVold = Uold*Uold+Vold*Vold+Wold*Wold;
				Pold = (U5q[istart-1][j][k]*J[istart-1][j][k]-0.5*rhoold*VVold)*(K-1);
				Cold = sqrt(K*Pold/rhoold);
				Told = Pold/rhoold;

				rho_ = U1_[istart][j][k]*J[istart][j][k];
				U_ = U2_[istart][j][k]/U1_[istart][j][k];
				V_ = U3_[istart][j][k]/U1_[istart][j][k];
				W_ = U4_[istart][j][k]/U1_[istart][j][k];
				VV_ = U_*U_+V_*V_+W_*W_;
				P_ = (U5_[istart][j][k]*J[istart][j][k]-0.5*rho_*VV_)*(K-1);
				C_ = sqrt(K*P_/rho_);
				T_ = P_/rho_;

				rho = U1_[istart-1][j][k]*J[istart-1][j][k];
				U = U2_[istart-1][j][k]/U1_[istart-1][j][k];
				V = U3_[istart-1][j][k]/U1_[istart-1][j][k];
				W = U4_[istart-1][j][k]/U1_[istart-1][j][k];
				VV = U*U+V*V+W*W;
				P = (U5_[istart-1][j][k]*J[istart-1][j][k]-0.5*rho*VV)*(K-1);
				C = sqrt(K*P/rho);
				T = P/rho;
				
				beta = max(VV/C/C,e);

				U_p = 0.5*(beta+1)*U;

				C_p = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dx = (P_-P)/deltaXI;

				du_dx = (U_-U)/deltaXI;

				dv_dx = (V_-V)/deltaXI;

				dw_dx = (W_-W)/deltaXI;

				dT_dx = (T_-T)/deltaXI;

				L5 = (U_p-C_p)*(dp_dx-rho*(U_p+C_p-U)*du_dx);

				//L4 =  L5*(U_p-C_p-U)/(U_p+C_p-U);
				
				L4 = 0.25*(U_p-C_p-U)/(X_out*deltaXI)*(P-101300);
				
				if (U <= 0) {

					//L1 = U*(dT_dx+1/rho*dp_dx*(1-K)/K);
					L1 = U*dT_dx+1/rho*dp_dx*(1-K)/K;

					L2 = U*dw_dx;

					L3 = -U*dv_dx;

				}

				else {

					//L1 = L2 = L3 = 0;

					L1 = 0;
					
					L2 = 0.25*(U_p-C_p-U)/(X_out*deltaXI)*(W-0.);
					
					L3 = 0.25*(U_p-C_p-U)/(X_out*deltaXI)*(V-0.);
					
				}

				U = Uold-0.5/(rho*C_p)*(L4-L5)*deltaTau;

				P = Pold-0.5/(C_p)*(L4*(U_p+C_p-U)-L5*(U_p-C_p-U))*deltaTau;
				
				//P = Pold-rho*(U_p-C_p-U)*(U-Uold)-L4;

				V = Vold+L3*deltaTau;

				W = Wold-L2*deltaTau;

				T = Told-L1*deltaTau+1/rho*(K-1)/K*(P-Pold);

				rho = P/T;


				U1_[istart-1][j][k] = rho/J[istart-1][j][k];
				U2_[istart-1][j][k] = rho*U/J[istart-1][j][k];
				U3_[istart-1][j][k] = rho*V/J[istart-1][j][k];
				U4_[istart-1][j][k] = rho*W/J[istart-1][j][k];
				U5_[istart-1][j][k] = (P/(K-1)+0.5*rho*(U*U+V*V+W*W))/J[istart-1][j][k];

				U1q[istart-1][j][k] = U1_[istart-1][j][k];
				U2q[istart-1][j][k] = U2_[istart-1][j][k];
				U3q[istart-1][j][k] = U3_[istart-1][j][k];
				U4q[istart-1][j][k] = U4_[istart-1][j][k];
				U5q[istart-1][j][k] = U5_[istart-1][j][k];

			}
		}


	}


	if (myid == np-1) {

		iend = gend[myid];

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

				/* outlet */
				
				rhoold = U1q[iend+1][j][k]*J[iend+1][j][k];
				Uold = U2q[iend+1][j][k]/U1q[iend+1][j][k];
				Vold = U3q[iend+1][j][k]/U1q[iend+1][j][k];
				Wold = U4q[iend+1][j][k]/U1q[iend+1][j][k];     
				VVold = Uold*Uold+Vold*Vold+Wold*Wold;
				Pold = (U5q[iend+1][j][k]*J[iend+1][j][k]-0.5*rhoold*VVold)*(K-1);
				Cold = sqrt(K*Pold/rhoold);
				Told = Pold/rhoold;

				_rho = U1_[iend][j][k]*J[iend][j][k];
				_U = U2_[iend][j][k]/U1_[iend][j][k];
				_V = U3_[iend][j][k]/U1_[iend][j][k];
				_W = U4_[iend][j][k]/U1_[iend][j][k];
				_VV = _U*_U+_V*_V+_W*_W;
				_P = (U5_[iend][j][k]*J[iend][j][k]-0.5*_rho*_VV)*(K-1);
				_C = sqrt(K*_P/_rho);
				_T = _P/_rho;

				rho = U1_[iend+1][j][k]*J[iend+1][j][k];
				U = U2_[iend+1][j][k]/U1_[iend+1][j][k];
				V = U3_[iend+1][j][k]/U1_[iend+1][j][k];
				W = U4_[iend+1][j][k]/U1_[iend+1][j][k];
				VV = U*U+V*V+W*W;
				P = (U5_[iend+1][j][k]*J[iend+1][j][k]-0.5*rho*VV)*(K-1);
				C = sqrt(K*P/rho);
				T = P/rho;

				beta = max(VV/C/C,e);

				U_p = 0.5*(beta+1)*U;

				C_p = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dx = (P-_P)/deltaXI;

				du_dx = (U-_U)/deltaXI;

				dv_dx = (V-_V)/deltaXI;

				dw_dx = (W-_W)/deltaXI;

				dT_dx = (T-_T)/deltaXI;

				L4 = (U_p+C_p)*(dp_dx-rho*(U_p-C_p-U)*du_dx);

				//L5 = L4*(U_p+C_p-U)/(U_p-C_p-U);
				
				L5 = 0.25*(U_p+C_p-U)/(X_out*deltaXI)*(P-101300);
				
				if (U > 0) {

					//L1 = U*(dT_dx+1/rho*dp_dx*(1-K)/K);
					
					L1 = U*dT_dx+1/rho*dp_dx*(1-K)/K;

					L2 = U*dw_dx;

					L3 = -U*dv_dx;
					
				}

				else {

					//L1 = L2 = L3 = 0;
					
					L1 = 0;
					
					L2 = 0.25*(U_p-C_p-U)/(X_out*deltaXI)*(W-0.);
					
					L3 = 0.25*(U_p-C_p-U)/(X_out*deltaXI)*(V-0.);
					
				}
				
				
				U = Uold-0.5/(rho*C_p)*(L4-L5)*deltaTau;

				P = Pold-0.5/(C_p)*(L4*(U_p+C_p-U)-L5*(U_p-C_p-U))*deltaTau;
				
				//P = Pold-rho*(U_p+C_p-U)*(U-Uold)-L5;
				
				V = Vold+L3*deltaTau;

				W = Wold-L2*deltaTau;

				T = Told-L1*deltaTau+1/rho*(K-1)/K*(P-Pold);

				rho = P/T;

				U1_[iend+1][j][k] = rho/J[iend+1][j][k];
				U2_[iend+1][j][k] = rho*U/J[iend+1][j][k];
				U3_[iend+1][j][k] = rho*V/J[iend+1][j][k];
				U4_[iend+1][j][k] = rho*W/J[iend+1][j][k];
				U5_[iend+1][j][k] = (P/(K-1)+0.5*rho*(U*U+V*V+W*W))/J[iend+1][j][k];

				U1q[iend+1][j][k] = U1_[iend+1][j][k];
				U2q[iend+1][j][k] = U2_[iend+1][j][k];
				U3q[iend+1][j][k] = U3_[iend+1][j][k];
				U4q[iend+1][j][k] = U4_[iend+1][j][k];
				U5q[iend+1][j][k] = U5_[iend+1][j][k];    


			}
		}

	}


}

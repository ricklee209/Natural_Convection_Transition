




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

double (*U1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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

	double rho, U, V, W, VV, P, C, T, H;

	double _rho, _U, _V, _W, _VV, _P, _C, _T;

	double rho_, U_, V_, W_, VV_, P_, C_, T_;

	double rhoold, Uold, Vold, Wold, VVold, Pold, Cold, Told; 

	double beta;

	double deltaTau = deltaT/200.;
	
	double d11,d12,d13,d14,d15,
	       d21,d22,d23,d24,d25,
	       d31,d32,d33,d34,d35,
	       d41,d42,d43,d44,d45,
	       d51,d52,d53,d54,d55;
		   
	double KAL1, KAL2, KAL3, KAL4, KAL5;
		   
	double temp,temp2,temp3;
	
	
	if (myid == 0) {
		
		istart=3;

		for (j = 1; j <= nyy; j++) {
			for (k = 2; k < nzz; k++) {

				U1_[istart-1][j][k] = U1_[istart][j][k];
				U2_[istart-1][j][k] = -U2_[istart][j][k];
				U3_[istart-1][j][k] = -U3_[istart][j][k];
				U4_[istart-1][j][k] = -U4_[istart][j][k];
				U5_[istart-1][j][k] = U5_[istart][j][k];
				
			}
		}
		
	}
	
	

	
	if (myid == nproc-1) {
		
		iend = gend[myid];

		for (j = 1; j <= nyy; j++) {
			for (k = 2; k < nzz; k++) {

				U1_[iend+1][j][k] = U1_[iend][j][k];
				U2_[iend+1][j][k] = -U2_[iend][j][k];
				U3_[iend+1][j][k] = -U3_[iend][j][k];
				U4_[iend+1][j][k] = -U4_[iend][j][k];
				U5_[iend+1][j][k] = U5_[iend][j][k];

			}
		}

	}
	
	
	

}






#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;


#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Y_boundary_condition
(
 // ============================================================================ //
 int myid,

 double deltaT,

 double e,

 double heat_flux,

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

 double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

 double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]

 // ============================================================================ //
 )

{

#include "ijk.h"
#include "Viscous_terms.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	 double rho,U,V,W,VV,P,C,T,h,H;
	 double mu_E;
	 double Tb,TTu,TTb;


	 double L1, L2, L3, L4, L5;

	 double dp_dy, du_dy, dv_dy, dw_dy, dT_dy;

	 double V_p, C_p;

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



	 istart = 3;
	 iend = gend[myid];

	 //#pragma omp parallel for private(k,rho,U,V,W,VV,P,temp,T)
	 for (i = istart; i <= iend; i++) {


		  for (k = 2; k <= nz; k++) {

			   U1_[i][nyy][k] = U1_[i][ny][k];
			   U2_[i][nyy][k] = -U2_[i][ny][k];
			   U3_[i][nyy][k] = -U3_[i][ny][k];
			   U4_[i][nyy][k] = -U4_[i][ny][k];
			   U5_[i][nyy][k] = U5_[i][ny][k];
		  }


		  for (k = 2; k <= nz; k++) {

			   U1_[i][1][k] = U1_[i][2][k];
			   U2_[i][1][k] = -U2_[i][2][k];
			   U3_[i][1][k] = -U3_[i][2][k];
			   U4_[i][1][k] = -U4_[i][2][k];
			   U5_[i][1][k] = U5_[i][2][k];

		  }	


		  if (i+gstart[myid] > nx_inlet+2 && i+gstart[myid] < (X_out-nx_outlet+3)) {

			   for (k = 2; k <= nz; k++) {

					rho = U1_[i][2][k]*J[i][2][k];
					U = U2_[i][2][k]/U1_[i][2][k];
					V = U3_[i][2][k]/U1_[i][2][k];
					W = U4_[i][2][k]/U1_[i][2][k];     
					VV = U*U+V*V+W*W;
					P = (U5_[i][2][k]*J[i][2][k]-0.5*rho*VV)*(K-1);
					temp = P/rho/R;

					mu_E = mu_L*pow((temp/298.0592),1.5)*(298.0592+110.)/(temp+110.);

					lambda_L = mu_E*Cv*K/Pr_L;

					T = heat_flux/lambda_L/etdy[i][1][k]*deltaET+temp;

					rho = P/R/T;

					U1_[i][1][k] = rho/J[i][1][k];
					U2_[i][1][k] = -rho*U/J[i][1][k];
					U3_[i][1][k] = -rho*V/J[i][1][k];
					U4_[i][1][k] = -rho*W/J[i][1][k];
					U5_[i][1][k] = (P/(K-1)+0.5*rho*VV)/J[i][1][k];

			   }


		  }

	 }





}

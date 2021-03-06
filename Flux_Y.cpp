



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

void Flux_Y
(
// ============================================================================ //
int myid,

double Ep,

double e,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*J_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"


	double rho,U,V,W,VV,P,C,T,h,H;
	double u,v,w;
	double temp,temp1,temp2,temp3,temp4,temp5;
	double beta,S,_S_;

	double m11,m12,m13,m14,m15,m22,m23,m24,m25,m33,m34,m35,m44,m45,m55;
	double thedac1,thedac2,thedac3,thedad1,thedad2,thedad3;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2;
	double d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35;
	double d41,d42,d43,d44,d45,d51,d52,d53,d54,d55, Fav1,Fav2,Fav3,Fav4,Fav5;

/**** MUSCL 5th-order ****/
// ============================================================================ //
	
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;						  ////	
//// ============================================ ////
		iend = gend[myid];						  ////
//// ============================================ ////
	
	

	for (i = istart; i <= iend; i++) {
#pragma omp parallel for private(k,temp)
		for ( j = 4; j <= ny-2; j++) {
			for (k = 2; k <= nz; k++) {
			
				temp = 1./(2/J[i][j-2][k]-13/J[i][j-1][k]+47/J[i][j][k]+27/J[i][j+1][k]-3/J[i][j+2][k]);

				ML1[i-1][j][k-1] = temp*(2*U1_[i][j-2][k]-13*U1_[i][j-1][k]+47*U1_[i][j][k]+27*U1_[i][j+1][k]-3*U1_[i][j+2][k]);
				ML2[i-1][j][k-1] = temp*(2*U2_[i][j-2][k]-13*U2_[i][j-1][k]+47*U2_[i][j][k]+27*U2_[i][j+1][k]-3*U2_[i][j+2][k]);
				ML3[i-1][j][k-1] = temp*(2*U3_[i][j-2][k]-13*U3_[i][j-1][k]+47*U3_[i][j][k]+27*U3_[i][j+1][k]-3*U3_[i][j+2][k]);
				ML4[i-1][j][k-1] = temp*(2*U4_[i][j-2][k]-13*U4_[i][j-1][k]+47*U4_[i][j][k]+27*U4_[i][j+1][k]-3*U4_[i][j+2][k]);
				ML5[i-1][j][k-1] = temp*(2*U5_[i][j-2][k]-13*U5_[i][j-1][k]+47*U5_[i][j][k]+27*U5_[i][j+1][k]-3*U5_[i][j+2][k]);

				temp = 1./(-3/J[i][j-2][k]+27/J[i][j-1][k]+47/J[i][j][k]-13/J[i][j+1][k]+2/J[i][j+2][k]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i][j-2][k]+27*U1_[i][j-1][k]+47*U1_[i][j][k]-13*U1_[i][j+1][k]+2*U1_[i][j+2][k]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i][j-2][k]+27*U2_[i][j-1][k]+47*U2_[i][j][k]-13*U2_[i][j+1][k]+2*U2_[i][j+2][k]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i][j-2][k]+27*U3_[i][j-1][k]+47*U3_[i][j][k]-13*U3_[i][j+1][k]+2*U3_[i][j+2][k]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i][j-2][k]+27*U4_[i][j-1][k]+47*U4_[i][j][k]-13*U4_[i][j+1][k]+2*U4_[i][j+2][k]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i][j-2][k]+27*U5_[i][j-1][k]+47*U5_[i][j][k]-13*U5_[i][j+1][k]+2*U5_[i][j+2][k]);
					
				
			}
		}
	}
	

	
	
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
		
			
		
			j = 1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// ====================================================================================//
		
			
			j = 2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// ====================================================================================//		
			j = 3; 
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//	
		
			j = ny-1;
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//



			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
// ====================================================================================//

			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
			MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
				
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	for (i = istart; i <= iend; i++) {
	
		if (i+gstart[myid] > nx_inlet+2 && i+gstart[myid] < (X_out-nx_outlet)+3) {
		
			for (k = 2; k <= nz; k++) {
		
			
			j = Y_out-ny_abs+1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// ====================================================================================//
		
			
			j = Y_out-ny_abs+2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// ====================================================================================//		
			j = Y_out-ny_abs+3;
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//	
				
		}
		
		
		
		
		
			for (k = 2; k <= nz; k++) {
			
	// ------------------------------------------------------------------------------------ //			
	// ---------------------------------- adiabatic wall ---------------------------------- //

	
	// ====================================================================================//		
				j = Y_out-ny_abs; 
				
				temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

				ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
				ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
				ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
				ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
				ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
							
				temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
							
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
				
// ====================================================================================//	
							
				
				j = Y_out-ny_abs+1;
				
				temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
						
				ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
				ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
				ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
				ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
				ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
				
				
				temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
				
				MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
				MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
				MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
				MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
				MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
				
				
	// ====================================================================================//	
			
				j = Y_out-ny_abs+1;
				
				temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
				MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
				MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
				MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
				MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
				MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
				
				
			}
			
		}
	}
		
	
	
	



//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

	/*---Y fluxes---*/
	for (i = istart; i <= iend; i++) {

	#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ETx,ETy,ETz,\
	_rho,_u,_v,_w,_U,_V,_W,__V,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,V__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_V_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)

		for (j = 1; j < nyy; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_v[i][j][k];
				xiy=xidy_v[i][j][k];
				xiz=xidz_v[i][j][k];
				etx=etdx_v[i][j][k];
				ety=etdy_v[i][j][k];
				etz=etdz_v[i][j][k];          
				ztx=ztdx_v[i][j][k];
				zty=ztdy_v[i][j][k];
				ztz=ztdz_v[i][j][k];
				ETx=etx/(sqrt(etx*etx+ety*ety+etz*etz));
				ETy=ety/(sqrt(etx*etx+ety*ety+etz*etz));
				ETz=etz/(sqrt(etx*etx+ety*ety+etz*etz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				_U = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;

				_W = ztx*_u+zty*_v+ztz*_w;


				__V = ETx*_u+ETy*_v+ETz*_w;


				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_T = _P/_rho;
				_C = K*_P/_rho;
				_H = 0.5*_VV+_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				V__ = ETx*u_+ETy*v_+ETz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ =K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);

				/* flux varibale */
				rho = 0.5*(_rho+rho_);
				u = 0.5*(_u+u_);
				v = 0.5*(_v+v_);
				w = 0.5*(_w+w_);

				U = 0.5*(_U+U_);
				V = 0.5*(_V+V_);
				W = 0.5*(_W+W_); 

				_V_ = 0.5*(__V+V__);

				VV = u*u+v*v+w*w;
				H = 0.5*(_H+H_);
				C = (H-0.5*VV)*(K-1);
				P = rho*C/K;
				T = P/rho;

				

				/* jump dU */
				dU1 = P_-_P;
				dU2 = u_-_u;
				dU3 = v_-_v;
				dU4 = w_-_w;
				dU5 = T_-_T;

				
				/* preconditioning */
				beta = max(VV/C,e);

				/*V=etx*u+ety*v+etz*w;*/
				//S=sqrt(V*V*(beta-1)+4*beta*C*C*(etx*etx+ety*ety+etz*etz));
				S=sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*(etx*etx+ety*ety+etz*etz));
				/*_V_=ETx*u+ETy*v+ETz*w;*/
				_S_=sqrt(_V_*_V_*(beta-1)*(beta-1)+4*beta*C*(ETx*ETx+ETy*ETy+ETz*ETz));

				temp = 4*K*_S_*T*beta;
				temp1 = S+V+V*beta;
				temp2 = -S+V+V*beta;
				temp3 = _S_-_V_+_V_*beta;
				temp4 = _S_+_V_-_V_*beta;
				temp5 = _S_*_S_-(beta-1)*(beta-1)*_V_*_V_;
				tempu1 = 2*T*ETx*K*beta+u*_S_+u*(beta-1)*_V_;
				tempu2 = 2*T*ETx*K*beta-u*_S_+u*(beta-1)*_V_;
				tempv1 = 2*T*ETy*K*beta+v*_S_+v*(beta-1)*_V_;
				tempv2 = 2*T*ETy*K*beta-v*_S_+v*(beta-1)*_V_;
				tempw1 = 2*T*ETz*K*beta+w*_S_+w*(beta-1)*_V_;
				tempw2 = 2*T*ETz*K*beta-w*_S_+w*(beta-1)*_V_;
				tempuvw = 2*T*K*beta*(u*ETx+v*ETy+w*ETz);
				temph1 = tempuvw+H*_S_+H*_V_*beta-H*_V_; 
				temph2 = tempuvw-H*_S_+H*_V_*beta-H*_V_;


				d11 = 1/temp*(4*(K-1)*_S_*beta*fabs(V)+temp3*fabs(temp1)+temp4*fabs(temp2));
				d12 = -1/(2*temp)*(ETx*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d13 = -1/(2*temp)*(ETy*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d14 = -1/(2*temp)*(ETz*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d15 = -rho*fabs(V)/T;

				d21 = 1/temp*(u*_S_*(4*(K-1)*beta*fabs(V)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ETx*K*beta+u*_V_*(beta-1)));
				d22 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ETy*ETy+ETz*ETz)*fabs(V)+ETx*(temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1))));
				d23 = 1/(2*temp)*(rho*ETy*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d24 = 1/(2*temp)*(rho*ETz*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d25 = -u*rho*fabs(V)/T;

				d31 = 1/temp*(v*_S_*(4*(K-1)*beta*fabs(V)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ETy*K*beta+v*_V_*(beta-1)));
				d32 = 1/(2*temp)*(rho*ETx*(-8*T*K*_S_*beta*ETy*fabs(V)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d33 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ETx*ETx+ETz*ETz)*fabs(V)+ETy*(temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1))));
				d34 = 1/(2*temp)*(rho*ETz*(-8*T*K*_S_*beta*ETy*fabs(V)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d35 = -v*rho*fabs(V)/T;

				d41 = 1/temp*(w*_S_*(4*(K-1)*beta*fabs(V)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ETz*K*beta+w*_V_*(beta-1)));
				d42 = 1/(2*temp)*(rho*ETx*(-8*T*K*_S_*beta*ETz*fabs(V)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d43 = 1/(2*temp)*(rho*ETy*(-8*T*K*_S_*beta*ETz*fabs(V)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d44 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ETx*ETx+ETy*ETy)*fabs(V)+ETz*(temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1))));
				d45 = -w*rho*fabs(V)/T;

				d51 = 1/temp*(_S_*(4*beta*(H*K-H-T*K)*fabs(V)+H*(fabs(temp2)+fabs(temp1)))-(fabs(temp2)-fabs(temp1))*(tempuvw+H*_V_*beta-H*_V_));
				d52 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-v*ETx*ETy-w*ETx*ETz+u*(ETy*ETy+ETz*ETz)*fabs(V))+ETx*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d53 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*ETx*ETy-w*ETy*ETz+v*(ETx*ETx+ETz*ETz)*fabs(V))+ETy*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d54 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*ETx*ETz-v*ETy*ETz+w*(ETx*ETx+ETy*ETy)*fabs(V))+ETz*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d55 = (-H/T+K/(K-1))*rho*fabs(V);

				/* artificial viscosity */
				Fav1 = d11*dU1+d12*dU2+d13*dU3+d14*dU4+d15*dU5;
				Fav2 = d21*dU1+d22*dU2+d23*dU3+d24*dU4+d25*dU5;
				Fav3 = d31*dU1+d32*dU2+d33*dU3+d34*dU4+d35*dU5;
				Fav4 = d41*dU1+d42*dU2+d43*dU3+d44*dU4+d45*dU5;
				Fav5 = d51*dU1+d52*dU2+d53*dU3+d54*dU4+d55*dU5;

				/* inviscid fluxes */
				
				inFy1[i][j][k] = 0.5*((_rho*_V+rho_*V_)-Ep*Fav1)/J_v[i][j][k];
				inFy2[i][j][k] = 0.5*((_rho*_u*_V+rho_*u_*V_+etx*(_P+P_))-Ep*Fav2)/J_v[i][j][k];
				inFy3[i][j][k] = 0.5*((_rho*_v*_V+rho_*v_*V_+ety*(_P+P_))-Ep*Fav3)/J_v[i][j][k];
				inFy4[i][j][k] = 0.5*((_rho*_w*_V+rho_*w_*V_+etz*(_P+P_))-Ep*Fav4)/J_v[i][j][k];
				inFy5[i][j][k] = 0.5*((_V*(3.5*_P+0.5*_rho*_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5)/J_v[i][j][k];
				
				
				
				if (j==1 | j==ny) {

					inFy1[i][j][k] = 0;
					inFy2[i][j][k] = 0;
					inFy3[i][j][k] = 0.5*(ety*(_P+P_)-Ep*Fav3)/J_v[i][j][k];
					inFy4[i][j][k] = 0;
					inFy5[i][j][k] = 0;

				}
				
				
			}
		}
	}
	

}
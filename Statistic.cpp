




#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

void Statistic
(
// ============================================================================ //
int myid,
int step,
int iteration,
int statistic_step,

double obs1,
double obs2,
double obs3,
double obs4,
double obs5,

double e,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

	
#include "ijk.h"
#include "Viscous_terms.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	istart1 = gstart[myid];

	char LESdata[100];
	
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


// ============================================================================================================= //
	
	static double 
		Um[X_m][Y_m],Vm[X_m][Y_m],Wm[X_m][Y_m],Tm[X_m][Y_m],
		UUm[X_m][Y_m],VVm[X_m][Y_m],WWm[X_m][Y_m],TTm[X_m][Y_m],UVm[X_m][Y_m],
		UTm[X_m][Y_m],VTm[X_m][Y_m],
		ULm[X_m][Y_m],URm[X_m][Y_m],
		VLm[X_m][Y_m],VRm[X_m][Y_m],
		rhoRm[X_m][Y_m],rhoLm[X_m][Y_m],
		RUVm[X_m][Y_m],SXYm[X_m][Y_m],Roe_Dm[X_m][Y_m],dRLm[X_m][Y_m];

	// ============================ //

	
	int x_np = gcount[myid]+6; 
	
	// double (*MUL)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];
	// double (*MUR)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];

	// double (*MVL)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];
	// double (*MVR)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];

	// double (*MrhoL)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];
	// double (*MrhoR)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];

	// double (*MRUV)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];
	// double (*MSXY)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];

	// double (*MRoe_D)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];
	// double (*MdRL)[Y_m][Z_m] = new double[x_np][Y_m][Z_m];

// ============================================================================================================= //

	
/**** MUSCL 5th-order ****/
// ============================================================================ //
	
// //// ============================================ ////
		// if (myid ==0) istart = 2;		          ////	
		// else istart = 3;						  ////	
// //// ============================================ ////
		// iend = gend[myid];						  ////
// //// ============================================ ////
	
	

	// for (i = istart; i <= iend; i++) {
// #pragma omp parallel for private(k,temp)
		// for ( j = 4; j <= ny-2; j++) {
			// for (k = 2; k <= nz; k++) {
			
				// temp = 1./(2/J[i][j-2][k]-13/J[i][j-1][k]+47/J[i][j][k]+27/J[i][j+1][k]-3/J[i][j+2][k]);

				// ML1[i-1][j][k-1] = temp*(2*U1_[i][j-2][k]-13*U1_[i][j-1][k]+47*U1_[i][j][k]+27*U1_[i][j+1][k]-3*U1_[i][j+2][k]);
				// ML2[i-1][j][k-1] = temp*(2*U2_[i][j-2][k]-13*U2_[i][j-1][k]+47*U2_[i][j][k]+27*U2_[i][j+1][k]-3*U2_[i][j+2][k]);
				// ML3[i-1][j][k-1] = temp*(2*U3_[i][j-2][k]-13*U3_[i][j-1][k]+47*U3_[i][j][k]+27*U3_[i][j+1][k]-3*U3_[i][j+2][k]);
				// ML4[i-1][j][k-1] = temp*(2*U4_[i][j-2][k]-13*U4_[i][j-1][k]+47*U4_[i][j][k]+27*U4_[i][j+1][k]-3*U4_[i][j+2][k]);
				// ML5[i-1][j][k-1] = temp*(2*U5_[i][j-2][k]-13*U5_[i][j-1][k]+47*U5_[i][j][k]+27*U5_[i][j+1][k]-3*U5_[i][j+2][k]);

				// temp = 1./(-3/J[i][j-2][k]+27/J[i][j-1][k]+47/J[i][j][k]-13/J[i][j+1][k]+2/J[i][j+2][k]);
						
				// MR1[i-1][j-1][k-1] = temp*(-3*U1_[i][j-2][k]+27*U1_[i][j-1][k]+47*U1_[i][j][k]-13*U1_[i][j+1][k]+2*U1_[i][j+2][k]);
				// MR2[i-1][j-1][k-1] = temp*(-3*U2_[i][j-2][k]+27*U2_[i][j-1][k]+47*U2_[i][j][k]-13*U2_[i][j+1][k]+2*U2_[i][j+2][k]);
				// MR3[i-1][j-1][k-1] = temp*(-3*U3_[i][j-2][k]+27*U3_[i][j-1][k]+47*U3_[i][j][k]-13*U3_[i][j+1][k]+2*U3_[i][j+2][k]);
				// MR4[i-1][j-1][k-1] = temp*(-3*U4_[i][j-2][k]+27*U4_[i][j-1][k]+47*U4_[i][j][k]-13*U4_[i][j+1][k]+2*U4_[i][j+2][k]);
				// MR5[i-1][j-1][k-1] = temp*(-3*U5_[i][j-2][k]+27*U5_[i][j-1][k]+47*U5_[i][j][k]-13*U5_[i][j+1][k]+2*U5_[i][j+2][k]);
					
				
			// }
		// }
	// }
	
	
	// for (i = istart; i <= iend; i++) {
		// for (k = 2; k <= nz; k++) {
		
			// j = 1;
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			// ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			// ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			// ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			// ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			// ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// // ====================================================================================//
		
			
			// j = 2;
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			// ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			// ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			// ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			// ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			// ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			// MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			// MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			// MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			// MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			// MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// // ====================================================================================//		
			// j = 3; 
			
			// temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			// ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			// ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			// ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			// ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			// ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			// temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			// MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			// MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			// MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			// MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			// MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// // ====================================================================================//	
		
			// j = ny-1;
			
			// temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			// ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			// ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			// ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			// ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			// ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			// temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			// MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			// MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			// MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			// MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			// MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// // ====================================================================================//



			// j = ny;
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
					
			// ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			// ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			// ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			// ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			// ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
			
			// MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			// MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			// MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			// MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			// MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
// // ====================================================================================//

			// j = ny;
			
			// temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
			// MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			// MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			// MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			// MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			// MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
				
		// }
	// }
	
	
	
// //// ============================================ ////
			// istart = 2;							  ////	
// //// ============================================ ////
			// iend = gend[myid]-1;	     		  ////
// //// ============================================ ////

	// #pragma omp parallel for private(\
	// xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ETx,ETy,ETz,\
	// _rho,_u,_v,_w,_U,_V,_W,__V,_VV,_P,_T,_C,_H,\
	// rho_,u_,v_,w_,U_,V_,W_,V__,VV_,P_,T_,C_,H_,\
	// rho,u,v,w,U,V,W,_V_,VV,H,C,P,T,\
	// dU1,dU2,dU3,dU4,dU5,\
	// beta,S,_S_,\
	// temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	// d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	// d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	// Fav1,Fav2,Fav3,Fav4,Fav5,\
	// j,k\
	// )

	// /*---Y fluxes---*/
	// for (i = istart; i <= iend; i++) {
		// for (j = 1; j < nyy; j++) {
			// for (k = 1; k < nz; k++) {

				// xix=xidx_v[i][j][k];
				// xiy=xidy_v[i][j][k];
				// xiz=xidz_v[i][j][k];
				// etx=etdx_v[i][j][k];
				// ety=etdy_v[i][j][k];
				// etz=etdz_v[i][j][k];          
				// ztx=ztdx_v[i][j][k];
				// zty=ztdy_v[i][j][k];
				// ztz=ztdz_v[i][j][k];
				// ETx=etx/(sqrt(etx*etx+ety*ety+etz*etz));
				// ETy=ety/(sqrt(etx*etx+ety*ety+etz*etz));
				// ETz=etz/(sqrt(etx*etx+ety*ety+etz*etz));

				// /* lefr parameter */
				// _rho = ML1[i][j][k];
				// _u = ML2[i][j][k]/_rho;
				// _v = ML3[i][j][k]/_rho;
				// _w = ML4[i][j][k]/_rho;

				// _U = xix*_u+xiy*_v+xiz*_w;
				// _V = etx*_u+ety*_v+etz*_w;

				// _W = ztx*_u+zty*_v+ztz*_w;


				// __V = ETx*_u+ETy*_v+ETz*_w;


				// _VV = _u*_u+_v*_v+_w*_w;
				// _P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				// _T = _P/_rho;
				// _C = sqrt(K*_P/_rho);
				// _H = 0.5*_VV+_C*_C/(K-1);

				// /* right parameter */
				// rho_ = MR1[i][j][k];
				// u_ = MR2[i][j][k]/rho_;
				// v_ = MR3[i][j][k]/rho_;
				// w_ = MR4[i][j][k]/rho_;

				// U_ = xix*u_+xiy*v_+xiz*w_;
				// V_ = etx*u_+ety*v_+etz*w_;
				// W_ = ztx*u_+zty*v_+ztz*w_;

				// V__ = ETx*u_+ETy*v_+ETz*w_;

				// VV_ = u_*u_+v_*v_+w_*w_;
				// P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				// T_ = P_/rho_;
				// C_ = sqrt(K*P_/rho_);
				// H_ = 0.5*VV_+C_*C_/(K-1);

				// /* flux varibale */
				// rho = 0.5*(_rho+rho_);
				// u = 0.5*(_u+u_);
				// v = 0.5*(_v+v_);
				// w = 0.5*(_w+w_);

				// U = 0.5*(_U+U_);
				// V = 0.5*(_V+V_);
				// W = 0.5*(_W+W_); 

				// _V_ = 0.5*(__V+V__);

				// VV = u*u+v*v+w*w;
				// H = 0.5*(_H+H_);
				// C = sqrt((H-0.5*VV)*(K-1));
				// P = rho*C*C/K;
				// T = P/rho;


				// /* jump dU */
				// dU1 = P_-_P;
				// dU2 = u_-_u;
				// dU3 = v_-_v;
				// dU4 = w_-_w;
				// dU5 = T_-_T;

				// /* preconditioning */
				// if (VV/C/C <= e)
					// beta = e;
				// else
					// beta = VV/C/C;

				// /*V=etx*u+ety*v+etz*w;*/
				// //S=sqrt(V*V*(beta-1)+4*beta*C*C*(etx*etx+ety*ety+etz*etz));
				// S=sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C*(etx*etx+ety*ety+etz*etz));
				// /*_V_=ETx*u+ETy*v+ETz*w;*/
				// _S_=sqrt(_V_*_V_*(beta-1)*(beta-1)+4*beta*C*C*(ETx*ETx+ETy*ETy+ETz*ETz));

				// temp = 4*K*_S_*T*beta;
				// temp1 = S+V+V*beta;
				// temp2 = -S+V+V*beta;
				// temp3 = _S_-_V_+_V_*beta;
				// temp4 = _S_+_V_-_V_*beta;
				// temp5 = _S_*_S_-(beta-1)*(beta-1)*_V_*_V_;
				// tempu1 = 2*T*ETx*K*beta+u*_S_+u*(beta-1)*_V_;
				// tempu2 = 2*T*ETx*K*beta-u*_S_+u*(beta-1)*_V_;
				// tempv1 = 2*T*ETy*K*beta+v*_S_+v*(beta-1)*_V_;
				// tempv2 = 2*T*ETy*K*beta-v*_S_+v*(beta-1)*_V_;
				// tempw1 = 2*T*ETz*K*beta+w*_S_+w*(beta-1)*_V_;
				// tempw2 = 2*T*ETz*K*beta-w*_S_+w*(beta-1)*_V_;
				// tempuvw = 2*T*K*beta*(u*ETx+v*ETy+w*ETz);
				// temph1 = tempuvw+H*_S_+H*_V_*beta-H*_V_; 
				// temph2 = tempuvw-H*_S_+H*_V_*beta-H*_V_;


				// d21 = 1/temp*(u*_S_*(4*(K-1)*beta*fabs(V)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ETx*K*beta+u*_V_*(beta-1)));
				// d22 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ETy*ETy+ETz*ETz)*fabs(V)+ETx*(temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1))));
				// d23 = 1/(2*temp)*(rho*ETy*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				// d24 = 1/(2*temp)*(rho*ETz*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				// d25 = -u*rho*fabs(V)/T;


				// Fav2 = d21*dU1+d22*dU2+d23*dU3+d24*dU4+d25*dU5;


				// MUL[i][j][k] = _u;
				// MUR[i][j][k] = u_;

				// MVL[i][j][k] = _v;
				// MVR[i][j][k] = v_;

				// MrhoL[i][j][k] = _rho;
				// MrhoR[i][j][k] = rho_;

				// MRUV[i][j][k] = 0.5*(_rho*_u*_v+rho_*u_*v_);

				// MRoe_D[i][j][k] = 0.5*EpY[i][j][k]*Fav2;
				// MdRL[i][j][k] = dU2;
				
			// }
		// }
	// }




	// double 	irho,iU,iV,iW,iP,iT,
			// jrho,jU,jV,jW,jP,jT,
			// krho,kU,kV,kW,kP,kT,

			// rhoi,Ui,Vi,Wi,Pi,Ti,
			// rhoj,Uj,Vj,Wj,Pj,Tj,
			// rhok,Uk,Vk,Wk,Pk,Tk,

			// ijrho,ijU,ijV,ijW,ijP,ijT,
			// jkrho,jkU,jkV,jkW,jkP,jkT,
			// ikrho,ikU,ikV,ikW,ikP,ikT,

			// rhoij,Uij,Vij,Wij,Pij,Tij,
			// rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,
			// rhoik,Uik,Vik,Wik,Pik,Tik,

			// irhoj,iUj,iVj,iWj,iPj,iTj,
			// irhok,iUk,iVk,iWk,iPk,iTk,

			// jrhoi,jUi,jVi,jWi,jPi,jTi,
			// jrhok,jUk,jVk,jWk,jPk,jTk,

			// krhoj,kUj,kVj,kWj,kPj,kTj,
			// krhoi,kUi,kVi,kWi,kPi,kTi,

			// du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
			// du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,
			// duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,
			// dT_dx,dT_dy,dT_dz,

			// mu_E,mu_T, Pr_E;

	// double a1,a2,a3,a4,a5,a6,a7,
		   // b1,b2,b3,b4,b5,b6,b7,
		   // c1,c2,c3,c4,c5,c6,c7,
		   // d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,
		   // e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,
		   // f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;

	// double Ux,Uy,Uz,
		   // Vx,Vy,Vz,
		   // Wx,Wy,Wz;




	// double invXI = 1./(deltaXI);
	// double invET = 1./(deltaET);
	// double invZT = 1./(deltaZT);

	// double inv4XI = 1./(4*deltaXI);
	// double inv4ET = 1./(4*deltaET);
	// double inv4ZT = 1./(4*deltaZT);


// //// ============================ ////
	// istart = 3;				      ////	
// //// ============================ ////
	// iend = gend[myid]+1;	      ////
// //// ============================ ////

// #pragma omp parallel for private(\
	// xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
	// a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6,c7,\
	// d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,\
	// f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,\
	// rho,U,V,W,VV,P,T,\
	// irho,iU,iV,iW,iP,iT,\
	// jrho,jU,jV,jW,jP,jT,\
	// krho,kU,kV,kW,kP,kT,\
	// rhoi,Ui,Vi,Wi,Pi,Ti,\
	// rhoj,Uj,Vj,Wj,Pj,Tj,\
	// rhok,Uk,Vk,Wk,Pk,Tk,\
	// ijrho,ijU,ijV,ijW,ijP,ijT,\
	// jkrho,jkU,jkV,jkW,jkP,jkT,\
	// ikrho,ikU,ikV,ikW,ikP,ikT,\
	// rhoij,Uij,Vij,Wij,Pij,Tij,\
	// rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,\
	// rhoik,Uik,Vik,Wik,Pik,Tik,\
	// irhoj,iUj,iVj,iWj,iPj,iTj,\
	// irhok,iUk,iVk,iWk,iPk,iTk,\
	// jrhoi,jUi,jVi,jWi,jPi,jTi,\
	// jrhok,jUk,jVk,jWk,jPk,jTk,\
	// krhoj,kUj,kVj,kWj,kPj,kTj,\
	// krhoi,kUi,kVi,kWi,kPi,kTi,\
	// Ux,Uy,Uz,\
	// Vx,Vy,Vz,\
	// Wx,Wy,Wz,\
	// du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,\
	// du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,\
	// duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,\
	// dT_dx,dT_dy,dT_dz,mu_E,Pr_E,\
	// _j,_k,__j,__k,j,k\
	// )

	// for (i = istart; i <= iend; i++) {
		// for (__j = 0,_j = 1, j = 2; j <= nyy; __j++,_j++, j++) {
			// for (__k = 0,_k = 1, k = 2; k <= nzz; __k++,_k++, k++) {


				// // ==== 0 0 0 ==== //
				// rho = U1_[i][j][k];
				// U = U2_[i][j][k]/rho;
				// V = U3_[i][j][k]/rho;
				// W = U4_[i][j][k]/rho;     
				// VV = U*U+V*V+W*W;
				// P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
				// T = P/(rho*R)/J[i][j][k];


				// // ==== -1 1 0 ==== //
				// irhoj = U1_[i-1][j+1][k];
				// iUj = U2_[i-1][j+1][k]/irhoj;
				// iVj = U3_[i-1][j+1][k]/irhoj;
				// iWj = U4_[i-1][j+1][k]/irhoj;
				// iPj = (U5_[i-1][j+1][k]-0.5*irhoj*(iUj*iUj+iVj*iVj+iWj*iWj))*(K-1)*J[i-1][j+1][k];
				// iTj = iPj/(irhoj*R)/J[i-1][j+1][k];

				// // ==== -1 0 1 ==== //
				// irhok = U1_[i-1][j][k+1];
				// iUk = U2_[i-1][j][k+1]/irhok;
				// iVk = U3_[i-1][j][k+1]/irhok;
				// iWk = U4_[i-1][j][k+1]/irhok;
				// iPk = (U5_[i-1][j][k+1]-0.5*irhok*(iUk*iUk+iVk*iVk+iWk*iWk))*(K-1)*J[i-1][j][k+1];
				// iTk = iPk/(irhok*R)/J[i-1][j][k+1];

				// // ==== 1 -1 0 ==== //
				// jrhoi = U1_[i+1][j-1][k];
				// jUi = U2_[i+1][j-1][k]/jrhoi;
				// jVi = U3_[i+1][j-1][k]/jrhoi;
				// jWi = U4_[i+1][j-1][k]/jrhoi;
				// jPi = (U5_[i+1][j-1][k]-0.5*jrhoi*(jUi*jUi+jVi*jVi+jWi*jWi))*(K-1)*J[i+1][j-1][k];
				// jTi = jPi/(jrhoi*R)/J[i+1][j-1][k];

				// // ==== 1 0 -1 ==== //
				// krhoi = U1_[i+1][j][k-1];
				// kUi = U2_[i+1][j][k-1]/krhoi;
				// kVi = U3_[i+1][j][k-1]/krhoi;
				// kWi = U4_[i+1][j][k-1]/krhoi;
				// kPi = (U5_[i+1][j][k-1]-0.5*krhoi*(kUi*kUi+kVi*kVi+kWi*kWi))*(K-1)*J[i+1][j][k-1];
				// kTi = kPi/(krhoi*R)/J[i+1][j][k-1];


				// // ==== 0 -1 1 ==== //
				// jrhok = U1_[i][j-1][k+1];
				// jUk = U2_[i][j-1][k+1]/jrhok;
				// jVk = U3_[i][j-1][k+1]/jrhok;
				// jWk = U4_[i][j-1][k+1]/jrhok;
				// jPk = (U5_[i][j-1][k+1]-0.5*jrhok*(jUk*jUk+jVk*jVk+jWk*jWk))*(K-1)*J[i][j-1][k+1];
				// jTk = jPk/(jrhok*R)/J[i][j-1][k+1];


				// // ==== 0 1 -1 ==== //
				// krhoj = U1_[i][j+1][k-1];
				// kUj = U2_[i][j+1][k-1]/krhoj;
				// kVj = U3_[i][j+1][k-1]/krhoj;
				// kWj = U4_[i][j+1][k-1]/krhoj;
				// kPj = (U5_[i][j+1][k-1]-0.5*krhoj*(kUj*kUj+kVj*kVj+kWj*kWj))*(K-1)*J[i][j+1][k-1];
				// kTj = kPj/(krhoj*R)/J[i][j+1][k-1];

				// // ==== -1 0 0 ==== //
				// irho = U1_[i-1][j][k];
				// iU = U2_[i-1][j][k]/irho;
				// iV = U3_[i-1][j][k]/irho;
				// iW = U4_[i-1][j][k]/irho;
				// iP = (U5_[i-1][j][k]-0.5*irho*(iU*iU+iV*iV+iW*iW))*(K-1)*J[i-1][j][k];
				// iT = iP/(irho*R)/J[i-1][j][k];

				// // ==== 0 -1 0 ==== //
				// jrho = U1_[i][j-1][k];
				// jU = U2_[i][j-1][k]/jrho;
				// jV = U3_[i][j-1][k]/jrho;
				// jW = U4_[i][j-1][k]/jrho;
				// jP = (U5_[i][j-1][k]-0.5*jrho*(jU*jU+jV*jV+jW*jW))*(K-1)*J[i][j-1][k];
				// jT = jP/(jrho*R)/J[i][j-1][k];

				// // ==== 0 0 -1 ==== //
				// krho = U1_[i][j][k-1];
				// kU = U2_[i][j][k-1]/krho;
				// kV = U3_[i][j][k-1]/krho;
				// kW = U4_[i][j][k-1]/krho;
				// kP = (U5_[i][j][k-1]-0.5*krho*(kU*kU+kV*kV+kW*kW))*(K-1)*J[i][j][k-1];
				// kT = kP/(krho*R)/J[i][j][k-1];

				// // ==== -1 -1 0 ==== //
				// ijrho = U1_[i-1][j-1][k];
				// ijU = U2_[i-1][j-1][k]/ijrho;
				// ijV = U3_[i-1][j-1][k]/ijrho;
				// ijW = U4_[i-1][j-1][k]/ijrho;
				// ijP = (U5_[i-1][j-1][k]-0.5*ijrho*(ijU*ijU+ijV*ijV+ijW*ijW))*(K-1)*J[i-1][j-1][k];
				// ijT = ijP/(ijrho*R)/J[i-1][j-1][k];


				// // ==== 0 -1 -1 ==== //
				// jkrho = U1_[i][j-1][k-1];
				// jkU = U2_[i][j-1][k-1]/jkrho;
				// jkV = U3_[i][j-1][k-1]/jkrho;
				// jkW = U4_[i][j-1][k-1]/jkrho;
				// jkP = (U5_[i][j-1][k-1]-0.5*jkrho*(jkU*jkU+jkV*jkV+jkW*jkW))*(K-1)*J[i][j-1][k-1];
				// jkT = jkP/(jkrho*R)/J[i][j-1][k-1];


				// // ==== -1 0 -1 ==== //
				// ikrho = U1_[i-1][j][k-1];
				// ikU = U2_[i-1][j][k-1]/ikrho;
				// ikV = U3_[i-1][j][k-1]/ikrho;
				// ikW = U4_[i-1][j][k-1]/ikrho;
				// ikP = (U5_[i-1][j][k-1]-0.5*ikrho*(ikU*ikU+ikV*ikV+ikW*ikW))*(K-1)*J[i-1][j][k-1];
				// ikT = ikP/(ikrho*R)/J[i-1][j][k-1];
				



				// // ==== 1 0 0 ==== //
				// rhoi = U1_[i+1][j][k];
				// Ui = U2_[i+1][j][k]/rhoi;
				// Vi = U3_[i+1][j][k]/rhoi;
				// Wi = U4_[i+1][j][k]/rhoi;
				// Pi = (U5_[i+1][j][k]-0.5*rhoi*(Ui*Ui+Vi*Vi+Wi*Wi))*(K-1)*J[i+1][j][k];
				// Ti = Pi/(rhoi*R)/J[i+1][j][k];

				// // ==== 0 1 0 ==== //
				// rhoj = U1_[i][j+1][k];
				// Uj = U2_[i][j+1][k]/rhoj;
				// Vj = U3_[i][j+1][k]/rhoj;
				// Wj = U4_[i][j+1][k]/rhoj;
				// Pj = (U5_[i][j+1][k]-0.5*rhoj*(Uj*Uj+Vj*Vj+Wj*Wj))*(K-1)*J[i][j+1][k];
				// Tj = Pj/(rhoj*R)/J[i][j+1][k];

				// // ==== 0 0 1  ==== //
				// rhok = U1_[i][j][k+1];
				// Uk = U2_[i][j][k+1]/rhok;
				// Vk = U3_[i][j][k+1]/rhok;
				// Wk = U4_[i][j][k+1]/rhok;
				// Pk = (U5_[i][j][k+1]-0.5*rhok*(Uk*Uk+Vk*Vk+Wk*Wk))*(K-1)*J[i][j][k+1];
				// Tk = Pk/(rhok*R)/J[i][j][k+1];
				


				// xix=xidx_v[i-1][j-1][k-1];
				// xiy=xidy_v[i-1][j-1][k-1];
				// xiz=xidz_v[i-1][j-1][k-1];
				// etx=etdx_v[i-1][j-1][k-1];
				// ety=etdy_v[i-1][j-1][k-1];
				// etz=etdz_v[i-1][j-1][k-1];          
				// ztx=ztdx_v[i-1][j-1][k-1];
				// zty=ztdy_v[i-1][j-1][k-1];
				// ztz=ztdz_v[i-1][j-1][k-1];


				// /* derivatives of velocity */
				
				// /* X-direction */
				// du_dx = (Ui+jUi-iU-ijU)*inv4XI;
				// dv_dx = (Vi+jVi-iV-ijV)*inv4XI;
				// dw_dx = (Wi+jWi-iW-ijW)*inv4XI;

				// du2_dx = (Ui*Ui+jUi*jUi-iU*iU-ijU*ijU)*inv4XI;
				// dv2_dx = (Vi*Vi+jVi*jVi-iV*iV-ijV*ijV)*inv4XI;
				// dw2_dx = (Wi*Wi+jWi*jWi-iW*iW-ijW*ijW)*inv4XI;


				// duv_dx = (Ui*Vi+jUi*jVi-iU*iV-ijU*ijV)*inv4XI;
				// dvw_dx = (Vi*Wi+jVi*jWi-iV*iW-ijV*ijW)*inv4XI;
				// duw_dx = (Ui*Wi+jUi*jWi-iU*iW-ijU*ijW)*inv4XI;


				// dT_dx = (Ti+jTi-iT-ijT)*inv4XI;


				// /* Y-direction */
				// du_dy = (U-jU)*invET;
				// dv_dy = (V-jV)*invET;
				// dw_dy = (W-jW)*invET;


				// du2_dy = (U*U-jU*jU)*invET;
				// dv2_dy = (V*V-jV*jV)*invET;
				// dw2_dy = (W*W-jW*jW)*invET;


				// duv_dy = (U*V-jU*jV)*invET;
				// dvw_dy = (V*W-jV*jW)*invET;
				// duw_dy = (U*W-jU*jW)*invET;


				// dT_dy = (T-jT)*invET;


				// /* Z-direction */
				// du_dz = (Uk+jUk-kU-jkU)*inv4ZT;
				// dv_dz = (Vk+jVk-kV-jkV)*inv4ZT;
				// dw_dz = (Wk+jWk-kW-jkW)*inv4ZT;

				// du2_dz = (Uk*Uk+jUk*jUk-kU*kU-jkU*jkU)*inv4ZT;
				// dv2_dz = (Vk*Vk+jVk*jVk-kV*kV-jkV*jkV)*inv4ZT;
				// dw2_dz = (Wk*Wk+jWk*jWk-kW*kW-jkW*jkW)*inv4ZT;


				// duv_dz = (Uk*Vk+jUk*jVk-kU*kV-jkU*jkV)*inv4ZT;
				// dvw_dz = (Vk*Wk+jVk*jWk-kV*kW-jkV*jkW)*inv4ZT;
				// duw_dz = (Uk*Wk+jUk*jWk-kU*kW-jkU*jkW)*inv4ZT;


				// dT_dz = (Tk+jTk-kT-jkT)*inv4ZT;

				// /* viscous*/
				// mu_E = mu_L;
				// Pr_E = Pr_L;


				// b1 = (4./3)*etx*etx+ety*ety+etz*etz;
				// b5 = (1./3)*etx*ety;
				// b7 = (1./3)*etx*etz;

				// d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
				// d5 = xix*ety-(2./3)*xiy*etx;
				// d6 = xix*etz-(2./3)*xiz*etx;
				
				// f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
				// f7 = ety*ztx-(2./3)*etx*zty;
				// f9 = etz*ztx-(2./3)*etx*ztz;
				

				// /* Y viscid fluxes */
				// ML2[i][j][k] = mu_E*(d1*du_dx+d5*dv_dx+d6*dw_dx+
					// b1*du_dy+b5*dv_dy+b7*dw_dy+
					// f1*du_dz+f7*dv_dz+f9*dw_dz)/J_v[i-1][j-1][k-1];


				// /**** ES ****/
				// MSXY[i][j-1][k] = ML2[i][j][k];
				// /**** ES ****/

			// }
		// }
	// }





// //// ============================================ ////
		// istart = 2;	         					  ////	
// //// ============================================ ////
		// iend = gend[myid]-1;	    			  ////
// //// ============================================ ////


		// for (i = istart; i <= iend; i++) {
			// for (j = 1; j < nyy; j++) {
				// for (k = 1; k < nz; k++) {

					// ULm[i][j] = ULm[i][j]+MUL[i][j][k];
					// URm[i][j] = URm[i][j]+MUR[i][j][k];

					// VLm[i][j] = VLm[i][j]+MVL[i][j][k];
					// VRm[i][j] = VRm[i][j]+MVR[i][j][k];

					// rhoLm[i][j] = rhoLm[i][j]+MrhoL[i][j][k];
					// rhoRm[i][j] = rhoRm[i][j]+MrhoR[i][j][k];

					// RUVm[i][j] = RUVm[i][j]+MRUV[i][j][k];

					// Roe_Dm[i][j] = Roe_Dm[i][j]+MRoe_D[i][j][k];
					// dRLm[i][j] = dRLm[i][j]+MdRL[i][j][k];

				// }
			// }
		// }


// //// ============================================ ////
		 // istart = 3;		     	              ////	
// //// ============================================ ////
		// iend = gend[myid];				    	  ////
// //// ============================================ ////

		// for (i = istart; i <= iend; i++) {
			// for (j = 1; j < nyy; j++) {
				// for (k = 2; k <= nz; k++) {

					// SXYm[i][j] = SXYm[i][j]+MSXY[i][j][k];

				// }
			// }
		// }


// // =================== //
	// istart = 3;        //
	// iend = gend[myid]; //
// // =================== //

	// for (i = istart; i <= iend; i++) {

			// ii = i+istart1;

			// if ((ii-2-nx_inlet)*deltaXI <= obs1 && (ii-1-nx_inlet)*deltaXI >= obs1) {

				// FILE *fptr;
				// sprintf(LESdata,"Shear_stress_obs1.dat");
				// fptr = fopen(LESdata,"a");

				// fprintf(fptr,"%f\t%f\n",SXYm[i][1]/(nz-1),SXYm[i][ny]/(nz-1));

				// fclose(fptr);

			// }


			// if ((ii-2-nx_inlet)*deltaXI <= obs2 && (ii-1-nx_inlet)*deltaXI >= obs2) {

				// FILE *fptr;
				// sprintf(LESdata,"Shear_stress_obs2.dat");
				// fptr = fopen(LESdata,"a");

				// fprintf(fptr,"%f\t%f\n",SXYm[i][1]/(nz-1),SXYm[i][ny]/(nz-1));

				// fclose(fptr);

			// }

			// if ((ii-2-nx_inlet)*deltaXI <= obs3 && (ii-1-nx_inlet)*deltaXI >= obs3) {

				// FILE *fptr;
				// sprintf(LESdata,"Shear_stress_obs3.dat");
				// fptr = fopen(LESdata,"a");

				// fprintf(fptr,"%f\t%f\n",SXYm[i][1]/(nz-1),SXYm[i][ny]/(nz-1));

				// fclose(fptr);

			// }

			// if ((ii-2-nx_inlet)*deltaXI <= obs4 && (ii-1-nx_inlet)*deltaXI >= obs4) {

				// FILE *fptr;
				// sprintf(LESdata,"Shear_stress_obs4.dat");
				// fptr = fopen(LESdata,"a");

				// fprintf(fptr,"%f\t%f\n",SXYm[i][1]/(nz-1),SXYm[i][ny]/(nz-1));

				// fclose(fptr);

			// }

			// if ((ii-2-nx_inlet)*deltaXI <= obs5 && (ii-1-nx_inlet)*deltaXI >= obs5) {

				// FILE *fptr;
				// sprintf(LESdata,"Shear_stress_obs5.dat");
				// fptr = fopen(LESdata,"a");

				// fprintf(fptr,"%f\t%f\n",SXYm[i][1]/(nz-1),SXYm[i][ny]/(nz-1));

				// fclose(fptr);

			// }

	// }





// // ====================================== //
	// if ( (step%statistic_step) == 0 ) {   //   
// // ====================================== //

// double inv = 1./statistic_step/(nz-1);


// // =================== //
	// istart = 3;        //
	// iend = gend[myid]; //
// // =================== //


			// for (i = istart; i <= iend; i++) {
				// for (j = 2; j < nyy; j++) {

					// ULm[i][j] = ULm[i][j]*inv;
					// URm[i][j] = URm[i][j]*inv;

					// VLm[i][j] = VLm[i][j]*inv;
					// VRm[i][j] = VRm[i][j]*inv;

					// rhoLm[i][j] = rhoLm[i][j]*inv;
					// rhoRm[i][j] = rhoRm[i][j]*inv;

					// RUVm[i][j] = RUVm[i][j]*inv;
					// SXYm[i][j] = SXYm[i][j]*inv;

					// Roe_Dm[i][j] = Roe_Dm[i][j]*inv;
					// dRLm[i][j] = dRLm[i][j]*inv;

				// }
			// }

			

// // =================== //
	// istart = 3;        //
	// iend = gend[myid]; //
// // =================== //

	// for (i = istart; i <= iend; i++) {

		// ii = i+istart1;

		// if ((ii-2-nx_inlet)*deltaXI <= obs1 && (ii-1-nx_inlet)*deltaXI >= obs1) {

			// FILE *fptrND;
			// sprintf(LESdata,"ND_obs1_""%0.5d"".dat",step);
			// fptrND = fopen(LESdata,"w");

			// for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm[i][j]);}

			// fclose(fptrND);

			// }     // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //

		// if ((ii-2-nx_inlet)*deltaXI <= obs2 && (ii-1-nx_inlet)*deltaXI >= obs2) {

			// FILE *fptrND;
			// sprintf(LESdata,"ND_obs2_""%0.5d"".dat",step);
			// fptrND = fopen(LESdata,"w");

			// for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm[i][j]);}

			// fclose(fptrND);

			// }     // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //

		// if ((ii-2-nx_inlet)*deltaXI <= obs3 && (ii-1-nx_inlet)*deltaXI >= obs3) {

			// FILE *fptrND;
			// sprintf(LESdata,"ND_obs3_""%0.5d"".dat",step);
			// fptrND = fopen(LESdata,"w");

			// for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm[i][j]);}

			// fclose(fptrND);

			// }     // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //

		// if ((ii-2-nx_inlet)*deltaXI <= obs4 && (ii-1-nx_inlet)*deltaXI >= obs4) {

			// FILE *fptrND;
			// sprintf(LESdata,"ND_obs4_""%0.5d"".dat",step);
			// fptrND = fopen(LESdata,"w");

			// for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm[i][j]);}

			// fclose(fptrND);

			// }     // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //

		// if ((ii-2-nx_inlet)*deltaXI <= obs5 && (ii-1-nx_inlet)*deltaXI >= obs5) {

			// FILE *fptrND;
			// sprintf(LESdata,"ND_obs5_""%0.5d"".dat",step);
			// fptrND = fopen(LESdata,"w");

			// for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm[i][j]); }
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm[i][j]);}
			// for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm[i][j]);}

			// fclose(fptrND);

			// }     // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //

// // ========================= //
		// }					 //
// // ========================= //




// // =================== //
	// istart = 3;        //
	// iend = gend[myid]; //
// // =================== //

			// for (i = istart; i <= iend; i++) {
// #pragma omp parallel for private(j)
				// for (j = 1; j <= ny; j++) {

					// ULm[i][j] = 0;
					// URm[i][j] = 0;

					// VLm[i][j] = 0;
					// VRm[i][j] = 0;

					// rhoLm[i][j] = 0;
					// rhoRm[i][j] = 0;

					// RUVm[i][j] = 0;
					// SXYm[i][j] = 0;

					// Roe_Dm[i][j] = 0;
					// dRLm[i][j] = 0;

				// }
			// }

// // ====================================== //
	// }                                     //
// // ====================================== //



	// delete [] MUL;
	// delete [] MUR;

	// delete [] MVL;
	// delete [] MVR;

	// delete [] MrhoL;
	// delete [] MrhoR;

	// delete [] MRUV;
	// delete [] MSXY;

	// delete [] MRoe_D;
	// delete [] MdRL;



/**** mean profile and turbulenct intensties ****/




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {

					rho = U1_[i][j][k];
					U = U2_[i][j][k]/rho;
					V = U3_[i][j][k]/rho;
					W = U4_[i][j][k]/rho;     
					VV = U*U+V*V+W*W;
					P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
					T = P/(rho*R)/J[i][j][k];

					Um[i][j] = Um[i][j]+U;
					Vm[i][j] = Vm[i][j]+V;
					Wm[i][j] = Wm[i][j]+W;
					Tm[i][j] = Tm[i][j]+T;

					
					UUm[i][j] = UUm[i][j]+U*U;
					VVm[i][j] = VVm[i][j]+V*V;
					WWm[i][j] = WWm[i][j]+W*W;
					TTm[i][j] = TTm[i][j]+T*T;

					UVm[i][j] = UVm[i][j]+U*V;

					UTm[i][j] = UTm[i][j]+U*T;
					VTm[i][j] = VTm[i][j]+V*T;

				}
			}

		}



// ======================================================== //
	if ( (step%statistic_step) == 0) {					    //   
// ======================================================== //
		
		double inv = 1./statistic_step/(nz-1);


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {

				Um[i][j] = Um[i][j]*inv;
				Vm[i][j] = Vm[i][j]*inv;
				Wm[i][j] = Wm[i][j]*inv;
				Tm[i][j] = Tm[i][j]*inv;

				UUm[i][j] = UUm[i][j]*inv;
				VVm[i][j] = VVm[i][j]*inv;
				WWm[i][j] = WWm[i][j]*inv;
				TTm[i][j] = TTm[i][j]*inv;

				UVm[i][j] = UVm[i][j]*inv;
				UTm[i][j] = UTm[i][j]*inv;
				VTm[i][j] = VTm[i][j]*inv;

			}
		}


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //


	for (i = istart; i <= iend; i++) {

		ii = i+istart1;

		if ((ii-2-nx_inlet)*deltaXI <= obs1 && (ii-1-nx_inlet)*deltaXI >= obs1) {

			FILE *fptr;
			sprintf(LESdata,"st_obs1_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs2 && (ii-1-nx_inlet)*deltaXI >= obs2) {

			FILE *fptr;
			sprintf(LESdata,"st_obs2_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs1 && (ii-1)*deltaXI >= obs1) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs3 && (ii-1-nx_inlet)*deltaXI >= obs3) {

			FILE *fptr;
			sprintf(LESdata,"st_obs3_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs3 && (ii-1)*deltaXI >= obs3) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs4 && (ii-1-nx_inlet)*deltaXI >= obs4) {

			FILE *fptr;
			sprintf(LESdata,"st_obs4_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }

			fclose(fptr);

		}    // ---- if ((ii-2)*deltaXI <= obs4 && (ii-1)*deltaXI >= obs4) ---- //


		if ((ii-2-nx_inlet)*deltaXI <= obs5 && (ii-1-nx_inlet)*deltaXI >= obs5) {

			FILE *fptr;
			sprintf(LESdata,"st_obs5_""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Tm[i][j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",TTm[i][j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UTm[i][j]); }
			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VTm[i][j]); }

			fclose(fptr);

		}    // ---- if ((i+1)*deltaXI <= obs5 && (i+2)*deltaXI >= obs5) ---- //


// ========================= //
		}					 //
// ========================= //


/**** mean temperature output end ****/



		
	int nx_out = X_out;
	int ny_out = Y_out;

	double (*Tm_out)[Y_out] = new double[X_out+1][Y_out];
	double (*Um_out)[Y_out] = new double[X_out+1][Y_out];

	ii = 2;


//// ========================= ////
	  istart = gstart[myid];   ////	
//// ========================= ////
	  iend = gend0[myid];      ////
//// ========================= ////


	  for (i = istart; i <= iend; i++) {

		  ii = ii+1;

		  for (j = 0; j < ny-1; j++) { 

			  Tm_out[i][j] =Tm[ii][j+2]; 
			  Um_out[i][j] =Um[ii][j+2]; 

		  }
	  }


	  MPI_Comm comm;
	  comm=MPI_COMM_WORLD;
	  MPI_Status istat[8];

	  if (myid > 0) {

		  istart=gstart[myid];
		  icount=gcount[myid]*Y_out;
		  idest=0;

		  itag = 210;
		  MPI_Send((void *)&Tm_out[istart][0], icount, MPI_DOUBLE, idest, itag, comm);
		  itag = 220;
		  MPI_Send((void *)&Um_out[istart][0], icount, MPI_DOUBLE, idest, itag, comm);

	  }

	  else {

		  for ( isrc=1; isrc < nproc; isrc++ ) {

			  istart=gstart[isrc];
			  icount=gcount[isrc]*Y_out;

			  itag = 210;
			  MPI_Recv((void *)&Tm_out[istart][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
			  itag = 220;
			  MPI_Recv((void *)&Um_out[istart][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
			  
		  }

	  }


	  if (myid == 0) {

			char LESdata[100];
			FILE *fptr;
			sprintf(LESdata,"Tave""%0.5d"".bin",step);
			fptr = fopen(LESdata,"wb");

			if (fptr != NULL)
			{

				fwrite(Tm_out,sizeof(double),X_out*Y_out,fptr);
				
				fclose(fptr);

			}

			else printf("File opening Failure\n");
			
			
			sprintf(LESdata,"Uave""%0.5d"".bin",step);
			fptr = fopen(LESdata,"wb");

			if (fptr != NULL)
			{

				fwrite(Um_out,sizeof(double),X_out*Y_out,fptr);
				
				fclose(fptr);

			}

			else printf("File opening Failure\n");

		}    // ---- if (myid == 0) ---- //


	  

	// =================== //
		istart = 3;        //
		iend = gend[myid]; //
	// =================== //

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {

				Um[i][j] = 0;
				Vm[i][j] = 0;
				Wm[i][j] = 0;
				Tm[i][j] = 0;

				UUm[i][j] = 0;
				VVm[i][j] = 0;
				WWm[i][j] = 0;
				TTm[i][j] = 0;

				UVm[i][j] = 0;
				UTm[i][j] = 0;
				VTm[i][j] = 0;

			}
		}


		



	delete [] Tm_out;
	delete [] Um_out;



// ======================================================== //
	}  												        //
// ======================================================== //

/**** mean profile and turbulenct intensties - end ****/







}

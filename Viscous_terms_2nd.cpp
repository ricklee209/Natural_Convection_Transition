




#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>

#include "Resolution.h"

void Viscous_terms
(
// ============================================================================ //
int myid,

double (*U1_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U2_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U3_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U4_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U5_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*vF2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*vF3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*vF4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],
double (*vF5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*J)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],

double (*xidx)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidy)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidz)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*etdx)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdy)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdz)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ztdx)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdy)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdz)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xdxi)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdet)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdzt)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ydxi)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydet)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydzt)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*zdxi)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdet)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdzt)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],

double (*J_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xidx_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidy_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidz_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*etdx_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdy_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdz_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ztdx_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdy_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdz_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xdxi_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdet_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdzt_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ydxi_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydet_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydzt_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*zdxi_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdet_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdzt_u)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],

double (*J_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xidx_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidy_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidz_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*etdx_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdy_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdz_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ztdx_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdy_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdz_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xdxi_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdet_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdzt_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ydxi_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydet_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydzt_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*zdxi_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdet_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdzt_v)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],

double (*J_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xidx_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidy_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xidz_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*etdx_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdy_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*etdz_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ztdx_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdy_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ztdz_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*xdxi_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdet_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*xdzt_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*ydxi_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydet_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ydzt_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*zdxi_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdet_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*zdzt_w)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory],


/**** LR = vFx2_2 ****/
double (*LR1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LR2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LR3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LR4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LR5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
/**** LR = vFx2_2-end ****/


/**** LL = vFx2_1 ****/
double (*LL1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LL2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LL3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LL4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*LL5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
/**** LL = vFx2_1-end ****/


/**** MR = vFy2_2 ****/
double (*MR1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*MR2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*MR3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*MR4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*MR5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
/**** MR = vFy2_2 ****/


/**** ML = vFy2_1 ****/
double (*ML1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ML2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ML3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ML4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*ML5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
/**** ML = vFy2_1-end ****/


/**** NR = vFz2_2 ****/
double (*NR1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NR2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NR3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NR4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NR5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
/**** NR = vFz2_2-end ****/


/**** NL = vFz2_1 ****/
double (*NL1)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NL2)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NL3)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NL4)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*NL5)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory] 
/**** NL = vFz2_1-end ****/
// ============================================================================ //
)

{

	#include "ijk.h"
	#include "Viscous_terms.h"

	#include "MPI_prm.h"
	#include "Mpi_division.h"

	double rho,U,V,W,VV,P,C,T,h,H;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;

	double 	irho,iirho,iU,iiU,iV,iiV,iW,iiW,iP,iiP,iT,iiT,
			jrho,jjrho,jU,jjU,jV,jjV,jW,jjW,jP,jjP,jT,jjT,
			krho,kkrho,kU,kkU,kV,kkV,kW,kkW,kP,kkP,kT,kkT,

			du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
			du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,
			duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,
			dT_dx,dT_dy,dT_dz,

			mu_E,mu_T, Pr_E;

	double a1,a2,a3,a4,a5,a6,a7,
		   b1,b2,b3,b4,b5,b6,b7,
		   c1,c2,c3,c4,c5,c6,c7,
		   d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,
		   e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,
		   f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;
	

	double invXI = 1./(deltaXI);
	double invET = 1./(deltaET);
	double invZT = 1./(deltaZT);
	
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = gstart[myid];		     	  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]-1;  ////
		else iend = gend[myid]+1;				  ////
//// ============================================ ////
		
#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
	a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6,c7,\
	d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,\
	f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,\
	rho,U,V,W,VV,P,T,\
	irho,iU,iV,iW,iP,iT,\
	jrho,jU,jV,jW,jP,jT,\
	krho,kU,kV,kW,kP,kT,\
	du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,\
	du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,\
	duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,\
	dT_dx,dT_dy,dT_dz,mu_E,Pr_E,\
	_j,_k,__j,__k,j,k\
	)
	
		for (i = istart; i <= iend; i++) {
			for (__j = 0,_j = 1, j = 2; j <= nyy; __j++,_j++, j++) {
				for (__k = 0,_k = 1, k = 2; k <= nzz; __k++,_k++, k++) {


					rho = U1_[i][j][k];
					U = U2_[i][j][k]/rho;
					V = U3_[i][j][k]/rho;
					W = U4_[i][j][k]/rho;     
					VV = U*U+V*V+W*W;
					P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
					T = P/(rho*R)/J[i][j][k];


					/* forward lower */
					jrho = U1_[i][_j][k];
					jU = U2_[i][_j][k]/jrho;
					jV = U3_[i][_j][k]/jrho;
					jW = U4_[i][_j][k]/jrho;
					jP = (U5_[i][_j][k]-0.5*jrho*(jU*jU+jV*jV+jW*jW))*(K-1)*J[i][_j][k];
					jT = jP/(jrho*R)/J[i][_j][k];

					/* forward farther */
					krho = U1_[i][j][_k];
					kU = U2_[i][j][_k]/krho;
					kV = U3_[i][j][_k]/krho;
					kW = U4_[i][j][_k]/krho;
					kP = (U5_[i][j][_k]-0.5*krho*(kU*kU+kV*kV+kW*kW))*(K-1)*J[i][j][_k];
					kT = kP/(krho*R)/J[i][j][_k];

					/* backward cells */
					/* backward */
					irho = U1_[i-1][j][k];
					iU = U2_[i-1][j][k]/irho;
					iV = U3_[i-1][j][k]/irho;
					iW = U4_[i-1][j][k]/irho;
					iP = (U5_[i-1][j][k]-0.5*irho*(iU*iU+iV*iV+iW*iW))*(K-1)*J[i-1][j][k];
					iT = iP/(irho*R)/J[i-1][j][k];


					/* viscous*/
					mu_E = mu_L;
					Pr_E = Pr_L;


					/* derivatives of velocity */
					/* X-direction */
					du_dx = (U-iU)*invXI;
					dv_dx = (V-iV)*invXI;
					dw_dx = (W-iW)*invXI;


					du2_dx = (U*U-iU*iU)*invXI;
					dv2_dx = (V*V-iV*iV)*invXI;
					dw2_dx = (W*W-iW*iW)*invXI;


					duv_dx = (U*V-iU*iV)*invXI;
					dvw_dx = (V*W-iV*iW)*invXI;
					duw_dx = (U*W-iU*iW)*invXI;


					dT_dx = (T-iT)*invXI;


					/* Y-direction */
					du_dy = (U-jU)*invET;
					dv_dy = (V-jV)*invET;
					dw_dy = (W-jW)*invET;

					du2_dy = (U*U-jU*jU)*invET;
					dv2_dy = (V*V-jV*jV)*invET;
					dw2_dy = (W*W-jW*jW)*invET;

					duv_dy = (U*V-jU*jV)*invET;
					dvw_dy = (V*W-jV*jW)*invET;
					duw_dy = (U*W-jU*jW)*invET;


					dT_dy = (T-jT)*invET;


					/* Z-direction */
					du_dz = (U-kU)*invZT;
					dv_dz = (V-kV)*invZT;
					dw_dz = (W-kW)*invZT;

					du2_dz = (U*U-kU*kU)*invZT;
					dv2_dz = (V*V-kV*kV)*invZT;
					dw2_dz = (W*W-kW*kW)*invZT;

					duv_dz = (U*V-kU*kV)*invZT;
					dvw_dz = (V*W-kV*kW)*invZT;
					duw_dz = (U*W-kU*kW)*invZT;


					dT_dz = (T-kT)*invZT;


					xix=xidx_u[i-1][j-1][k-1];
					xiy=xidy_u[i-1][j-1][k-1];
					xiz=xidz_u[i-1][j-1][k-1];
					etx=etdx_u[i-1][j-1][k-1];
					ety=etdy_u[i-1][j-1][k-1];
					etz=etdz_u[i-1][j-1][k-1];          
					ztx=ztdx_u[i-1][j-1][k-1];
					zty=ztdy_u[i-1][j-1][k-1];
					ztz=ztdz_u[i-1][j-1][k-1];
					

					a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
					a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
					a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
					a4 = xix*xix+xiy*xiy+xiz*xiz;
					a5 = (1./3)*xix*xiy;
					a6 = (1./3)*xiy*xiz;
					a7 = (1./3)*xix*xiz;

					b1 = (4./3)*etx*etx+ety*ety+etz*etz;
					b2 = etx*etx+(4./3)*ety*ety+etz*etz;
					b3 = etx*etx+ety*ety+(4./3)*etz*etz;
					b4 = etx*etx+ety*ety+etz*etz;
					b5 = (1./3)*etx*ety;
					b6 = (1./3)*ety*etz;
					b7 = (1./3)*etx*etz;

					c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
					c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
					c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
					c4 = ztx*ztx+zty*zty+ztz*ztz;
					c5 = (1./3)*ztx*zty;
					c6 = (1./3)*zty*ztz;
					c7 = (1./3)*ztx*ztz;

					d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
					d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
					d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
					d4 = xix*etx+xiy*ety+xiz*etz;
					d5 = xix*ety-(2./3)*xiy*etx;
					d6 = xix*etz-(2./3)*xiz*etx;
					d7 = xiy*etx-(2./3)*xix*ety;
					d8 = xiy*etz-(2./3)*xiz*ety;
					d9 = xiz*etx-(2./3)*xix*etz;
					d10 = xiz*ety-(2./3)*xiy*etz;

					e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
					e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
					e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
					e4 = xix*ztx+xiy*zty+xiz*ztz;
					e5 = xix*zty-(2./3)*xiy*ztx;
					e6 = xix*ztz-(2./3)*xiz*ztx;
					e7 = xiy*ztx-(2./3)*xix*zty;
					e8 = xiy*ztz-(2./3)*xiz*zty;
					e9 = xiz*ztx-(2./3)*xix*ztz;
					e10 = xiz*zty-(2./3)*xiy*ztz;

					f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
					f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
					f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
					f4 = etx*ztx+ety*zty+etz*ztz;
					f5 = etx*zty-(2./3)*ety*ztx;
					f6 = etx*ztz-(2./3)*etz*ztx;
					f7 = ety*ztx-(2./3)*etx*zty;
					f8 = ety*ztz-(2./3)*etz*zty;
					f9 = etz*ztx-(2./3)*etx*ztz;
					f10 = etz*zty-(2./3)*ety*ztz;


					/* X viscid fluxes */
					LL2[i][j][k] = mu_E*(a1*du_dx+a5*dv_dx+a7*dw_dx+
						d1*du_dy+d7*dv_dy+d9*dw_dy+
						e1*du_dz+e7*dv_dz+e9*dw_dz)/J_u[i-1][j-1][k-1];

					LL3[i][j][k] = mu_E*(a5*du_dx+a2*dv_dx+a6*dw_dx+
						d5*du_dy+d2*dv_dy+d10*dw_dy+
						e5*du_dz+e2*dv_dz+e10*dw_dz)/J_u[i-1][j-1][k-1];

					LL4[i][j][k] = mu_E*(a7*du_dx+a6*dv_dx+a3*dw_dx+
						d6*du_dy+d8*dv_dy+d3*dw_dy+
						e6*du_dz+e8*dv_dz+e3*dw_dz)/J_u[i-1][j-1][k-1];

					LL5[i][j][k] = mu_E*(0.5*a1*du2_dx+0.5*a2*dv2_dx+0.5*a3*dw2_dx+a5*duv_dx+a6*dvw_dx+a7*duw_dx+Cv*K*a4*dT_dx/(Pr_E)+
						0.5*d1*du2_dy+0.5*d2*dv2_dy+0.5*d3*dw2_dy+d5*V*du_dy+d6*W*du_dy+d7*U*dv_dy+d8*W*dv_dy+d9*U*dw_dy+d10*V*dw_dy+Cv*K*d4*dT_dy/(Pr_E)+
						0.5*e1*du2_dz+0.5*e2*dv2_dz+0.5*e3*dw2_dz+e5*V*du_dz+e6*W*du_dz+e7*U*dv_dz+e8*W*dv_dz+e9*U*dw_dz+e10*V*dw_dz+Cv*K*e4*dT_dz/(Pr_E))/J_u[i-1][j-1][k-1];


					xix=xidx_v[i-1][j-1][k-1];
					xiy=xidy_v[i-1][j-1][k-1];
					xiz=xidz_v[i-1][j-1][k-1];
					etx=etdx_v[i-1][j-1][k-1];
					ety=etdy_v[i-1][j-1][k-1];
					etz=etdz_v[i-1][j-1][k-1];          
					ztx=ztdx_v[i-1][j-1][k-1];
					zty=ztdy_v[i-1][j-1][k-1];
					ztz=ztdz_v[i-1][j-1][k-1];
					

					a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
					a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
					a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
					a4 = xix*xix+xiy*xiy+xiz*xiz;
					a5 = (1./3)*xix*xiy;
					a6 = (1./3)*xiy*xiz;
					a7 = (1./3)*xix*xiz;

					b1 = (4./3)*etx*etx+ety*ety+etz*etz;
					b2 = etx*etx+(4./3)*ety*ety+etz*etz;
					b3 = etx*etx+ety*ety+(4./3)*etz*etz;
					b4 = etx*etx+ety*ety+etz*etz;
					b5 = (1./3)*etx*ety;
					b6 = (1./3)*ety*etz;
					b7 = (1./3)*etx*etz;

					c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
					c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
					c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
					c4 = ztx*ztx+zty*zty+ztz*ztz;
					c5 = (1./3)*ztx*zty;
					c6 = (1./3)*zty*ztz;
					c7 = (1./3)*ztx*ztz;

					d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
					d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
					d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
					d4 = xix*etx+xiy*ety+xiz*etz;
					d5 = xix*ety-(2./3)*xiy*etx;
					d6 = xix*etz-(2./3)*xiz*etx;
					d7 = xiy*etx-(2./3)*xix*ety;
					d8 = xiy*etz-(2./3)*xiz*ety;
					d9 = xiz*etx-(2./3)*xix*etz;
					d10 = xiz*ety-(2./3)*xiy*etz;

					e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
					e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
					e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
					e4 = xix*ztx+xiy*zty+xiz*ztz;
					e5 = xix*zty-(2./3)*xiy*ztx;
					e6 = xix*ztz-(2./3)*xiz*ztx;
					e7 = xiy*ztx-(2./3)*xix*zty;
					e8 = xiy*ztz-(2./3)*xiz*zty;
					e9 = xiz*ztx-(2./3)*xix*ztz;
					e10 = xiz*zty-(2./3)*xiy*ztz;

					f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
					f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
					f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
					f4 = etx*ztx+ety*zty+etz*ztz;
					f5 = etx*zty-(2./3)*ety*ztx;
					f6 = etx*ztz-(2./3)*etz*ztx;
					f7 = ety*ztx-(2./3)*etx*zty;
					f8 = ety*ztz-(2./3)*etz*zty;
					f9 = etz*ztx-(2./3)*etx*ztz;
					f10 = etz*zty-(2./3)*ety*ztz;




					/* Y viscid fluxes */
					ML2[i][j][k] = mu_E*(d1*du_dx+d5*dv_dx+d6*dw_dx+
						b1*du_dy+b5*dv_dy+b7*dw_dy+
						f1*du_dz+f7*dv_dz+f9*dw_dz)/J_v[i-1][j-1][k-1];

					ML3[i][j][k] = mu_E*(d7*du_dx+d2*dv_dx+d8*dw_dx+
						b5*du_dy+b2*dv_dy+b6*dw_dy+
						f5*du_dz+f2*dv_dz+f10*dw_dz)/J_v[i-1][j-1][k-1];

					ML4[i][j][k] = mu_E*(d9*du_dx+d10*dv_dx+d3*dw_dx+
						b7*du_dy+b6*dv_dy+b3*dw_dy+
						f6*du_dz+f8*dv_dz+f3*dw_dz)/J_v[i-1][j-1][k-1];

					ML5[i][j][k] = mu_E*(0.5*b1*du2_dy+0.5*b2*dv2_dy+0.5*b3*dw2_dy+b5*duv_dy+b6*dvw_dy+b7*duw_dy+Cv*K*b4*dT_dy/(Pr_E)+
						0.5*d1*du2_dx+0.5*d2*dv2_dx+0.5*d3*dw2_dx+d5*U*dv_dx+d6*U*dw_dx+d7*V*du_dx+d8*V*dw_dx+d9*W*du_dx+d10*W*dv_dx+Cv*K*d4*dT_dx/(Pr_E)+
						0.5*f1*du2_dz+0.5*f2*dv2_dz+0.5*f3*dw2_dz+f5*V*du_dz+f6*W*du_dz+f7*U*dv_dz+f8*W*dv_dz+f9*U*dw_dz+f10*V*dw_dz+Cv*K*f4*dT_dz/(Pr_E))/J_v[i-1][j-1][k-1];

					xix=xidx_w[i-1][j-1][k-1];
					xiy=xidy_w[i-1][j-1][k-1];
					xiz=xidz_w[i-1][j-1][k-1];
					etx=etdx_w[i-1][j-1][k-1];
					ety=etdy_w[i-1][j-1][k-1];
					etz=etdz_w[i-1][j-1][k-1];          
					ztx=ztdx_w[i-1][j-1][k-1];
					zty=ztdy_w[i-1][j-1][k-1];
					ztz=ztdz_w[i-1][j-1][k-1];
					

					a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
					a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
					a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
					a4 = xix*xix+xiy*xiy+xiz*xiz;
					a5 = (1./3)*xix*xiy;
					a6 = (1./3)*xiy*xiz;
					a7 = (1./3)*xix*xiz;

					b1 = (4./3)*etx*etx+ety*ety+etz*etz;
					b2 = etx*etx+(4./3)*ety*ety+etz*etz;
					b3 = etx*etx+ety*ety+(4./3)*etz*etz;
					b4 = etx*etx+ety*ety+etz*etz;
					b5 = (1./3)*etx*ety;
					b6 = (1./3)*ety*etz;
					b7 = (1./3)*etx*etz;

					c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
					c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
					c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
					c4 = ztx*ztx+zty*zty+ztz*ztz;
					c5 = (1./3)*ztx*zty;
					c6 = (1./3)*zty*ztz;
					c7 = (1./3)*ztx*ztz;

					d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
					d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
					d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
					d4 = xix*etx+xiy*ety+xiz*etz;
					d5 = xix*ety-(2./3)*xiy*etx;
					d6 = xix*etz-(2./3)*xiz*etx;
					d7 = xiy*etx-(2./3)*xix*ety;
					d8 = xiy*etz-(2./3)*xiz*ety;
					d9 = xiz*etx-(2./3)*xix*etz;
					d10 = xiz*ety-(2./3)*xiy*etz;

					e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
					e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
					e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
					e4 = xix*ztx+xiy*zty+xiz*ztz;
					e5 = xix*zty-(2./3)*xiy*ztx;
					e6 = xix*ztz-(2./3)*xiz*ztx;
					e7 = xiy*ztx-(2./3)*xix*zty;
					e8 = xiy*ztz-(2./3)*xiz*zty;
					e9 = xiz*ztx-(2./3)*xix*ztz;
					e10 = xiz*zty-(2./3)*xiy*ztz;

					f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
					f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
					f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
					f4 = etx*ztx+ety*zty+etz*ztz;
					f5 = etx*zty-(2./3)*ety*ztx;
					f6 = etx*ztz-(2./3)*etz*ztx;
					f7 = ety*ztx-(2./3)*etx*zty;
					f8 = ety*ztz-(2./3)*etz*zty;
					f9 = etz*ztx-(2./3)*etx*ztz;
					f10 = etz*zty-(2./3)*ety*ztz;

					/* Z viscid fluxes */
					NL2[i][j][k] = mu_E*(e1*du_dx+e5*dv_dx+e6*dw_dx+
						f1*du_dy+f5*dv_dy+f6*dw_dy+
						c1*du_dz+c5*dv_dz+c7*dw_dz)/J_w[i-1][j-1][k-1];

					NL3[i][j][k] = mu_E*(e7*du_dx+e2*dv_dx+e8*dw_dx+
						f7*du_dy+f2*dv_dy+f8*dw_dy+
						c5*du_dz+c2*dv_dz+c6*dw_dz)/J_w[i-1][j-1][k-1];

					NL4[i][j][k] = mu_E*(e9*du_dx+e10*dv_dx+e3*dw_dx+
						f9*du_dy+f10*dv_dy+f3*dw_dy+
						c7*du_dz+c6*dv_dz+c3*dw_dz)/J_w[i-1][j-1][k-1];

					NL5[i][j][k] = mu_E*(0.5*c1*du2_dz+0.5*c2*dv2_dz+0.5*c3*dw2_dz+c5*duv_dz+c6*dvw_dz+c7*duw_dz+Cv*K*c4*dT_dz/(Pr_E)+
						0.5*e1*du2_dx+0.5*e2*dv2_dx+0.5*e3*dw2_dx+e5*U*dv_dx+e6*U*dw_dx+e7*V*du_dx+e8*V*dw_dx+e9*W*du_dx+e10*W*dv_dx+Cv*K*e4*dT_dx/(Pr_E)+
						0.5*f1*du2_dy+0.5*f2*dv2_dy+0.5*f3*dw2_dy+f5*U*dv_dy+f6*U*dw_dy+f7*V*du_dy+f8*V*dw_dy+f9*W*du_dy+f10*W*dv_dy+Cv*K*f4*dT_dy/(Pr_E))/J_w[i-1][j-1][k-1];
					
				
					}
				}
			}


//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = gstart[myid];			      ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]-2;  ////
		else iend = gend[myid];		    		  ////
//// ============================================ ////

		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) {
				for (k = 2, k_ = 3; k <= nz; k++, k_++) {


					vF2[i][j][k] = (LL2[i+1][j][k]-LL2[i][j][k])*invXI+
						(ML2[i][j+1][k]-ML2[i][j][k])*invET+
						(NL2[i][j][k_]-NL2[i][j][k])*invZT;

					vF3[i][j][k] = (LL3[i+1][j][k]-LL3[i][j][k])*invXI+
						(ML3[i][j+1][k]-ML3[i][j][k])*invET+
						(NL3[i][j][k_]-NL3[i][j][k])*invZT;

					vF4[i][j][k] = (LL4[i+1][j][k]-LL4[i][j][k])*invXI+
						(ML4[i][j+1][k]-ML4[i][j][k])*invET+
						(NL4[i][j][k_]-NL4[i][j][k])*invZT;

					vF5[i][j][k] = (LL5[i+1][j][k]-LL5[i][j][k])*invXI+
						(ML5[i][j+1][k]-ML5[i][j][k])*invET+
						(NL5[i][j][k_]-NL5[i][j][k])*invZT;

				}
			}
		}
		
		
	

}
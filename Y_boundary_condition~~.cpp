




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

	
	
	istart = 3;
	iend = gend[myid];

//#pragma omp parallel for private(k,rho,U,V,W,VV,P,temp,T)
	for (i = istart; i <= iend; i++) {
	
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

				
				// TTu = TTb = 0;
				
				// for (j = 2; j <= ny-ny_abs; j++) {
				
					// rho = U1_[i][j][k]*J[i][j][k];
					// U = U2_[i][j][k]/U1_[i][j][k];
					// V = U3_[i][j][k]/U1_[i][j][k];
					// W = U4_[i][j][k]/U1_[i][j][k];     
					// VV = U*U+V*V+W*W;
					// P = (U5_[i][j][k]*J[i][j][k]-0.5*rho*VV)*(K-1);
					// temp = P/rho/R;
					
					// mu_E = mu_L*pow((temp/298.0592),1.5)*(298.0592+110.)/(temp+110.);
				
					// TTu = TTu + Cp*mu_E*temp*deltaET/etdy[i][j][k];
					
					// TTb = TTb + Cp*mu_E*deltaET/etdy[i][j][k];
					
				// }
				
				// Tb = TTu/TTb;
				
				// mu_E = mu_L*pow((Tb/298.0592),1.5)*(298.0592+110.)/(Tb+110.);
				
				// lambda_L = mu_E*Cv*K/Pr_L;

				// T = heat_flux/lambda_L/etdy[i][1][k]*deltaET+Tb;
				
				// U = U2_[i][2][k]/U1_[i][2][k];
				// V = U3_[i][2][k]/U1_[i][2][k];
				// W = U4_[i][2][k]/U1_[i][2][k];     
				// VV = U*U+V*V+W*W;
				// P = (U5_[i][2][k]*J[i][2][k]-0.5*U1_[i][2][k]*J[i][2][k]*VV)*(K-1);
				
				// rho = P/R/T;
				
				
				// U1_[i][1][k] = rho/J[i][1][k];
				// U2_[i][1][k] = -rho*U/J[i][1][k];
				// U3_[i][1][k] = -rho*V/J[i][1][k];
				// U4_[i][1][k] = -rho*W/J[i][1][k];
				// U5_[i][1][k] = (P/(K-1)+0.5*rho*VV)/J[i][1][k];

			}
			
			
			for (k = 2; k <= nz; k++) {
			
				jj = ny-ny_abs+1;
				
				rho = U1_[i][jj-1][k]*J[i][jj-1][k];
				U = U2_[i][jj-1][k]/U1_[i][jj-1][k];
				V = U3_[i][jj-1][k]/U1_[i][jj-1][k];
				W = U4_[i][jj-1][k]/U1_[i][jj-1][k];
				VV = U*U+V*V+W*W;
				P = (U5_[i][jj-1][k]*J[i][jj-1][k]-0.5*rho*VV)*(K-1);
			
			
				U1_[i][jj][k] = U1_[i][jj-1][k]*J[i][jj-1][k]/J[i][jj][k];
				U2_[i][jj][k] = 0;
				U3_[i][jj][k] = 0;
				U4_[i][jj][k] = 0;
				U5_[i][jj][k] = (P/(K-1))/J[i][jj][k];

			}
			
			for (k = 2; k <= nz; k++) {
			
				jj = ny-ny_abs+2;
				
				rho = U1_[i][jj+2][k]*J[i][jj+2][k];
				U = U2_[i][jj+2][k]/U1_[i][jj+2][k];
				V = U3_[i][jj+2][k]/U1_[i][jj+2][k];
				W = U4_[i][jj+2][k]/U1_[i][jj+2][k];
				VV = U*U+V*V+W*W;
				P = (U5_[i][jj+2][k]*J[i][jj+2][k]-0.5*rho*VV)*(K-1);
			
			
				U1_[i][jj][k] = U1_[i][jj+1][k]*J[i][jj+1][k]/J[i][jj][k];
				U2_[i][jj][k] = 0;
				U3_[i][jj][k] = 0;
				U4_[i][jj][k] = 0;
				U5_[i][jj][k] = (P/(K-1))/J[i][jj][k];

			}
			
			
			
		}
	}
	
	
	
	int nx_inp = nx/2;
	
	int nx_outp = nx/2;
	
	
	if( (gstart[myid]) < nx_inp ) {
	
	
//// ============================================ ////
		istart =  2;	
//// ============================================ ////
		if ( (gend0[myid]) < nx_inp )
			iend = gend[myid];
		else
			iend = nx_inp-gstart[myid]+2;
//// ============================================ ////

		for (i = istart ; i <= iend; i++) {
		
			for (k = 2; k <= nz; k++) {

				
				// U1_[i][ny][k] = U1_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U2_[i][ny][k] = U2_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U3_[i][ny][k] = U3_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U4_[i][ny][k] = U4_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U5_[i][ny][k] = U5_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				
				rhoold = U1q[i][nyy][k]*J[i][nyy][k];
				Uold = U2q[i][nyy][k]/U1q[i][nyy][k];
				Vold = U3q[i][nyy][k]/U1q[i][nyy][k];
				Wold = U4q[i][nyy][k]/U1q[i][nyy][k];
				VVold = Uold*Uold+Vold*Vold+Wold*Wold;
				Pold = (U5q[i][nyy][k]*J[i][nyy][k]-0.5*rhoold*VVold)*(K-1);
				Cold = sqrt(K*Pold/rhoold);
				Told = Pold/rhoold;
	
	
				_rho = U1_[i][ny][k]*J[i][ny][k];
				_U = U2_[i][ny][k]/U1_[i][ny][k];
				_V = U3_[i][ny][k]/U1_[i][ny][k];
				_W = U4_[i][ny][k]/U1_[i][ny][k];
				_VV = _U*_U+_V*_V+_W*_W;
				_P = (U5_[i][ny][k]*J[i][ny][k]-0.5*_rho*_VV)*(K-1);
				_C = sqrt(K*_P/_rho);
				_T = _P/_rho;

				rho = U1_[i][nyy][k]*J[i][nyy][k];
				U = U2_[i][nyy][k]/U1_[i][nyy][k];
				V = U3_[i][nyy][k]/U1_[i][nyy][k];
				W = U4_[i][nyy][k]/U1_[i][nyy][k];
				VV = U*U+V*V+W*W;
				P = (U5_[i][nyy][k]*J[i][nyy][k]-0.5*rho*VV)*(K-1);
				C = sqrt(K*P/rho);
				T = P/rho;
				H = 0.5*VV+C*C/(K-1);

				beta = max(VV/C/C,e);
				
				// /* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				// temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				// temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				// temp3 = H-H*K-VV+K*(T+VV);

				// d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d12 = (-4*U*temp2)/temp;
				// d13 = (-4*V*temp2)/temp;
				// d14 = (-4*W*temp2)/temp;
				// d15 = (4*temp2)/temp;

				// d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d23 = (-4*U*V*temp2)/temp;
				// d24 = (-4*U*W*temp2)/temp;
				// d25 = (4*U*temp2)/temp;

				// d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d32 = (-4*U*V*temp2)/temp;
				// d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d34 = (-4*V*W*temp2)/temp;
				// d35 = (4*V*temp2)/temp;

				// d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d42 = (-4*U*W*temp2)/temp;
				// d43 = (-4*V*W*temp2)/temp;
				// d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d45 = (4*W*temp2)/temp;

				// d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d52 = (-4*H*U*temp2)/temp;
				// d53 = (-4*H*V*temp2)/temp;
				// d54 = (-4*H*W*temp2)/temp;
				// d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				

				V_p = 0.5*(beta+1)*V;

				C_p = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dy = (P-_P)/deltaET*etdy[i][nyy][k];

				du_dy = (U-_U)/deltaET*etdy[i][nyy][k];

				dv_dy = (V-_V)/deltaET*etdy[i][nyy][k];

				dw_dy = (W-_W)/deltaET*etdy[i][nyy][k];

				dT_dy = (T-_T)/deltaET*etdy[i][nyy][k];

				L4 = (V_p+C_p)*(dp_dy-rho*(V_p-C_p-V)*dv_dy);

				L5 = L4*(V_p+C_p-V)/(V_p-C_p-V);
				
				//L5 = 0.25*(V_p+C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(P-101300);
				
				if (V > 0) {

					L1 = V*(dT_dy+1/rho*dp_dy*(1-K)/K);

					L2 = V*du_dy;

					L3 = -V*dw_dy;
					
				}

				else {

					L1 = L2 = L3 = 0;
					
					// L1 = 0;
					
					// L2 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(U-0.);
					
					// L3 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(W-0.);
					
				}
				
				
				V = Vold-0.5/(rho*C_p)*(L4-L5)*deltaTau;

				P = Pold-0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V))*deltaTau;

				U = Uold-L2*deltaTau;

				W = Wold+L3*deltaTau;

				T = Told-L1*deltaTau+1/rho*(K+1)/K*(P-Pold);

				rho = P/T;
				
				U1_[i][nyy][k] = rho/J[i][nyy][k];
				U2_[i][nyy][k] = rho*U/J[i][nyy][k];
				U3_[i][nyy][k] = rho*V/J[i][nyy][k];
				U4_[i][nyy][k] = rho*W/J[i][nyy][k];
				U5_[i][nyy][k] = (P/(K-1)+0.5*rho*(U*U+V*V+W*W))/J[i][nyy][k];

				U1q[i][nyy][k] = U1_[i][nyy][k];
				U2q[i][nyy][k] = U2_[i][nyy][k];
				U3q[i][nyy][k] = U3_[i][nyy][k];
				U4q[i][nyy][k] = U4_[i][nyy][k];
				U5q[i][nyy][k] = U5_[i][nyy][k];
				
				
				// KAL1 = 0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V));
				
				// KAL2 = L2;
				
				// KAL3 = 0.5/(rho*C_p)*(L4-L5);
				
				// KAL4 = -L3;
				
				// KAL5 = L1+1./rho*(K-1)/K*KAL1;

				// KAL1 = KAL1/J[i][nyy][k]+(3*U1_[i][nyy][k]-4*U1[i][nyy][k]+U1q[i][nyy][k])/(2*deltaT);
				// KAL2 = KAL2/J[i][nyy][k]+(3*U2_[i][nyy][k]-4*U2[i][nyy][k]+U2q[i][nyy][k])/(2*deltaT);
				// KAL3 = KAL3/J[i][nyy][k]+(3*U3_[i][nyy][k]-4*U3[i][nyy][k]+U3q[i][nyy][k])/(2*deltaT);
				// KAL4 = KAL4/J[i][nyy][k]+(3*U4_[i][nyy][k]-4*U4[i][nyy][k]+U4q[i][nyy][k])/(2*deltaT);
				// KAL5 = KAL5/J[i][nyy][k]+(3*U5_[i][nyy][k]-4*U5[i][nyy][k]+U5q[i][nyy][k])/(2*deltaT);
				
				// U1_[i][nyy][k] = U1_[i][nyy][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				// U2_[i][nyy][k] = U2_[i][nyy][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				// U3_[i][nyy][k] = U3_[i][nyy][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				// U4_[i][nyy][k] = U4_[i][nyy][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				// U5_[i][nyy][k] = U5_[i][nyy][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
				
				
			}
		
		}
		
	} // ---- if( (gstart[myid]) < nx_inlet ) ---- //
	
	
	
	
	
	
	
	if (gend0[myid] > (X_out-nx_outp-1)) {
	
	
//// ============================================ ////

		if ( (gstart[myid]) > (X_out-nx_outp-1))
			istart =  3;	
		else 
			istart =  (X_out-nx_outp)-gstart[myid]+3;

//// ============================================ ////
		iend = gend[myid];			    		  ////
//// ============================================ ////


		for (i = istart ; i <= iend; i++) {
				
			for (k = 2; k <= nz; k++) {
			
				rhoold = U1q[i][nyy][k]*J[i][nyy][k];
				Uold = U2q[i][nyy][k]/U1q[i][nyy][k];
				Vold = U3q[i][nyy][k]/U1q[i][nyy][k];
				Wold = U4q[i][nyy][k]/U1q[i][nyy][k];
				VVold = Uold*Uold+Vold*Vold+Wold*Wold;
				Pold = (U5q[i][nyy][k]*J[i][nyy][k]-0.5*rhoold*VVold)*(K-1);
				Cold = sqrt(K*Pold/rhoold);
				Told = Pold/rhoold;
	
				_rho = U1_[i][ny][k]*J[i][ny][k];
				_U = U2_[i][ny][k]/U1_[i][ny][k];
				_V = U3_[i][ny][k]/U1_[i][ny][k];
				_W = U4_[i][ny][k]/U1_[i][ny][k];
				_VV = _U*_U+_V*_V+_W*_W;
				_P = (U5_[i][ny][k]*J[i][ny][k]-0.5*_rho*_VV)*(K-1);
				_C = sqrt(K*_P/_rho);
				_T = _P/_rho;

				rho = U1_[i][nyy][k]*J[i][nyy][k];
				U = U2_[i][nyy][k]/U1_[i][nyy][k];
				V = U3_[i][nyy][k]/U1_[i][nyy][k];
				W = U4_[i][nyy][k]/U1_[i][nyy][k];
				VV = U*U+V*V+W*W;
				P = (U5_[i][nyy][k]*J[i][nyy][k]-0.5*rho*VV)*(K-1);
				C = sqrt(K*P/rho);
				T = P/rho;
				H = 0.5*VV+C*C/(K-1);

				beta = max(VV/C/C,e);
				
				// /* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				// temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				// temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				// temp3 = H-H*K-VV+K*(T+VV);

				// d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d12 = (-4*U*temp2)/temp;
				// d13 = (-4*V*temp2)/temp;
				// d14 = (-4*W*temp2)/temp;
				// d15 = (4*temp2)/temp;

				// d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d23 = (-4*U*V*temp2)/temp;
				// d24 = (-4*U*W*temp2)/temp;
				// d25 = (4*U*temp2)/temp;

				// d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d32 = (-4*U*V*temp2)/temp;
				// d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d34 = (-4*V*W*temp2)/temp;
				// d35 = (4*V*temp2)/temp;

				// d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d42 = (-4*U*W*temp2)/temp;
				// d43 = (-4*V*W*temp2)/temp;
				// d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d45 = (4*W*temp2)/temp;

				// d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d52 = (-4*H*U*temp2)/temp;
				// d53 = (-4*H*V*temp2)/temp;
				// d54 = (-4*H*W*temp2)/temp;
				// d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				

				V_p = 0.5*(beta+1)*V;

				C_p = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dy = (P-_P)/deltaET*etdy[i][nyy][k];

				du_dy = (U-_U)/deltaET*etdy[i][nyy][k];

				dv_dy = (V-_V)/deltaET*etdy[i][nyy][k];

				dw_dy = (W-_W)/deltaET*etdy[i][nyy][k];

				dT_dy = (T-_T)/deltaET*etdy[i][nyy][k];

				L4 = (V_p+C_p)*(dp_dy-rho*(V_p-C_p-V)*dv_dy);

				L5 = L4*(V_p+C_p-V)/(V_p-C_p-V);
				
				//L5 = 0.25*(V_p+C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(P-101300);
				
				if (V > 0) {

					L1 = V*(dT_dy+1/rho*dp_dy*(1-K)/K);

					L2 = V*du_dy;

					L3 = -V*dw_dy;
					
				}

				else {

					L1 = L2 = L3 = 0;
					
					// L1 = 0;
					
					// L2 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(U-0.);
					
					// L3 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(W-0.);
					
				}
				
				
				V = Vold-0.5/(rho*C_p)*(L4-L5)*deltaTau;

				P = Pold-0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V))*deltaTau;

				U = Uold-L2*deltaTau;

				W = Wold+L3*deltaTau;

				T = Told-L1*deltaTau+1/rho*(K+1)/K*(P-Pold);

				rho = P/T;
				
				U1_[i][nyy][k] = rho/J[i][nyy][k];
				U2_[i][nyy][k] = rho*U/J[i][nyy][k];
				U3_[i][nyy][k] = rho*V/J[i][nyy][k];
				U4_[i][nyy][k] = rho*W/J[i][nyy][k];
				U5_[i][nyy][k] = (P/(K-1)+0.5*rho*(U*U+V*V+W*W))/J[i][nyy][k];

				U1q[i][nyy][k] = U1_[i][nyy][k];
				U2q[i][nyy][k] = U2_[i][nyy][k];
				U3q[i][nyy][k] = U3_[i][nyy][k];
				U4q[i][nyy][k] = U4_[i][nyy][k];
				U5q[i][nyy][k] = U5_[i][nyy][k];
				
				
				// KAL1 = 0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V));
				
				// KAL2 = L2;
				
				// KAL3 = 0.5/(rho*C_p)*(L4-L5);
				
				// KAL4 = -L3;
				
				// KAL5 = L1+1./rho*(K-1)/K*KAL1;

				
				// KAL1 = KAL1/J[i][nyy][k]+(3*U1_[i][nyy][k]-4*U1[i][nyy][k]+U1q[i][nyy][k])/(2*deltaT);
				// KAL2 = KAL2/J[i][nyy][k]+(3*U2_[i][nyy][k]-4*U2[i][nyy][k]+U2q[i][nyy][k])/(2*deltaT);
				// KAL3 = KAL3/J[i][nyy][k]+(3*U3_[i][nyy][k]-4*U3[i][nyy][k]+U3q[i][nyy][k])/(2*deltaT);
				// KAL4 = KAL4/J[i][nyy][k]+(3*U4_[i][nyy][k]-4*U4[i][nyy][k]+U4q[i][nyy][k])/(2*deltaT);
				// KAL5 = KAL5/J[i][nyy][k]+(3*U5_[i][nyy][k]-4*U5[i][nyy][k]+U5q[i][nyy][k])/(2*deltaT);
				
				// U1_[i][nyy][k] = U1_[i][nyy][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				// U2_[i][nyy][k] = U2_[i][nyy][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				// U3_[i][nyy][k] = U3_[i][nyy][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				// U4_[i][nyy][k] = U4_[i][nyy][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				// U5_[i][nyy][k] = U5_[i][nyy][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
				
				
			}
				
				
		}
	}
	
	
	// istart = 2;
	// iend = gend[myid]+1;

	// for (i = istart; i <= iend; i++) {
		// for (k = 2; k <= nz; k++) {
			
			// U1_[i][nyy][k] = U1_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			// U2_[i][nyy][k] = U2_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			// U3_[i][nyy][k] = U3_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			// U4_[i][nyy][k] = U4_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			// U5_[i][nyy][k] = U5_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];	
				
		// }
	// }
	
	
	

	
	
	


	if( (gstart[myid]) < nx_inlet ) {
	
//// ============================================ ////
		istart =  3;	
//// ============================================ ////
		if ( (gend0[myid]) < nx_inlet )
			iend = gend[myid];
		else
			iend = nx_inlet-gstart[myid]+2;
//// ============================================ ////

		for (i = istart ; i <= iend; i++) {
			for (k = 2; k <= nz; k++) {

// ------------------------------------------------------------------------- //			
// ------------------------------- Symmetric ------------------------------- // 
				
				U1_[i][1][k] = U1_[i][2][k]*J[i][2][k]/J[i][1][k];
				U2_[i][1][k] = U2_[i][2][k]*J[i][2][k]/J[i][1][k];
				U3_[i][1][k] = -U3_[i][2][k]*J[i][2][k]/J[i][1][k];
				U4_[i][1][k] = U4_[i][2][k]*J[i][2][k]/J[i][1][k];
				U5_[i][1][k] = U5_[i][2][k]*J[i][2][k]/J[i][1][k];
				
				U1_[i][0][k] = U1_[i][3][k]*J[i][3][k]/J[i][0][k];
				U2_[i][0][k] = U2_[i][3][k]*J[i][3][k]/J[i][0][k];
				U3_[i][0][k] = -U3_[i][3][k]*J[i][3][k]/J[i][0][k];
				U4_[i][0][k] = U4_[i][3][k]*J[i][3][k]/J[i][0][k];
				U5_[i][0][k] = U5_[i][3][k]*J[i][3][k]/J[i][0][k];
				
// ------------------------------- Symmetric ------------------------------- // 				
// ------------------------------------------------------------------------- //	
				
				
				
				// U1_[i][ny][k] = U1_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U2_[i][ny][k] = U2_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U3_[i][ny][k] = U3_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U4_[i][ny][k] = U4_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				// U5_[i][ny][k] = U5_[i][ny-ny_abs][k]*J[i][ny-ny_abs][k]/J[i][ny][k];
				
				
				
				// _rho = U1_[i][ny][k]*J[i][ny][k];
				// _U = U2_[i][ny][k]/U1_[i][ny][k];
				// _V = U3_[i][ny][k]/U1_[i][ny][k];
				// _W = U4_[i][ny][k]/U1_[i][ny][k];
				// _VV = _U*_U+_V*_V+_W*_W;
				// _P = (U5_[i][ny][k]*J[i][ny][k]-0.5*_rho*_VV)*(K-1);
				// _C = sqrt(K*_P/_rho);
				// _T = _P/_rho;

				// rho = U1_[i][nyy][k]*J[i][nyy][k];
				// U = U2_[i][nyy][k]/U1_[i][nyy][k];
				// V = U3_[i][nyy][k]/U1_[i][nyy][k];
				// W = U4_[i][nyy][k]/U1_[i][nyy][k];
				// VV = U*U+V*V+W*W;
				// P = (U5_[i][nyy][k]*J[i][nyy][k]-0.5*rho*VV)*(K-1);
				// C = sqrt(K*P/rho);
				// T = P/rho;
				// H = 0.5*VV+C*C/(K-1);

				// beta = max(VV/C/C,e);
				
				// /* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				// temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				// temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				// temp3 = H-H*K-VV+K*(T+VV);

				// d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d12 = (-4*U*temp2)/temp;
				// d13 = (-4*V*temp2)/temp;
				// d14 = (-4*W*temp2)/temp;
				// d15 = (4*temp2)/temp;

				// d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d23 = (-4*U*V*temp2)/temp;
				// d24 = (-4*U*W*temp2)/temp;
				// d25 = (4*U*temp2)/temp;

				// d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d32 = (-4*U*V*temp2)/temp;
				// d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d34 = (-4*V*W*temp2)/temp;
				// d35 = (4*V*temp2)/temp;

				// d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d42 = (-4*U*W*temp2)/temp;
				// d43 = (-4*V*W*temp2)/temp;
				// d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d45 = (4*W*temp2)/temp;

				// d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d52 = (-4*H*U*temp2)/temp;
				// d53 = (-4*H*V*temp2)/temp;
				// d54 = (-4*H*W*temp2)/temp;
				// d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				

				// V_p = 0.5*(beta+1)*V;

				// C_p = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C);

				// dp_dy = (P-_P)/deltaET*etdy[i][nyy][k];

				// du_dy = (U-_U)/deltaET*etdy[i][nyy][k];

				// dv_dy = (V-_V)/deltaET*etdy[i][nyy][k];

				// dw_dy = (W-_W)/deltaET*etdy[i][nyy][k];

				// dT_dy = (T-_T)/deltaET*etdy[i][nyy][k];

				// L4 = (V_p+C_p)*(dp_dy-rho*(V_p-C_p-V)*dv_dy);

				// L5 = L4*(V_p+C_p-V)/(V_p-C_p-V);
				
				// //L5 = 0.25*(V_p+C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(P-101300);
				
				// if (V > 0) {

					// L1 = V*(dT_dy+1/rho*dp_dy*(1-K)/K);

					// L2 = V*du_dy;

					// L3 = -V*dw_dy;
					
				// }

				// else {

					// L1 = L2 = L3 = 0;
					
					// // L1 = 0;
					
					// // L2 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(U-0.);
					
					// // L3 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(W-0.);
					
				// }
				
				// KAL1 = 0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V));
				
				// KAL2 = L2;
				
				// KAL3 = 0.5/(rho*C_p)*(L4-L5);
				
				// KAL4 = -L3;
				
				// KAL5 = L1+1./rho*(K-1)/K*KAL1;

				// KAL1 = KAL1/J[i][nyy][k]+(3*U1_[i][nyy][k]-4*U1[i][nyy][k]+U1q[i][nyy][k])/(2*deltaT);
				// KAL2 = KAL2/J[i][nyy][k]+(3*U2_[i][nyy][k]-4*U2[i][nyy][k]+U2q[i][nyy][k])/(2*deltaT);
				// KAL3 = KAL3/J[i][nyy][k]+(3*U3_[i][nyy][k]-4*U3[i][nyy][k]+U3q[i][nyy][k])/(2*deltaT);
				// KAL4 = KAL4/J[i][nyy][k]+(3*U4_[i][nyy][k]-4*U4[i][nyy][k]+U4q[i][nyy][k])/(2*deltaT);
				// KAL5 = KAL5/J[i][nyy][k]+(3*U5_[i][nyy][k]-4*U5[i][nyy][k]+U5q[i][nyy][k])/(2*deltaT);
				
				// U1_[i][nyy][k] = U1_[i][nyy][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				// U2_[i][nyy][k] = U2_[i][nyy][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				// U3_[i][nyy][k] = U3_[i][nyy][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				// U4_[i][nyy][k] = U4_[i][nyy][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				// U5_[i][nyy][k] = U5_[i][nyy][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
				
				
			}
		}
		
	}

	
	if (gend0[myid] > (X_out-nx_outlet-1)) {
	
//// ============================================ ////

		if ( (gstart[myid]) > (X_out-nx_outlet-1))
			istart =  3;	
		else 
			istart =  (X_out-nx_outlet)-gstart[myid]+3;

//// ============================================ ////
		iend = gend[myid];			    		  ////
//// ============================================ ////

		for (i = istart ; i <= iend; i++) {
			for (k = 2; k <= nz; k++) {
	
// ------------------------------------------------------------------------- //	
// ------------------------------- Symmetric ------------------------------- //
	
				U1_[i][1][k] = U1_[i][2][k]*J[i][2][k]/J[i][1][k];
				U2_[i][1][k] = U2_[i][2][k]*J[i][2][k]/J[i][1][k];
				U3_[i][1][k] = -U3_[i][2][k]*J[i][2][k]/J[i][1][k];
				U4_[i][1][k] = U4_[i][2][k]*J[i][2][k]/J[i][1][k];
				U5_[i][1][k] = U5_[i][2][k]*J[i][2][k]/J[i][1][k];

				U1_[i][0][k] = U1_[i][3][k]*J[i][3][k]/J[i][0][k];
				U2_[i][0][k] = U2_[i][3][k]*J[i][3][k]/J[i][0][k];
				U3_[i][0][k] = -U3_[i][3][k]*J[i][3][k]/J[i][0][k];
				U4_[i][0][k] = U4_[i][3][k]*J[i][3][k]/J[i][0][k];
				U5_[i][0][k] = U5_[i][3][k]*J[i][3][k]/J[i][0][k];

// ------------------------------- Symmetric ------------------------------- // 				
// ------------------------------------------------------------------------- //	
				
				// _rho = U1_[i][ny][k]*J[i][ny][k];
				// _U = U2_[i][ny][k]/U1_[i][ny][k];
				// _V = U3_[i][ny][k]/U1_[i][ny][k];
				// _W = U4_[i][ny][k]/U1_[i][ny][k];
				// _VV = _U*_U+_V*_V+_W*_W;
				// _P = (U5_[i][ny][k]*J[i][ny][k]-0.5*_rho*_VV)*(K-1);
				// _C = sqrt(K*_P/_rho);
				// _T = _P/_rho;

				// rho = U1_[i][nyy][k]*J[i][nyy][k];
				// U = U2_[i][nyy][k]/U1_[i][nyy][k];
				// V = U3_[i][nyy][k]/U1_[i][nyy][k];
				// W = U4_[i][nyy][k]/U1_[i][nyy][k];
				// VV = U*U+V*V+W*W;
				// P = (U5_[i][nyy][k]*J[i][nyy][k]-0.5*rho*VV)*(K-1);
				// C = sqrt(K*P/rho);
				// T = P/rho;
				// H = 0.5*VV+C*C/(K-1);

				// beta = max(VV/C/C,e);
				
				// /* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				// temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				// temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				// temp3 = H-H*K-VV+K*(T+VV);

				// d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d12 = (-4*U*temp2)/temp;
				// d13 = (-4*V*temp2)/temp;
				// d14 = (-4*W*temp2)/temp;
				// d15 = (4*temp2)/temp;

				// d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d23 = (-4*U*V*temp2)/temp;
				// d24 = (-4*U*W*temp2)/temp;
				// d25 = (4*U*temp2)/temp;

				// d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d32 = (-4*U*V*temp2)/temp;
				// d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d34 = (-4*V*W*temp2)/temp;
				// d35 = (4*V*temp2)/temp;

				// d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d42 = (-4*U*W*temp2)/temp;
				// d43 = (-4*V*W*temp2)/temp;
				// d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				// d45 = (4*W*temp2)/temp;

				// d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				// d52 = (-4*H*U*temp2)/temp;
				// d53 = (-4*H*V*temp2)/temp;
				// d54 = (-4*H*W*temp2)/temp;
				// d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				

				// V_p = 0.5*(beta+1)*V;

				// C_p = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C);

				// dp_dy = (P-_P)/deltaET*etdy[i][nyy][k];

				// du_dy = (U-_U)/deltaET*etdy[i][nyy][k];

				// dv_dy = (V-_V)/deltaET*etdy[i][nyy][k];

				// dw_dy = (W-_W)/deltaET*etdy[i][nyy][k];

				// dT_dy = (T-_T)/deltaET*etdy[i][nyy][k];

				// L4 = (V_p+C_p)*(dp_dy-rho*(V_p-C_p-V)*dv_dy);

				// L5 = L4*(V_p+C_p-V)/(V_p-C_p-V);
				
				// //L5 = 0.25*(V_p+C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(P-101300);
				
				// if (V > 0) {

					// L1 = V*(dT_dy+1/rho*dp_dy*(1-K)/K);

					// L2 = V*du_dy;

					// L3 = -V*dw_dy;
					
				// }

				// else {

					// L1 = L2 = L3 = 0;
					
					// // L1 = 0;
					
					// // L2 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(U-0.);
					
					// // L3 = 0.25*(V_p-C_p-V)/(Y_out*deltaET)*etdy[i][nyy][k]*(W-0.);
					
				// }
				
				// KAL1 = 0.5/(C_p)*(L4*(V_p+C_p-V)-L5*(V_p-C_p-V));
				
				// KAL2 = L2;
				
				// KAL3 = 0.5/(rho*C_p)*(L4-L5);
				
				// KAL4 = -L3;
				
				// KAL5 = L1+1./rho*(K-1)/K*KAL1;

				
				// KAL1 = KAL1/J[i][nyy][k]+(3*U1_[i][nyy][k]-4*U1[i][nyy][k]+U1q[i][nyy][k])/(2*deltaT);
				// KAL2 = KAL2/J[i][nyy][k]+(3*U2_[i][nyy][k]-4*U2[i][nyy][k]+U2q[i][nyy][k])/(2*deltaT);
				// KAL3 = KAL3/J[i][nyy][k]+(3*U3_[i][nyy][k]-4*U3[i][nyy][k]+U3q[i][nyy][k])/(2*deltaT);
				// KAL4 = KAL4/J[i][nyy][k]+(3*U4_[i][nyy][k]-4*U4[i][nyy][k]+U4q[i][nyy][k])/(2*deltaT);
				// KAL5 = KAL5/J[i][nyy][k]+(3*U5_[i][nyy][k]-4*U5[i][nyy][k]+U5q[i][nyy][k])/(2*deltaT);
				
				// U1_[i][nyy][k] = U1_[i][nyy][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				// U2_[i][nyy][k] = U2_[i][nyy][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				// U3_[i][nyy][k] = U3_[i][nyy][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				// U4_[i][nyy][k] = U4_[i][nyy][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				// U5_[i][nyy][k] = U5_[i][nyy][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
				
				
			}
		}
	}	
		
		
		
}






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

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

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
				H = 0.5*VV+C*C/(K-1);
				
				beta = max(VV/C/C,e);
				
				/* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				temp3 = H-H*K-VV+K*(T+VV);

				d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				d12 = (-4*U*temp2)/temp;
				d13 = (-4*V*temp2)/temp;
				d14 = (-4*W*temp2)/temp;
				d15 = (4*temp2)/temp;

				d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d23 = (-4*U*V*temp2)/temp;
				d24 = (-4*U*W*temp2)/temp;
				d25 = (4*U*temp2)/temp;

				d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				d32 = (-4*U*V*temp2)/temp;
				d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d34 = (-4*V*W*temp2)/temp;
				d35 = (4*V*temp2)/temp;

				d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				d42 = (-4*U*W*temp2)/temp;
				d43 = (-4*V*W*temp2)/temp;
				d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d45 = (4*W*temp2)/temp;

				d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				d52 = (-4*H*U*temp2)/temp;
				d53 = (-4*H*V*temp2)/temp;
				d54 = (-4*H*W*temp2)/temp;
				d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				
				U_p = 0.5*(beta+1)*U;

				C_p = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dx = (P_-P)/deltaXI;

				du_dx = (U_-U)/deltaXI;

				dv_dx = (V_-V)/deltaXI;

				dw_dx = (W_-W)/deltaXI;

				dT_dx = (T_-T)/deltaXI;

				L5 = (U_p-C_p)*(dp_dx-rho*(U_p+C_p-U)*du_dx);

				L4 =  L5*(U_p-C_p-U)/(U_p+C_p-U);
				
				if (U <= 0) {

					L1 = U*(dT_dx+1/rho*dp_dx*(1-K)/K);

					L2 = U*dw_dx;

					L3 = -U*dv_dx;

				}

				else {

					L1 = L2 = L3 = 0;
										
				}

				KAL1 = 0.5/(C_p)*(L4*(U_p+C_p-U)-L5*(U_p-C_p-U));
				
				KAL2 = 0.5/(rho*C_p)*(L4-L5);
				
				KAL3 = -L3;
				
				KAL4 = L2;
				
				KAL5 = L1+1./rho*(K-1)/K*KAL1;
				
				KAL1 = KAL1/J[istart-1][j][k]+(3*U1_[istart-1][j][k]-4*U1[istart-1][j][k]+U1q[istart-1][j][k])/(2*deltaT);
				KAL2 = KAL2/J[istart-1][j][k]+(3*U2_[istart-1][j][k]-4*U2[istart-1][j][k]+U2q[istart-1][j][k])/(2*deltaT);
				KAL3 = KAL3/J[istart-1][j][k]+(3*U3_[istart-1][j][k]-4*U3[istart-1][j][k]+U3q[istart-1][j][k])/(2*deltaT);
				KAL4 = KAL4/J[istart-1][j][k]+(3*U4_[istart-1][j][k]-4*U4[istart-1][j][k]+U4q[istart-1][j][k])/(2*deltaT);
				KAL5 = KAL5/J[istart-1][j][k]+(3*U5_[istart-1][j][k]-4*U5[istart-1][j][k]+U5q[istart-1][j][k])/(2*deltaT);
				
	
				U1_[istart-1][j][k] = U1_[istart-1][j][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				U2_[istart-1][j][k] = U2_[istart-1][j][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				U3_[istart-1][j][k] = U3_[istart-1][j][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				U4_[istart-1][j][k] = U4_[istart-1][j][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				U5_[istart-1][j][k] = U5_[istart-1][j][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
			}
		}


	}


	if (myid == np-1) {

		iend = gend[myid];

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {

				/* outlet */

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
				H = 0.5*VV+C*C/(K-1);

				beta = max(VV/C/C,e);
				
				/* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
				temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
				temp2 = (K-1)*(beta-1)*deltaT*deltaT;
				temp3 = H-H*K-VV+K*(T+VV);

				d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
				d12 = (-4*U*temp2)/temp;
				d13 = (-4*V*temp2)/temp;
				d14 = (-4*W*temp2)/temp;
				d15 = (4*temp2)/temp;

				d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
				d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d23 = (-4*U*V*temp2)/temp;
				d24 = (-4*U*W*temp2)/temp;
				d25 = (4*U*temp2)/temp;

				d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
				d32 = (-4*U*V*temp2)/temp;
				d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d34 = (-4*V*W*temp2)/temp;
				d35 = (4*V*temp2)/temp;

				d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
				d42 = (-4*U*W*temp2)/temp;
				d43 = (-4*V*W*temp2)/temp;
				d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
				d45 = (4*W*temp2)/temp;

				d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
				d52 = (-4*H*U*temp2)/temp;
				d53 = (-4*H*V*temp2)/temp;
				d54 = (-4*H*W*temp2)/temp;
				d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;

				

				U_p = 0.5*(beta+1)*U;

				C_p = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*C);

				dp_dx = (P-_P)/deltaXI;

				du_dx = (U-_U)/deltaXI;

				dv_dx = (V-_V)/deltaXI;

				dw_dx = (W-_W)/deltaXI;

				dT_dx = (T-_T)/deltaXI;

				L4 = (U_p+C_p)*(dp_dx-rho*(U_p-C_p-U)*du_dx);

				L5 = L4*(U_p+C_p-U)/(U_p-C_p-U);
				
				if (U > 0) {

					L1 = U*(dT_dx+1/rho*dp_dx*(1-K)/K);

					L2 = U*dw_dx;

					L3 = -U*dv_dx;
					
				}

				else {

					L1 = L2 = L3 = 0;
					
				}
				
				
				KAL1 = 0.5/(C_p)*(L4*(U_p+C_p-U)-L5*(U_p-C_p-U));
				
				KAL2 = 0.5/(rho*C_p)*(L4-L5);
				
				KAL3 = -L3;
				
				KAL4 = L2;
				
				KAL5 = L1+1./rho*(K-1)/K*KAL1;
				
				
				KAL1 = KAL1/J[iend+1][j][k]+(3*U1_[iend+1][j][k]-4*U1[iend+1][j][k]+U1q[iend+1][j][k])/(2*deltaT);
				KAL2 = KAL2/J[iend+1][j][k]+(3*U2_[iend+1][j][k]-4*U2[iend+1][j][k]+U2q[iend+1][j][k])/(2*deltaT);
				KAL3 = KAL3/J[iend+1][j][k]+(3*U3_[iend+1][j][k]-4*U3[iend+1][j][k]+U3q[iend+1][j][k])/(2*deltaT);
				KAL4 = KAL4/J[iend+1][j][k]+(3*U4_[iend+1][j][k]-4*U4[iend+1][j][k]+U4q[iend+1][j][k])/(2*deltaT);
				KAL5 = KAL5/J[iend+1][j][k]+(3*U5_[iend+1][j][k]-4*U5[iend+1][j][k]+U5q[iend+1][j][k])/(2*deltaT);
				
				U1_[iend+1][j][k] = U1_[iend+1][j][k]-deltaTau*(d11*KAL1+d12*KAL2+d13*KAL3+d14*KAL4+d15*KAL5);
				U2_[iend+1][j][k] = U2_[iend+1][j][k]-deltaTau*(d21*KAL1+d22*KAL2+d23*KAL3+d24*KAL4+d25*KAL5);
				U3_[iend+1][j][k] = U3_[iend+1][j][k]-deltaTau*(d31*KAL1+d32*KAL2+d33*KAL3+d34*KAL4+d35*KAL5);
				U4_[iend+1][j][k] = U4_[iend+1][j][k]-deltaTau*(d41*KAL1+d42*KAL2+d43*KAL3+d44*KAL4+d45*KAL5);
				U5_[iend+1][j][k] = U5_[iend+1][j][k]-deltaTau*(d51*KAL1+d52*KAL2+d53*KAL3+d54*KAL4+d55*KAL5);
				
			}
		}

	}

	

}

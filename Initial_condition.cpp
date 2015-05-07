



#include <mpi.h>
#include <stdlib.h> 
#include <omp.h>
#include <stdio.h>


#include "Resolution.h"

extern int X_np;

void Initial_condition
(
// ======================================================== //
int myid,

int switch_initial,

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
// ======================================================== //
)
{
	
#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"



// ------------------------------------------------------------------------------------- //
// ------------------------------- Initial condition = 0 ------------------------------- //
	double rho, U, V, W, VV, P, T;

//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid];					  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (j = 0; j <= nyyy; j++) {
			for (k = 0; k <= nzzz; k++) {  

				T = 298.0592;
				U = 0;
				V = 0;
				W = 0;    
				VV = U*U+V*V+W*W;
				P = 101300;
				
				rho = P/T/R;
				
				U1[i][j][k] = rho/J[i][j][k];
				U2[i][j][k] = rho*U/J[i][j][k];
				U3[i][j][k] = rho*V/J[i][j][k];
				U4[i][j][k] = rho*W/J[i][j][k];
				U5[i][j][k] = (P/(K-1)+0.5*rho*VV)/J[i][j][k];
				
				U1q[i][j][k] = rho/J[i][j][k];
				U2q[i][j][k] = rho*U/J[i][j][k];
				U3q[i][j][k] = rho*V/J[i][j][k];
				U4q[i][j][k] = rho*W/J[i][j][k];
				U5q[i][j][k] = (P/(K-1)+0.5*rho*VV)/J[i][j][k];

				U1_[i][j][k] = rho/J[i][j][k];
				U2_[i][j][k] = rho*U/J[i][j][k];
				U3_[i][j][k] = rho*V/J[i][j][k];
				U4_[i][j][k] = rho*W/J[i][j][k];
				U5_[i][j][k] = (P/(K-1)+0.5*rho*VV)/J[i][j][k];
				
			}
		}
	}
	
	
// ------------------------------- Initial condition = 0 ------------------------------- //
// ------------------------------------------------------------------------------------- //





// ------------------------------------------------------------------------------------- //
// ------------------------------- Initial condition = 1 ------------------------------- //


	if (switch_initial == 1) {

		
		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];

		MPI_Offset dis;
		
		int Nsize = gcount[myid]*Y_out*Z_out;


		double (*U1out) =  new double[Nsize];
		double (*U2out) =  new double[Nsize];
		double (*U3out) =  new double[Nsize];
		double (*U4out) =  new double[Nsize];
		double (*U5out) =  new double[Nsize];

		double (*U1out_old) =  new double[Nsize];
		double (*U2out_old) =  new double[Nsize];
		double (*U3out_old) =  new double[Nsize];
		double (*U4out_old) =  new double[Nsize];
		double (*U5out_old) =  new double[Nsize];
		
		
		MPI_Offset gs = static_cast<MPI_Offset>(gstart[myid]);
		MPI_Offset ge0 = static_cast<MPI_Offset>(gend0[nproc-1]);
		MPI_Offset one = 1;
		MPI_Offset ge = ge0+one;
		MPI_Offset Yo = static_cast<MPI_Offset>(Y_out);
		MPI_Offset Zo = static_cast<MPI_Offset>(Z_out);
		MPI_Offset si = 8; // sizeof(double) //
		
		
		char data[100];

		sprintf(data,"initial.bin");
		
		MPI_File fh_initial;

		MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_initial ) ; 

		dis = gs * Yo * Zo* si;
		MPI_File_read_at_all(fh_initial, dis, U1out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
		 
		dis = ( (ge*Yo*Zo*1)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U2out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*2)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U3out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
		
		dis = ( (ge*Yo*Zo*3)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U4out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*4)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U5out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


		dis = ( (ge*Yo*Zo*5)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U1out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*6)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U2out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*7)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U3out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*8)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U4out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		dis = ( (ge*Yo*Zo*9)+(gs*Yo*Zo) )*si;
		MPI_File_read_at_all(fh_initial, dis, U5out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
		
		
		MPI_File_close( &fh_initial );

	// =================== //
		istart = 3;        //
		iend = gend[myid]; //
	// =================== //

		  for (i = istart; i <= iend; i++) {

	#pragma omp parallel for private(ii,k)
			  for (j = 2; j <= ny; j++) { 
				  for (k = 2; k <= nz; k++) {

					 ii = (i-istart)*Y_out*Z_out+(j-2)*Z_out+k-2;
					 
					 U1[i][j][k] = U1out[ii]/J[i][j][k];
					 U2[i][j][k] = U2out[ii]/J[i][j][k];
					 U3[i][j][k] = U3out[ii]/J[i][j][k];
					 U4[i][j][k] = U4out[ii]/J[i][j][k];
					 U5[i][j][k] = U5out[ii]/J[i][j][k];

					 U1q[i][j][k] = U1out_old[ii]/J[i][j][k];
					 U2q[i][j][k] = U2out_old[ii]/J[i][j][k];
					 U3q[i][j][k] = U3out_old[ii]/J[i][j][k];
					 U4q[i][j][k] = U4out_old[ii]/J[i][j][k];
					 U5q[i][j][k] = U5out_old[ii]/J[i][j][k];

					 U1_[i][j][k] = U1[i][j][k];
					 U2_[i][j][k] = U2[i][j][k];
					 U3_[i][j][k] = U3[i][j][k];
					 U4_[i][j][k] = U4[i][j][k];
					 U5_[i][j][k] = U5[i][j][k];
					 
				 }
			 }
			
		 }
		 
		 
		 delete [] U1out;
		 delete [] U2out;
		 delete [] U3out;
		 delete [] U4out;
		 delete [] U5out;

		 delete [] U1out_old;
		 delete [] U2out_old;
		 delete [] U3out_old;
		 delete [] U4out_old;
		 delete [] U5out_old;
			
	}  // ---- if (switch_initial == 1) ---- //
		

// ------------------------------- Initial condition = 1 ------------------------------- //
// ------------------------------------------------------------------------------------- //



	
}
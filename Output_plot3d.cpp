



#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

void Output_plot3d
	(
	// ============================================================================ //
	int step,

	int myid,
	
	int switch_output,

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

	double (*X_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
	double (*Y_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
	double (*Z_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

	double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "MPI_prm.h"
#include "Mpi_division.h"
#include "prm.h"


	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Offset disp;

	double temp = 1.0;    // mach, alpha, reyn, time //
	
	int last_count, next_count;

	int Nblock = 1;
	
	int nx_out = X_out;
	int ny_out = Y_out-ny_abs;
	int nz_out = Z_out;

	double tempX,tempY,tempZ;
	
	MPI_Offset dis;

	MPI_Offset Nsize;
	
	
	MPI_Offset  x_gdisp[np];


	istart = 0;
	
	
	
	x_gdisp[0] = 0; 
	
	for (i = 0; i < np; i++) { 

		if (i < (np-1)) {
			
			x_gdisp[i+1] = x_gdisp[i]+(MPI_Offset)gcount[i];  // ---- how many grids in X-direction ---- //
			
		}
		
	}
	
	
		// printf("%d\t%d\t%d\t%d\t%d\n",myid,gstart[myid],gend[myid],gcount[myid],x_gdisp[myid]);

	
	
	double (*Yout)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];

	double gamma1 = 2.8;
	double gamma2 = 2.5;

	for (i = 0; i < nx_out; i++) {
		for (j = 0; j < ny_out; j++) { 
			for (k = 0; k < nz_out; k++) { 

				if (j < Y_out-ny_abs)
					Yout[i][j][k] = (0.5*high)*(1-1./tanh(gamma1)*tanh(gamma1*(1-2*(j+0.5)*deltaET)));
				else {

					Yout[i][j][k] = (0.5*high)*(1-1./tanh(gamma2)*tanh(gamma2*(1-2*(j-(ny_out-ny_abs)+0.5)*deltaET)))+high;

				}

			}
		}
	}  


	
	
	
	char LESdata[100];

	if ( switch_output == 1 ) {

		sprintf(LESdata,"P3D_LES""%0.5d"".x",1);
		MPI_File fh0;


		MPI_File_open( MPI_COMM_WORLD, LESdata,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh0 ) ; 

		if (myid == np-1) {

			FILE *fptr_xyz;

			fptr_xyz = fopen(LESdata,"wb");

			fwrite(&Nblock, sizeof(int), 1,fptr_xyz);

			fwrite(&nx_out, sizeof(int), 1,fptr_xyz);
			fwrite(&ny_out, sizeof(int), 1,fptr_xyz);
			fwrite(&nz_out, sizeof(int), 1,fptr_xyz);
			
			fclose(fptr_xyz);

		}

		
		float (*Grid_X) = new float[gcount[myid]*ny_out*nz_out];
		float (*Grid_Y) = new float[gcount[myid]*ny_out*nz_out];
		float (*Grid_Z) = new float[gcount[myid]*ny_out*nz_out];

//// ============================== ////
		 istart = 3;		        ////	
//// ============================== ////
		iend = gend[myid];	        ////
//// ============================== ////

		
		icount = -1;
	
		for (k = 2; k <= nz; k++) { 
			for (j = 2; j <= ny; j++) { 
				for (i = istart; i <= iend; i++) {
				
					icount = icount+1;

					Grid_X[icount] = deltaXI*(i-istart+0.5)+x_gdisp[myid]*deltaXI;
					
				}
			}
		}

		
		icount = -1;
	
		for (k = 2; k <= nz; k++) { 
			for (j = 2; j <= ny; j++) { 
				for (i = istart; i <= iend; i++) {

					icount = icount+1;

					Grid_Y[icount] = Yout[i-istart][j-2][k-2];

				}
			}
		}

		
		icount = -1;
	
		for (k = 2; k <= nz; k++) { 
			for (j = 2; j <= ny; j++) { 
				for (i = istart; i <= iend; i++) {

					icount = icount+1;

					Grid_Z[icount] = deltaZT*(k-1.5);

				}
			}
		}

		
		disp = (MPI_Offset)sizeof(int)*4 + x_gdisp[myid]*(MPI_Offset)ny_out*(MPI_Offset)nz_out*(MPI_Offset)sizeof(float);

		MPI_File_write_at_all(fh0, disp, Grid_X, gcount[myid]*ny_out*nz_out, MPI_FLOAT, MPI_STATUS_IGNORE);
		
		
		disp = (MPI_Offset)sizeof(int)*4 + (MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*(MPI_Offset)sizeof(float)+\
		       x_gdisp[myid]*(MPI_Offset)ny_out*(MPI_Offset)nz_out*(MPI_Offset)sizeof(float);

		MPI_File_write_at_all(fh0, disp, Grid_Y, gcount[myid]*ny_out*nz_out, MPI_FLOAT, MPI_STATUS_IGNORE);
		
		
		disp = (MPI_Offset)sizeof(int)*4 + 2*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*(MPI_Offset)sizeof(float)+\
		       x_gdisp[myid]*(MPI_Offset)ny_out*(MPI_Offset)nz_out*(MPI_Offset)sizeof(float);

		MPI_File_write_at_all(fh0, disp, Grid_Z, gcount[myid]*ny_out*nz_out, MPI_FLOAT, MPI_STATUS_IGNORE);

		
		delete []Grid_X;
		delete []Grid_Y;
		delete []Grid_Z;

		MPI_File_close( &fh0 );

	}

	

	delete [] Yout;
	
	
	
	
	
	
	float (*Solution) = new float[gcount[myid]*ny_out*nz_out*5];

	
	sprintf(LESdata,"qP3D_LES""%0.5d"".q",step);
	MPI_File fh1;
	
	MPI_File_open( MPI_COMM_WORLD, LESdata,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1 ) ; 

	if (myid == 0) {

		FILE *fptr_solution;

		fptr_solution = fopen(LESdata,"wb");

		fwrite(&Nblock, sizeof(int), 1,fptr_solution);

		fwrite(&nx_out, sizeof(int), 1,fptr_solution);
		fwrite(&ny_out, sizeof(int), 1,fptr_solution);
		fwrite(&nz_out, sizeof(int), 1,fptr_solution);
			
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);

		fclose(fptr_solution);

	}

	icount = -1;

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U1[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U2[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U3[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U4[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U5[i][j][k]*J[i][j][k];

			}
		}
	}

	disp = (MPI_Offset)sizeof(int)*8 + x_gdisp[myid]*(MPI_Offset)ny_out*(MPI_Offset)nz_out*5*(MPI_Offset)sizeof(float);

	MPI_File_write_at_all(fh1, disp, Solution, gcount[myid]*ny_out*nz_out*5, MPI_FLOAT, MPI_STATUS_IGNORE);

	MPI_File_close( &fh1 );
	


	
// -------------------------------------- Previsous time step -------------------------------------- //

	
	sprintf(LESdata,"QqP3D_LES""%0.5d"".q",step);
	MPI_File fh2;
	
	MPI_File_open( MPI_COMM_WORLD, LESdata,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh2 ) ; 

	if (myid == 0) {

		FILE *fptr_solution;

		fptr_solution = fopen(LESdata,"wb");

		fwrite(&Nblock, sizeof(int), 1,fptr_solution);

		fwrite(&nx_out, sizeof(int), 1,fptr_solution);
		fwrite(&ny_out, sizeof(int), 1,fptr_solution);
		fwrite(&nz_out, sizeof(int), 1,fptr_solution);
			
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);
		fwrite(&temp, sizeof(double), 1,fptr_solution);

		fclose(fptr_solution);

	}

	icount = -1;

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U1q[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U2q[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U3q[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U4q[i][j][k]*J[i][j][k];

			}
		}
	}

	for (k = 2; k <= nz; k++) { 
		for (j = 2; j <= ny; j++) { 
			for (i = istart; i <= iend; i++) {
			
				icount = icount + 1;
				Solution[icount] = U5q[i][j][k]*J[i][j][k];

			}
		}
	}

	
	disp = (MPI_Offset)sizeof(int)*8 + x_gdisp[myid]*(MPI_Offset)ny_out*(MPI_Offset)nz_out*5*(MPI_Offset)sizeof(float);

	
	MPI_File_write_at_all(fh2, disp, Solution, gcount[myid]*ny_out*nz_out*5, MPI_FLOAT, MPI_STATUS_IGNORE);

	MPI_File_close( &fh2 );

	
	delete [] Solution;
	


}




#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

void Output
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


	double temp = 1.0;    // mach, alpha, reyn, time //

	int Nblock = 1;
	int nx_out = X_out;
	int ny_out = Y_out-ny_abs;
	int nz_out = Z_out;

	double tempX,tempY,tempZ;

	MPI_Offset dis;

	//unsigned long long dis;

	MPI_Offset Nsize;


	if ( switch_output == 1 ) {


		//double (*Xout)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
		double (*Yout)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
		//double (*Zout)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];


		if (myid == 0) {

			double gamma1 = 2.8;
			double gamma2 = 2.5;

			//#pragma omp parallel for private(j,k)
			for (i = 0; i < nx_out; i++) {
				for (j = 0; j < ny_out; j++) { 
					for (k = 0; k < nz_out; k++) { 

						//Xout[i][j][k] = deltaXI*(i+0.5);

						if (j < Y_out-ny_abs)
							Yout[i][j][k] = (0.5*high)*(1-1./tanh(gamma1)*tanh(gamma1*(1-2*(j+0.5)*deltaET)));
						else {

							Yout[i][j][k] = (0.5*high)*(1-1./tanh(gamma2)*tanh(gamma2*(1-2*(j-(ny_out-ny_abs)+0.5)*deltaET)))+high;

							//if(i==0 && k == 0) printf("%f\n",Yout[i][j][k]);

						}


						//Zout[i][j][k] = deltaZT*(k+0.5);


					}
				}
			}  



			char LESdata[100];
			FILE *fptr;
			sprintf(LESdata,"P3D_LES""%0.5d"".x",1);
			fptr = fopen(LESdata,"wb");

			fwrite(&Nblock, sizeof(int), 1,fptr);

			fwrite(&nx_out, sizeof(int), 1,fptr);
			fwrite(&ny_out, sizeof(int), 1,fptr);
			fwrite(&nz_out, sizeof(int), 1,fptr);

			for (k = 0; k < nz_out; k++) { 
				for (j = 0; j < ny_out; j++) { 
					for (i = 0; i < nx_out; i++) {

						//fwrite(&Xout[i][j][k],sizeof(double),1,fptr);

						tempX = deltaXI*(i+0.5);

						fwrite(&tempX,sizeof(double),1,fptr);

					}
				}
			}

			for (k = 0; k < nz_out; k++) { 
				for (j = 0; j < ny_out; j++) { 
					for (i = 0; i < nx_out; i++) {

						fwrite(&Yout[i][j][k],sizeof(double),1,fptr);

					}
				}
			}

			for (k = 0; k < nz_out; k++) { 
				for (j = 0; j < ny_out; j++) { 
					for (i = 0; i < nx_out; i++) {

						//fwrite(&Zout[i][j][k],sizeof(double),1,fptr);

						tempZ = deltaZT*(k+0.5);

						fwrite(&tempZ,sizeof(double),1,fptr);

					}
				}
			}

			fclose(fptr);

		}    // ---- if (myid == 0) ---- //


		//delete [] Xout;
		delete [] Yout;
		//delete [] Zout;
		
	}



	double (*U1out) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U2out) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U3out) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U4out) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U5out) =  new double[gcount[myid]*Y_out*Z_out];

	double (*U1out_old) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U2out_old) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U3out_old) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U4out_old) =  new double[gcount[myid]*Y_out*Z_out];
	double (*U5out_old) =  new double[gcount[myid]*Y_out*Z_out];

	// =================== //
	istart = 3;        //
	iend = gend[myid]; //
	// =================== //

	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(ii,k)
		for (j = 2; j <= ny; j++) { 
			for (k = 2; k <= nz; k++) {

				ii = (i-istart)*Y_out*Z_out+(j-2)*Z_out+k-2;

				U1out[ii] = U1[i][j][k]*J[i][j][k]; 
				U2out[ii] = U2[i][j][k]*J[i][j][k]; 
				U3out[ii] = U3[i][j][k]*J[i][j][k];  
				U4out[ii] = U4[i][j][k]*J[i][j][k]; 
				U5out[ii] = U5[i][j][k]*J[i][j][k]; 

				U1out_old[ii] = U1q[i][j][k]*J[i][j][k]; 
				U2out_old[ii] = U2q[i][j][k]*J[i][j][k];  
				U3out_old[ii] = U3q[i][j][k]*J[i][j][k]; 
				U4out_old[ii] = U4q[i][j][k]*J[i][j][k];   
				U5out_old[ii] = U5q[i][j][k]*J[i][j][k]; 


				//if(myid == 0) printf("%f\t%f\t%f\t%f\n",U1out[ii],U2out[ii],U1out_old[ii],U2out_old[ii]);


			}
		}
	}


	Nsize = gcount[myid]*Y_out*Z_out;

	char LESdata[100];
	sprintf(LESdata,"LES""%0.5d"".bin",step);

	MPI_File fh0;

	MPI_Offset gs = static_cast<MPI_Offset>(gstart[myid]);
	MPI_Offset ge0 = static_cast<MPI_Offset>(gend0[nproc-1]);
	MPI_Offset one = 1;
	MPI_Offset ge = ge0+one;
	MPI_Offset Xo = static_cast<MPI_Offset>(X_out);
	MPI_Offset Yo = static_cast<MPI_Offset>(Y_out);
	MPI_Offset Zo = static_cast<MPI_Offset>(Z_out);
	MPI_Offset si = 8; // sizeof(double) //


	MPI_File_open( MPI_COMM_WORLD, LESdata,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh0 ) ; 


	dis = gs * Yo * Zo* si;

	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U1out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*1)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U2out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*2)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U3out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*3)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U4out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*4)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U5out, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*5)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U1out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*6)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U2out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*7)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U3out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*8)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U4out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	dis = ( (ge*Yo*Zo*9)+(gs*Yo*Zo) )*si;
	MPI_File_set_view(fh0, dis, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fh0, U5out_old, Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


	MPI_File_close( &fh0 );


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

	// double (*U1s) =  new double[gcount[myid]*ny_out*nz_out];
	// double (*U2s) =  new double[gcount[myid]*ny_out*nz_out];
	// double (*U3s) =  new double[gcount[myid]*ny_out*nz_out];
	// double (*U4s) =  new double[gcount[myid]*ny_out*nz_out];
	// double (*U5s) =  new double[gcount[myid]*ny_out*nz_out];


	// // =================== //
	// istart = 3;        //
	// iend = gend[myid]; //
	// // =================== //

// #pragma omp parallel for private(ii,j,i)

	
	// for (k = 2; k <= nz; k++) {
		// for (j = 2; j <= ny_out+1; j++) { 
			// for (i = istart; i <= iend; i++) {

				// ii = (k-2)*ny_out*(iend-istart+1)+(j-2)*(iend-istart+1)+(i-istart);


				// U1s[ii] = U1[i][j][k]*J[i][j][k]; 
				// U2s[ii] = U2[i][j][k]*J[i][j][k]; 
				// U3s[ii] = U3[i][j][k]*J[i][j][k];  
				// U4s[ii] = U4[i][j][k]*J[i][j][k]; 
				// U5s[ii] = U5[i][j][k]*J[i][j][k]; 

				// //if(myid == 0 && k == 5 && i == 50) printf("%d\t%f\t%f\t%f\n",j,U1s[ii],U1[i][j][k],J[i][j][k]);

				// //if(myid == 0) printf("%d\t%d\t%d\t%d\t%f\n",i,j,k,ii,U1s[ii]);

				// //printf("%d\t%d\t%d\t%d\t%f\n",ii,j,k,myid,U1s[ii]);

			// }
		// }
	// }





	// Nsize = gcount[myid];

	
	// MPI_Offset ii_index;
	
	

	// sprintf(LESdata,"qP3D_LES""%0.5d"".q",step);

	// MPI_File fh_plot3d;

	// MPI_File_open( MPI_COMM_WORLD, LESdata,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh_plot3d ) ; 


		// if (myid == 0) {

			// char LESdata[100];
			// FILE *fptr;
			// sprintf(LESdata,"qP3D_LES""%0.5d"".q",step);
			// fptr = fopen(LESdata,"wb");

			// fwrite(&Nblock, sizeof(int), 1,fptr);

			// fwrite(&nx_out, sizeof(int), 1,fptr);
			// fwrite(&ny_out, sizeof(int), 1,fptr);
			// fwrite(&nz_out, sizeof(int), 1,fptr);

			// fwrite(&temp, sizeof(double), 1,fptr);
			// fwrite(&temp, sizeof(double), 1,fptr);
			// fwrite(&temp, sizeof(double), 1,fptr);
			// fwrite(&temp, sizeof(double), 1,fptr);

		// }    // ---- total bytes 48 ---- //

	

	// for (k = 0; k < nz_out; k++) { 
		// for (j = 0; j < ny_out; j++) { 

				// ii_index = k*(ny_out*gcount[myid])+j*gcount[myid];

				// dis = 48+gs*sizeof(double)+(k*ny_out+j)*nx_out*sizeof(double);

				// MPI_File_write_at_all(fh_plot3d, dis, &U1s[ii_index], Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


		// }
	// }


	// for (k = 0; k < nz_out; k++) { 
		// for (j = 0; j < ny_out; j++) { 

				// ii_index = k*(ny_out*gcount[myid])+j*gcount[myid];

				// dis = 48+gs*sizeof(double)+(k*ny_out+j)*nx_out*sizeof(double)+ 1 * Xo * Yo * Zo* si;

				// MPI_File_write_at_all(fh_plot3d, dis, &U2s[ii_index], Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


		// }
	// }



	// for (k = 0; k < nz_out; k++) { 
		// for (j = 0; j < ny_out; j++) { 

				// ii_index = k*(ny_out*gcount[myid])+j*gcount[myid];

				// dis = 48+gs*sizeof(double)+(k*ny_out+j)*nx_out*sizeof(double)+ 2 * Xo * Yo * Zo* si;;

				// MPI_File_write_at_all(fh_plot3d, dis, &U3s[ii_index], Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


		// }
	// }

	// for (k = 0; k < nz_out; k++) { 
		// for (j = 0; j < ny_out; j++) { 

				// ii_index = k*(ny_out*gcount[myid])+j*gcount[myid];

				// dis = 48+gs*sizeof(double)+(k*ny_out+j)*nx_out*sizeof(double)+ 3 * Xo * Yo * Zo* si;;

				// MPI_File_write_at_all(fh_plot3d, dis, &U4s[ii_index], Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

		// }
	// }

	// for (k = 0; k < nz_out; k++) { 
		// for (j = 0; j < ny_out; j++) { 

				// ii_index = k*(ny_out*gcount[myid])+j*gcount[myid];

				// dis = 48+gs*sizeof(double)+(k*ny_out+j)*nx_out*sizeof(double)+ 4 * Xo * Yo * Zo* si;;

				// MPI_File_write_at_all(fh_plot3d, dis, &U5s[ii_index], Nsize, MPI_DOUBLE, MPI_STATUS_IGNORE);


		// }
	// }



    // MPI_File_close( &fh_plot3d );



	// delete [] U1s;
	// delete [] U2s;
	// delete [] U3s;
	// delete [] U4s;
	// delete [] U5s;



}
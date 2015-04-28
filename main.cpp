




#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "main.h"
#include "MPI_prm.h"

int main(int argc, char **argv)
{   

//==== MPI start ====//

	int myid;

// ============================================ //
	MPI_Status istat[8];					    //
												//
	MPI_Comm comm;								//
												//
	MPI_Init (&argc, &argv);					//
												//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);		//
												//
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);		//
												//
	comm=MPI_COMM_WORLD;						//
// ============================================ //
	
#include "ijk.h"
	
#include "Resolution.h"

	int statistic_step = 10000;
	int iteration_end_step = 100;
	int output_step = 50;
	int count = 1001;
	int step;
	
	double e = 0.00002;
	double deltaT = 0.0008;
	double Ep = 1.;
	double Roe_criterion = 0.007;
	double E_high = 1.;
	double E_low = 0.1;

	double heat_flux = 104;

	double obs1 = 0.11;
	double obs2 = 0.82;
	double obs3 = 1.69;
	double obs4 = 2.63;
	double obs5 = 3.84;
	
	nproc = np;

	double er1, er2, er3, er4, er5;
	
	double E1 = 0;
	double E2 = 0;
	double E3 = 0;
	double E4 = 0;
	double E5 = 0;

#include "Mpi_division.h"

	
int X_np = gcount[myid]+6;    /**** How many cells in X-direction for each CPU ****/

#include "Array.h"    /**** allocate the memory ****/


	Grid(myid, X_point, Y_point, Z_point);


	Hyperplane(myid, Nijk, ijk);


	Jacobian(myid, J,
		     xidx,xidy,xidz,etdx,etdy,etdz,ztdx,ztdy,ztdz);
	

	Jacobian_u(myid, J_u,
		       xidx_u,xidy_u,xidz_u,etdx_u,etdy_u,etdz_u,ztdx_u,ztdy_u,ztdz_u);


	Jacobian_v(myid, J_v,
		       xidx_v,xidy_v,xidz_v,etdx_v,etdy_v,etdz_v,ztdx_v,ztdy_v,ztdz_v,
		       J,
		       xidx,xidy,xidz,etdx,etdy,etdz,ztdx,ztdy,ztdz);

	Jacobian_w(myid, J_w,
		       xidx_w,xidy_w,xidz_w,etdx_w,etdy_w,etdz_w,ztdx_w,ztdy_w,ztdz_w);

	
	
	Initial_condition(myid,
					  U1_,U2_,U3_,U4_,U5_,
					  U1, U2, U3, U4, U5,
					  U1q,U2q,U3q,U4q,U5q,
					  J);
					  
	

// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
	for (i = istart; i <= iend; i++) {
#pragma omp parallel for private(k)
			for (j = 2; j < nyy; j++) {
				for (k = 2; k < nzz; k++) {     

					U1_[i][j][k] = U1[i][j][k];
					U2_[i][j][k] = U2[i][j][k];
					U3_[i][j][k] = U3[i][j][k];
					U4_[i][j][k] = U4[i][j][k];
					U5_[i][j][k] = U5[i][j][k];
					
				}
			}
#pragma omp barrier
		}

//// ============================================ ////

	





// =============================================== //
	for (step = 1 ; step <= count; step++) {       //
// =============================================== //


// ============================================================================ //
 		for (int iteration = 1; iteration < 50000; iteration++) {               //
// ============================================================================ //		




////****  data transformation among nodes ****//// 
//// =========================================================================================================================================================== ////
			
			
			istart=3;
			iend = gend[myid];

			l_nbr = myid - 1;
			r_nbr = myid + 1;
			if(myid == 0) l_nbr=MPI_PROC_NULL;
			if(myid == nproc-1) r_nbr=MPI_PROC_NULL;

			icount = 3*Y_m*Z_m;

			
			itag=110;
			MPI_Sendrecv((void *)&U1_[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&U1_[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
			itag=120;
			MPI_Sendrecv((void *)&U2_[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&U2_[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
			itag=130;
			MPI_Sendrecv((void *)&U3_[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&U3_[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
			itag=140;
			MPI_Sendrecv((void *)&U4_[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&U4_[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
			itag=150;
			MPI_Sendrecv((void *)&U5_[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&U5_[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
			
			itag=160;
			MPI_Sendrecv((void *)&U1_[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&U1_[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
			itag=170;
			MPI_Sendrecv((void *)&U2_[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&U2_[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
			itag=180;
			MPI_Sendrecv((void *)&U3_[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&U3_[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
			itag=190;
			MPI_Sendrecv((void *)&U4_[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&U4_[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
			itag=200;
			MPI_Sendrecv((void *)&U5_[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&U5_[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
			
			MPI_Barrier(MPI_COMM_WORLD);
						
//// =========================================================================================================================================================== ////
////****  data transformation among nodes-end ****//// 
			

			Y_boundary_condition(myid,deltaT,e,heat_flux,U1_,U2_,U3_,U4_,U5_,U1,U2,U3,U4,U5,U1q,U2q,U3q,U4q,U5q,etdy,J);
			
			X_boundary_condition(myid,deltaT,e,U1_,U2_,U3_,U4_,U5_,U1,U2,U3,U4,U5,U1q,U2q,U3q,U4q,U5q,J);

			Z_boundary_condition(myid,U1_,U2_,U3_,U4_,U5_,J);

			
			Flux_X(myid, Ep, e,
				   U1_,U2_,U3_,U4_,U5_,MR1,MR2,MR3,MR4,MR5,ML1,ML2,ML3,ML4,ML5,
				   J,J_u,
				   xidx_u,xidy_u,xidz_u,etdx_u,etdy_u,etdz_u,ztdx_u,ztdy_u,ztdz_u,
			       inFx1,inFx2,inFx3,inFx4,inFx5, EpX);  

			Flux_Y(myid, Ep, e,
				   U1_,U2_,U3_,U4_,U5_,MR1,MR2,MR3,MR4,MR5,ML1,ML2,ML3,ML4,ML5,
				   J,J_v,
				   xidx_v,xidy_v,xidz_v,etdx_v,etdy_v,etdz_v,ztdx_v,ztdy_v,ztdz_v,
			       inFy1,inFy2,inFy3,inFy4,inFy5, EpY); 

			Flux_Z(myid, Ep, e,
				   U1_,U2_,U3_,U4_,U5_,MR1,MR2,MR3,MR4,MR5,ML1,ML2,ML3,ML4,ML5,
				   J,J_w,
				   xidx_w,xidy_w,xidz_w,etdx_w,etdy_w,etdz_w,ztdx_w,ztdy_w,ztdz_w,
			       inFz1,inFz2,inFz3,inFz4,inFz5, EpZ);  


			Viscous_terms(myid,
						  U1_,U2_,U3_,U4_,U5_,vF2,vF3,vF4,vF5,
						  J,
						  xidx,xidy,xidz,etdx,etdy,etdz,ztdx,ztdy,ztdz,
						  J_u,
						  xidx_u,xidy_u,xidz_u,etdx_u,etdy_u,etdz_u,ztdx_u,ztdy_u,ztdz_u,
						  J_v,
						  xidx_v,xidy_v,xidz_v,etdx_v,etdy_v,etdz_v,ztdx_v,ztdy_v,ztdz_v,
						  J_w,
						  xidx_w,xidy_w,xidz_w,etdx_w,etdy_w,etdz_w,ztdx_w,ztdy_w,ztdz_w,
						  LR1,LR2,LR3,LR4,LR5,LL1,LL2,LL3,LL4,LL5,
						  MR1,MR2,MR3,MR4,MR5,ML1,ML2,ML3,ML4,ML5,
						  NR1,NR2,NR3,NR4,NR5,NL1,NL2,NL3,NL4,NL5);

			
			Dual_time_stepping(myid, deltaT, e,
							   &er1,&er2,&er3,&er4,&er5,
							   Nijk,ijk,
							   U1_,U2_,U3_,U4_,U5_,
							   U1 ,U2 ,U3 ,U4 ,U5 ,
							   U1q,U2q,U3q,U4q,U5q,
							   J,
						       xidx,xidy,xidz,etdx,etdy,etdz,ztdx,ztdy,ztdz,
							   inFx1,inFx2,inFx3,inFx4,inFx5,
							   inFy1,inFy2,inFy3,inFy4,inFy5,
							   inFz1,inFz2,inFz3,inFz4,inFz5,
							   vF2,vF3,vF4,vF5,
							   MR1,MR2,MR3,MR4,MR5,
							   NR1,NR2,NR3,NR4,NR5,
							   NL1,NL2,NL3,NL4,NL5,
							   Residual1,Residual2,Residual3,Residual4,Residual5,
							   etdy_v,J_v
							   );
			
			

			/*
			DPLUSGS(myid, deltaT,
				    &er1,&er2,&er3,&er4,&er5,
					U1_,U2_,U3_,U4_,U5_,
					U1 ,U2 ,U3 ,U4 ,U5 ,
					U1q,U2q,U3q,U4q,U5q,
					J,
					xidx,xidy,xidz,etdx,etdy,etdz,ztdx,ztdy,ztdz,
					inFx1,inFx2,inFx3,inFx4,inFx5,
					inFy1,inFy2,inFy3,inFy4,inFy5,
					inFz1,inFz2,inFz3,inFz4,inFz5,
					vF2,vF3,vF4,vF5,
					MR1,MR2,MR3,MR4,MR5,
					NR1,NR2,NR3,NR4,NR5,
					NL1,NL2,NL3,NL4,NL5,
					Residual1,Residual2,Residual3,Residual4,Residual5);
			*/
			
			
	E1 = max(E1, er1);
	E2 = max(E2, er2);
	E3 = max(E3, er3);
	E4 = max(E4, er4);
	E5 = max(E5, er5);
	
	er1 = er1/E1;
	er2 = er2/E2;
	er3 = er3/E3;
	er4 = er4/E4;
	er5 = er5/E5;
	
	// if (myid == 0) 
		// printf("%d\t%f\t%f\t%f\t%f\t%f\n",iteration,er1,er2,er3,er4,er5);

		
// ================================================================================================ //
			if ((er1<0.001 & er2<0.001 & er3<0.001 & er4<0.001 & er5<0.001) | iteration == iteration_end_step) {    //
// ================================================================================================ //

				Driven_force(myid, U2,J);

				
				//// ============================================================================== ////
				if (myid == 0) {
		
					printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",step,iteration,er1,er2,er3,er4,er5);

					//printf("iteration = %d\n",iteration);


					// ====================================================================== //
					// char LESdata[100];
					// FILE *fptr;
					// sprintf(LESdata,"output.dat");
					// fptr = fopen(LESdata,"a");

					// fprintf(fptr,"%f\t%f\t%f\t%f\t%d\n\n",er2,er3,er4,er5,iteration);

					// fclose(fptr);

					// ====================================================================== //

				}
				//// ============================================================================= ////

			
			Statistic(myid,
					  step,iteration, statistic_step,
					  obs1, obs2, obs3, obs4, obs5,
					  e,
				      U1_,U2_,U3_,U4_,U5_,
				      J,
					  J_v,
					  xidx_v,xidy_v,xidz_v,etdx_v,etdy_v,etdz_v,ztdx_v,ztdy_v,ztdz_v,
					  MR1,MR2,MR3,MR4,MR5,ML1,ML2,ML3,ML4,ML5,EpY);
					  
			





			
//// ============================================ ////
		 istart = 3;		     	              ////	
//// ============================================ ////
		iend = gend[myid];				    	  ////
//// ============================================ ////

#pragma omp parallel for private(j,k)
		for (i = istart-1; i <= iend+1; i++) {
			for (j = 0; j <= nyyy; j++) { 
				for (k = 0; k <= nzzz; k++) {
			//for (j = 2; j < nyy; j++) { 
				//for (k = 2; k < nzz; k++) {


					U1q[i][j][k] = U1[i][j][k]; 
					U2q[i][j][k] = U2[i][j][k]; 
					U3q[i][j][k] = U3[i][j][k]; 
					U4q[i][j][k] = U4[i][j][k]; 
					U5q[i][j][k] = U5[i][j][k];

					U1[i][j][k] = U1_[i][j][k]; 
					U2[i][j][k] = U2_[i][j][k]; 
					U3[i][j][k] = U3_[i][j][k]; 
					U4[i][j][k] = U4_[i][j][k]; 
					U5[i][j][k] = U5_[i][j][k];

				}
			}
		}

		break;  /**** jump out the iteration ****/

			}
			else {

				//driven_force1 ();


// ================================================================================================ //
			}	/**** if-end ****/																    //
// ================================================================================================ //




// ============================================================================ //
	}    /**** iteration-end ****/   											//
// ============================================================================ //

//// ================================================================================================== ////
	if ((step%output_step) == 0) {

		Output(step,myid,U1,U2,U3,U4,U5,U1q,U2q,U3q,U4q,U5q,X_point, Y_point, Z_point, J);

	}

//// ================================================================================================== ////


// =============================================== //
	}    /**** step-end ****/                      // 
// =============================================== //


#include "Delete_Array.h"

// ============================================ //
MPI_Finalize();    /**** MPI-end ****/          //
// ============================================ //

}
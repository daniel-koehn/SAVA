/*------------------------------------------------------------------------
 *   This is function initproc.
 *   Dividing the FD grid into domains and assigning the
 *   PEs to these domains.
 *
 *   O. Hellwig
 *  ----------------------------------------------------------------------*/
 

#include "fd.h"

void initproc(void){

	/* extern variables */
	extern int	NXG[3];
	extern int	NX[3], POS[3], INDEX[6];
	extern int	NP, NPROC, NPROCX[3], MYID;
	extern int	PERIODIC[3];
	extern FILE	*FP;
	extern MPI_Comm	COMM_CART;
	
	/* local variables */
	const int	opt = 1;	/* allow for optimization (1...yes, 0...no) */
	int	   	MYID_CART;


	if (NPROC != NP)
		error("Number of processors specified in the parameter file \n and at command line (NP) differ !");

	/* create communicator for Cartesian topology */
	MPI_Cart_create(MPI_COMM_WORLD,3,NPROCX,PERIODIC,opt,&COMM_CART);

	/* new proc-ID in Cartesian grid */
	MPI_Comm_rank(COMM_CART,&MYID_CART);

	/* processor location in the 3D logical processor array */
	/* POS[0] ... x coordinate */
	/* POS[1] ... y coordinate */
	/* POS[2] ... z coordinate */
	MPI_Cart_coords(COMM_CART,MYID_CART,3,POS);

	/* neighbouring processes */
	/* INDEX[0] ... left  */
	/* INDEX[1] ... right */
	/* INDEX[2] ... back */
	/* INDEX[3] ... front */
	/* INDEX[4] ... upper */
	/* INDEX[5] ... lower */
	MPI_Cart_shift(COMM_CART,0,1,&INDEX[0],&INDEX[1]);
	MPI_Cart_shift(COMM_CART,1,1,&INDEX[2],&INDEX[3]);
	MPI_Cart_shift(COMM_CART,2,1,&INDEX[4],&INDEX[5]);

	/* length of the subarray on this processor */
	NX[0] = NXG[0]/NPROCX[0];
	NX[1] = NXG[1]/NPROCX[1];
	NX[2] = NXG[2]/NPROCX[2];

	if (NXG[0]%NPROCX[0])
		error(" NXG[0]%NPROCX[0] must be zero !");
	if (NXG[1]%NPROCX[1])
		error(" NXG[1]%NPROCX[1] must be zero !");
	if (NXG[2]%NPROCX[2])
		error(" NXG[2]%NPROCX[2] must be zero !");

	fprintf(FP,"\n **Message from initprocs (printed by PE %d):\n",MYID);
	fprintf(FP," Size of subarrays in gridpoints:\n");
	fprintf(FP," NX            = %d\n",NX[0]);
	fprintf(FP," NY            = %d\n",NX[1]);
	fprintf(FP," NZ (vertical) = %d\n",NX[2]);

	if (!(MYID)){
		fprintf(stdout,"\n");
		fprintf(stdout," Processor locations in the 3D logical processor array\n");
		fprintf(stdout," MYID \t POS(0): left,right \t POS(1): back,front \t POS(2): top,bottom\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	fprintf(stdout," %d \t\t %d: %d,%d \t\t %d: %d,%d \t\t %d: %d,%d \n",
		MYID,POS[0],INDEX[0],INDEX[1], POS[1],INDEX[2],INDEX[3], POS[2],INDEX[4],INDEX[5]);

	MPI_Barrier(MPI_COMM_WORLD);
}
